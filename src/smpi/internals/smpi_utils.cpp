/* Copyright (c) 2016-2021. The SimGrid Team.
 * All rights reserved.                                                     */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "smpi_utils.hpp"

#include "src/surf/xml/platf_private.hpp"
#include "xbt/log.h"
#include "xbt/parse_units.hpp"
#include "xbt/sysdep.h"
#include "xbt/file.hpp"
#include <boost/tokenizer.hpp>
#include "smpi_config.hpp"
#include "src/simix/smx_private.hpp"
#include <algorithm>
#include "private.hpp"

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(smpi_utils, smpi, "Logging specific to SMPI (utils)");

extern std::string surf_parsed_filename;
extern int surf_parse_lineno;

namespace simgrid {
namespace smpi {
namespace utils {

double total_benched_time=0;
unsigned long total_malloc_size=0;
unsigned long total_shared_size=0;
unsigned int total_shared_calls=0;
struct alloc_metadata_t {
  size_t size          = 0;
  unsigned int numcall = 0;
  int line             = 0;
  std::string file;
};

struct current_buffer_metadata_t {
  alloc_metadata_t alloc;
  std::string name;
};

alloc_metadata_t max_malloc;
F2C* current_handle = nullptr;
current_buffer_metadata_t current_buffer1;
current_buffer_metadata_t current_buffer2;

std::unordered_map<const void*, alloc_metadata_t> allocs;

std::vector<s_smpi_factor_t> parse_factor(const std::string& smpi_coef_string)
{
  std::vector<s_smpi_factor_t> smpi_factor;

  /** Setup the tokenizer that parses the string **/
  using Tokenizer = boost::tokenizer<boost::char_separator<char>>;
  boost::char_separator<char> sep(";");
  boost::char_separator<char> factor_separator(":");
  Tokenizer tokens(smpi_coef_string, sep);

  /**
   * Iterate over patterns like A:B:C:D;E:F;G:H
   * These will be broken down into:
   * A --> B, C, D
   * E --> F
   * G --> H
   */
  for (Tokenizer::iterator token_iter = tokens.begin(); token_iter != tokens.end(); ++token_iter) {
    XBT_DEBUG("token : %s", token_iter->c_str());
    Tokenizer factor_values(*token_iter, factor_separator);
    s_smpi_factor_t fact;
    xbt_assert(factor_values.begin() != factor_values.end(), "Malformed radical for smpi factor: '%s'",
               smpi_coef_string.c_str());
    unsigned int iteration = 0;
    for (Tokenizer::iterator factor_iter = factor_values.begin(); factor_iter != factor_values.end(); ++factor_iter) {
      iteration++;

      if (factor_iter == factor_values.begin()) { /* first element */
        try {
          fact.factor = std::stoi(*factor_iter);
        } catch (const std::invalid_argument&) {
          throw std::invalid_argument(std::string("Invalid factor in chunk ") + std::to_string(smpi_factor.size() + 1) +
                                      ": " + *factor_iter);
        }
      } else {
        try {
          fact.values.push_back(xbt_parse_get_time(surf_parsed_filename, surf_parse_lineno, *factor_iter, ""));
        } catch (const std::invalid_argument&) {
          throw std::invalid_argument(std::string("Invalid factor value ") + std::to_string(iteration) + " in chunk " +
                                      std::to_string(smpi_factor.size() + 1) + ": " + *factor_iter);
        }
      }
    }

    smpi_factor.push_back(fact);
    XBT_DEBUG("smpi_factor:\t%zu : %zu values, first: %f", fact.factor, smpi_factor.size(), fact.values[0]);
  }
  std::sort(smpi_factor.begin(), smpi_factor.end(), [](const s_smpi_factor_t &pa, const s_smpi_factor_t &pb) {
    return (pa.factor < pb.factor);
  });
  for (auto const& fact : smpi_factor) {
    XBT_DEBUG("smpi_factor:\t%zu : %zu values, first: %f", fact.factor, smpi_factor.size() ,fact.values[0]);
  }
  smpi_factor.shrink_to_fit();

  return smpi_factor;
}

void add_benched_time(double time){
  total_benched_time += time;
}

void account_malloc_size(size_t size, const std::string& file, int line, void* ptr)
{
  if (smpi_cfg_display_alloc()) {
    alloc_metadata_t metadata;
    metadata.size = size;
    metadata.line = line;
    metadata.numcall = 1;
    metadata.file    = file;
    allocs.emplace(ptr, metadata);

    total_malloc_size += size;
    if(size > max_malloc.size){
      max_malloc.size = size;
      max_malloc.line = line;
      max_malloc.numcall = 1;
      max_malloc.file    = file;
    } else if (size == max_malloc.size && max_malloc.line == line && max_malloc.file == file) {
      max_malloc.numcall++;
    }
  }
}

void account_shared_size(size_t size){
  if (smpi_cfg_display_alloc()) {
    total_shared_size += size;
    total_shared_calls++;
  }
}

void print_time_analysis(double global_time){
  if (simgrid::config::get_value<bool>("smpi/display-timing")) {
    XBT_INFO("Simulated time: %g seconds. \n\n"
             "The simulation took %g seconds (after parsing and platform setup)\n"
             "%g seconds were actual computation of the application",
             simgrid_get_clock(), global_time, total_benched_time);
    if (total_benched_time/global_time>=0.75)
      XBT_INFO("More than 75%% of the time was spent inside the application code.\n"
    "You may want to use sampling functions or trace replay to reduce this.");
  }
}

static void print_leaked_handles()
{
  // Put the leaked non-default handles in a vector to sort them by id
  std::vector<std::pair<unsigned int, smpi::F2C*>> handles;
  if (simgrid::smpi::F2C::lookup() != nullptr)
    std::copy_if(simgrid::smpi::F2C::lookup()->begin(), simgrid::smpi::F2C::lookup()->end(),
                 std::back_inserter(handles),
                 [](auto const& entry) { return entry.first >= simgrid::smpi::F2C::get_num_default_handles(); });
  if (handles.empty())
    return;

  auto max            = static_cast<unsigned long>(simgrid::config::get_value<int>("smpi/list-leaks"));
  std::string message = "Probable memory leaks in your code: SMPI detected %zu unfreed MPI handles :";
  if (max == 0)
    message += "\nHINT : Display types and addresses (n max) with --cfg=smpi/list-leaks:n.\n"
               "Running smpirun with -wrapper \"valgrind --leak-check=full\" can provide more information";
  XBT_INFO(message.c_str(), handles.size());
  if (max == 0)
    return;

  // we cannot trust F2C::lookup()->size() > F2C::get_num_default_handles() because some default handles are already
  // freed at this point
  bool display_advice = false;
  std::map<std::string, int, std::less<>> count;
  for (const auto& elem : handles) {
    std::string key = elem.second->name();
    if ((not xbt_log_no_loc) && (not elem.second->call_location().empty()))
      key += " at " + elem.second->call_location();
    else
      display_advice = true;
    auto result      = count.emplace(key, 1);
    if (result.second == false)
      result.first->second++;
  }
  if (display_advice)
    XBT_WARN("To get more information (location of allocations), compile your code with -trace-call-location flag of "
             "smpicc/f90");
  unsigned int i = 0;
  for (const auto& p : count) {
    if (p.second == 1)
      XBT_INFO("leaked handle of type %s", p.first.c_str());
    else
      XBT_INFO("%d leaked handles of type %s", p.second, p.first.c_str());
    i++;
    if (i == max)
      break;
  }
  if (max < count.size())
    XBT_INFO("(%lu more handle leaks hidden as you wanted to see only %lu of them)", count.size() - max, max);
}

static void print_leaked_buffers()
{
  if (allocs.empty())
    return;

  auto max            = static_cast<unsigned long>(simgrid::config::get_value<int>("smpi/list-leaks"));
  std::string message = "Probable memory leaks in your code: SMPI detected %zu unfreed buffers :";
  if (max == 0)
    message += "display types and addresses (n max) with --cfg=smpi/list-leaks:n.\nRunning smpirun with -wrapper "
               "\"valgrind --leak-check=full\" can provide more information";
  XBT_INFO(message.c_str(), allocs.size());

  if (max == 0)
    return;

  // gather by allocation origin (only one group reported in case of no-loc or if trace-call-location is not used)
  struct buff_leak {
    int count;
    size_t total_size;
    size_t min_size;
    size_t max_size;
  };
  std::map<std::string, struct buff_leak, std::less<>> leaks_aggreg;
  for (const auto& elem : allocs) {
    std::string key = "leaked allocations";
    if (not xbt_log_no_loc)
      key       = elem.second.file + ":" + std::to_string(elem.second.line) + " : " + key;
    auto result = leaks_aggreg.emplace(key, buff_leak{1, elem.second.size, elem.second.size, elem.second.size});
    if (result.second == false) {
      result.first->second.count++;
      result.first->second.total_size += elem.second.size;
      if (elem.second.size > result.first->second.max_size)
        result.first->second.max_size = elem.second.size;
      else if (elem.second.size < result.first->second.min_size)
        result.first->second.min_size = elem.second.size;
    }
  }
  // now we can order by total size.
  std::vector<std::pair<std::string, buff_leak>> leaks(leaks_aggreg.begin(), leaks_aggreg.end());
  std::sort(leaks.begin(), leaks.end(),
            [](auto const& a, auto const& b) { return a.second.total_size > b.second.total_size; });

  unsigned int i = 0;
  for (const auto& p : leaks) {
    if (p.second.min_size == p.second.max_size)
      XBT_INFO("%s of total size %zu, called %d times, each with size %zu", p.first.c_str(), p.second.total_size,
               p.second.count, p.second.min_size);
    else
      XBT_INFO("%s of total size %zu, called %d times, with minimum size %zu and maximum size %zu", p.first.c_str(),
               p.second.total_size, p.second.count, p.second.min_size, p.second.max_size);
    i++;
    if (i == max)
      break;
  }
  if (max < leaks_aggreg.size())
    XBT_INFO("(more buffer leaks hidden as you wanted to see only %lu of them)", max);
}

void print_memory_analysis()
{
  if (smpi_cfg_display_alloc()) {
    print_leaked_handles();
    print_leaked_buffers();

    if(total_malloc_size != 0)
      XBT_INFO("Memory Usage: Simulated application allocated %lu bytes during its lifetime through malloc/calloc calls.\n"
             "Largest allocation at once from a single process was %zu bytes, at %s:%d. It was called %u times during the whole simulation.\n"
             "If this is too much, consider sharing allocations for computation buffers.\n"
             "This can be done automatically by setting --cfg=smpi/auto-shared-malloc-thresh to the minimum size wanted size (this can alter execution if data content is necessary)\n",
             total_malloc_size, max_malloc.size, simgrid::xbt::Path(max_malloc.file).get_base_name().c_str(), max_malloc.line, max_malloc.numcall
      );
    else
      XBT_INFO("Allocations analysis asked, but 0 bytes were allocated through malloc/calloc calls intercepted by SMPI.\n"
               "Either code is using other ways of allocating memory, or it was built with SMPI_NO_OVERRIDE_MALLOC");
    if(total_shared_size != 0)
      XBT_INFO("%lu bytes were automatically shared between processes, in %u calls\n", total_shared_size, total_shared_calls);
  }
}

void set_current_handle(F2C* handle){
  current_handle=handle;
}

void print_current_handle(){
  if(current_handle){
    if(current_handle->call_location().empty())
      XBT_INFO("To get handle location information, pass -trace-call-location flag to smpicc/f90 as well");
    else
      XBT_INFO("Handle %s was allocated by a call at %s", current_handle->name().c_str(),
               (char*)(current_handle->call_location().c_str()));
  }
}

void set_current_buffer(int i, const char* name, const void* buf){
  //clear previous one
  if(i==1){
    if(not current_buffer1.name.empty()){
      current_buffer1.name="";
    }
    if(not current_buffer2.name.empty()){
      current_buffer2.name="";
    }
  }
  auto meta = allocs.find(buf);
  if (meta == allocs.end()) {
    XBT_DEBUG("Buffer %p was not allocated with malloc/calloc", buf);
    return;
  }
  if(i==1){
    current_buffer1.alloc = meta->second;
    current_buffer1.name = name;
  }else{
    current_buffer2.alloc=meta->second;
    current_buffer2.name=name;
  }
}

void print_buffer_info(){
    if(not current_buffer1.name.empty())
      XBT_INFO("Buffer %s was allocated from %s line %d, with size %zu", current_buffer1.name.c_str(), current_buffer1.alloc.file.c_str(), current_buffer1.alloc.line, current_buffer1.alloc.size);
    if(not current_buffer2.name.empty())
      XBT_INFO("Buffer %s was allocated from %s line %d, with size %zu", current_buffer2.name.c_str(), current_buffer2.alloc.file.c_str(), current_buffer2.alloc.line, current_buffer2.alloc.size);    
}

size_t get_buffer_size(const void* buf){
  auto meta = allocs.find(buf);
  if (meta == allocs.end()) {
    //we don't know this buffer (on stack or feature disabled), assume it's fine.
    return  std::numeric_limits<std::size_t>::max();
  }
  return meta->second.size;
}

void account_free(const void* ptr){
  if (smpi_cfg_display_alloc()) {
    allocs.erase(ptr);
  }
}

}
}
} // namespace simgrid
