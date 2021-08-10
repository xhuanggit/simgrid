/* Copyright (c) 2016-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "src/kernel/EngineImpl.hpp"
#include "mc/mc.h"
#include "simgrid/Exception.hpp"
#include "simgrid/kernel/Timer.hpp"
#include "simgrid/kernel/routing/NetPoint.hpp"
#include "simgrid/kernel/routing/NetZoneImpl.hpp"
#include "simgrid/s4u/Host.hpp"
#include "simgrid/sg_config.hpp"
#include "src/include/surf/surf.hpp" //get_clock() and surf_solve()
#include "src/kernel/resource/DiskImpl.hpp"
#include "src/mc/mc_record.hpp"
#include "src/mc/mc_replay.hpp"
#include "src/simix/smx_private.hpp"
#include "src/smpi/include/smpi_actor.hpp"
#include "src/surf/network_interface.hpp"
#include "src/surf/xml/platf.hpp" // FIXME: KILLME. There must be a better way than mimicking XML here

#include <boost/algorithm/string/predicate.hpp>
#ifndef _WIN32
#include <dlfcn.h>
#endif /* _WIN32 */

XBT_LOG_NEW_DEFAULT_CATEGORY(ker_engine, "Logging specific to Engine (kernel)");

namespace simgrid {
namespace kernel {

config::Flag<double> cfg_breakpoint{"debug/breakpoint",
                                    "When non-negative, raise a SIGTRAP after given (simulated) time", -1.0};
EngineImpl::~EngineImpl()
{
  while (not timer::kernel_timers().empty()) {
    delete timer::kernel_timers().top().second;
    timer::kernel_timers().pop();
  }

  /* Since hosts_ is a std::map, the hosts are destroyed in the lexicographic order, which ensures that the output is
   * reproducible.
   */
  while (not hosts_.empty())
    hosts_.begin()->second->destroy();

  /* Also delete the other data */
  delete netzone_root_;
  for (auto const& kv : netpoints_)
    delete kv.second;

  for (auto const& kv : links_)
    if (kv.second)
      kv.second->destroy();

  for (auto const& kv : mailboxes_)
    delete kv.second;

    /* Free the remaining data structures */
#if SIMGRID_HAVE_MC
  xbt_dynar_free(&actors_vector_);
  xbt_dynar_free(&dead_actors_vector_);
#endif
  /* clear models before freeing handle, network models can use external callback defined in the handle */
  models_prio_.clear();
}

void EngineImpl::load_platform(const std::string& platf)
{
  double start = xbt_os_time();
  if (boost::algorithm::ends_with(platf, ".so") or boost::algorithm::ends_with(platf, ".dylib")) {
#ifdef _WIN32
    xbt_die("loading platform through shared library isn't supported on windows");
#else
    void* handle = dlopen(platf.c_str(), RTLD_LAZY);
    xbt_assert(handle, "Impossible to open platform file: %s", platf.c_str());
    platf_handle_           = std::unique_ptr<void, std::function<int(void*)>>(handle, dlclose);
    using load_fct_t = void (*)(const simgrid::s4u::Engine&);
    auto callable           = (load_fct_t)dlsym(platf_handle_.get(), "load_platform");
    const char* dlsym_error = dlerror();
    xbt_assert(not dlsym_error, "Error: %s", dlsym_error);
    callable(*simgrid::s4u::Engine::get_instance());
#endif /* _WIN32 */
  } else {
    parse_platform_file(platf);
  }

  double end = xbt_os_time();
  XBT_DEBUG("PARSE TIME: %g", (end - start));
}

void EngineImpl::load_deployment(const std::string& file) const
{
  sg_platf_exit();
  sg_platf_init();

  surf_parse_open(file);
  surf_parse();
  surf_parse_close();
}

void EngineImpl::register_function(const std::string& name, const actor::ActorCodeFactory& code)
{
  registered_functions[name] = code;
}
void EngineImpl::register_default(const actor::ActorCodeFactory& code)
{
  default_function = code;
}

void EngineImpl::add_model(std::shared_ptr<resource::Model> model, const std::vector<resource::Model*>& dependencies)
{
  auto model_name = model->get_name();
  xbt_assert(models_prio_.find(model_name) == models_prio_.end(),
             "Model %s already exists, use model.set_name() to change its name", model_name.c_str());

  for (const auto dep : dependencies) {
    xbt_assert(models_prio_.find(dep->get_name()) != models_prio_.end(),
               "Model %s doesn't exists. Impossible to use it as dependency.", dep->get_name().c_str());
  }
  models_.push_back(model.get());
  models_prio_[model_name] = std::move(model);
}

void EngineImpl::add_split_duplex_link(const std::string& name, std::unique_ptr<resource::SplitDuplexLinkImpl> link)
{
  split_duplex_links_[name] = std::move(link);
}

/** Wake up all actors waiting for a Surf action to finish */
void EngineImpl::wake_all_waiting_actors() const
{
  for (auto const& model : models_) {
    XBT_DEBUG("Handling the failed actions (if any)");
    while (auto* action = model->extract_failed_action()) {
      XBT_DEBUG("   Handling Action %p", action);
      if (action->get_activity() != nullptr)
        activity::ActivityImplPtr(action->get_activity())->post();
    }
    XBT_DEBUG("Handling the terminated actions (if any)");
    while (auto* action = model->extract_done_action()) {
      XBT_DEBUG("   Handling Action %p", action);
      if (action->get_activity() == nullptr)
        XBT_DEBUG("probably vcpu's action %p, skip", action);
      else
        activity::ActivityImplPtr(action->get_activity())->post();
    }
  }
}
/**
 * @brief Executes the actors in actors_to_run.
 *
 * The actors in actors_to_run are run (in parallel if possible). On exit, actors_to_run is empty, and actors_that_ran
 * contains the list of actors that just ran.  The two lists are swapped so, be careful when using them before and after
 * a call to this function.
 */
void EngineImpl::run_all_actors()
{
  simix_global->get_context_factory()->run_all();

  actors_to_run_.swap(actors_that_ran_);
  actors_to_run_.clear();
}

actor::ActorImpl* EngineImpl::get_actor_by_pid(aid_t pid)
{
  auto item = actor_list_.find(pid);
  if (item != actor_list_.end())
    return item->second;

  // Search the trash
  for (auto& a : actors_to_destroy_)
    if (a.get_pid() == pid)
      return &a;
  return nullptr; // Not found, even in the trash
}
/** Execute all the tasks that are queued, e.g. `.then()` callbacks of futures. */
bool EngineImpl::execute_tasks()
{
  if (tasks.empty())
    return false;

  std::vector<xbt::Task<void()>> tasksTemp;
  do {
    // We don't want the callbacks to modify the vector we are iterating over:
    tasks.swap(tasksTemp);

    // Execute all the queued tasks:
    for (auto& task : tasksTemp)
      task();

    tasksTemp.clear();
  } while (not tasks.empty());

  return true;
}

void EngineImpl::remove_daemon(actor::ActorImpl* actor)
{
  auto it = daemons_.find(actor);
  xbt_assert(it != daemons_.end(), "The dying daemon is not a daemon after all. Please report that bug.");
  daemons_.erase(it);
}

void EngineImpl::add_actor_to_run_list_no_check(actor::ActorImpl* actor)
{
  XBT_DEBUG("Inserting [%p] %s(%s) in the to_run list", actor, actor->get_cname(), actor->get_host()->get_cname());
  actors_to_run_.push_back(actor);
}

void EngineImpl::add_actor_to_run_list(actor::ActorImpl* actor)
{
  if (std::find(begin(actors_to_run_), end(actors_to_run_), actor) != end(actors_to_run_)) {
    XBT_DEBUG("Actor %s is already in the to_run list", actor->get_cname());
  } else {
    XBT_DEBUG("Inserting [%p] %s(%s) in the to_run list", actor, actor->get_cname(), actor->get_host()->get_cname());
    actors_to_run_.push_back(actor);
  }
}
void EngineImpl::empty_trash()
{
  while (not actors_to_destroy_.empty()) {
    actor::ActorImpl* actor = &actors_to_destroy_.front();
    actors_to_destroy_.pop_front();
    XBT_DEBUG("Getting rid of %s (refcount: %d)", actor->get_cname(), actor->get_refcount());
    intrusive_ptr_release(actor);
  }
#if SIMGRID_HAVE_MC
  xbt_dynar_reset(dead_actors_vector_);
#endif
}

void EngineImpl::display_all_actor_status() const
{
  XBT_INFO("%zu actors are still running, waiting for something.", actor_list_.size());
  /*  List the actors and their state */
  XBT_INFO("Legend of the following listing: \"Actor <pid> (<name>@<host>): <status>\"");
  for (auto const& kv : actor_list_) {
    actor::ActorImpl* actor = kv.second;

    if (actor->waiting_synchro_) {
      const char* synchro_description = "unknown";

      if (boost::dynamic_pointer_cast<kernel::activity::ExecImpl>(actor->waiting_synchro_) != nullptr)
        synchro_description = "execution";

      if (boost::dynamic_pointer_cast<kernel::activity::CommImpl>(actor->waiting_synchro_) != nullptr)
        synchro_description = "communication";

      if (boost::dynamic_pointer_cast<kernel::activity::SleepImpl>(actor->waiting_synchro_) != nullptr)
        synchro_description = "sleeping";

      if (boost::dynamic_pointer_cast<kernel::activity::RawImpl>(actor->waiting_synchro_) != nullptr)
        synchro_description = "synchronization";

      if (boost::dynamic_pointer_cast<kernel::activity::IoImpl>(actor->waiting_synchro_) != nullptr)
        synchro_description = "I/O";

      XBT_INFO("Actor %ld (%s@%s): waiting for %s activity %#zx (%s) in state %d to finish", actor->get_pid(),
               actor->get_cname(), actor->get_host()->get_cname(), synchro_description,
               (xbt_log_no_loc ? (size_t)0xDEADBEEF : (size_t)actor->waiting_synchro_.get()),
               actor->waiting_synchro_->get_cname(), (int)actor->waiting_synchro_->state_);
    } else {
      XBT_INFO("Actor %ld (%s@%s) simcall %s", actor->get_pid(), actor->get_cname(), actor->get_host()->get_cname(),
               SIMIX_simcall_name(actor->simcall_));
    }
  }
}

void EngineImpl::run()
{
  if (MC_record_replay_is_active()) {
    mc::replay(MC_record_path());
    empty_trash();
    return;
  }

  double time = 0;

  do {
    XBT_DEBUG("New Schedule Round; size(queue)=%zu", actors_to_run_.size());

    if (cfg_breakpoint >= 0.0 && surf_get_clock() >= cfg_breakpoint) {
      XBT_DEBUG("Breakpoint reached (%g)", cfg_breakpoint.get());
      cfg_breakpoint = -1.0;
#ifdef SIGTRAP
      std::raise(SIGTRAP);
#else
      std::raise(SIGABRT);
#endif
    }

    execute_tasks();

    while (not actors_to_run_.empty()) {
      XBT_DEBUG("New Sub-Schedule Round; size(queue)=%zu", actors_to_run_.size());

      /* Run all actors that are ready to run, possibly in parallel */
      run_all_actors();

      /* answer sequentially and in a fixed arbitrary order all the simcalls that were issued during that sub-round */

      /* WARNING, the order *must* be fixed or you'll jeopardize the simulation reproducibility (see RR-7653) */

      /* Here, the order is ok because:
       *
       *   Short proof: only maestro adds stuff to the actors_to_run array, so the execution order of user contexts do
       *   not impact its order.
       *
       *   Long proof: actors remain sorted through an arbitrary (implicit, complex but fixed) order in all cases.
       *
       *   - if there is no kill during the simulation, actors remain sorted according by their PID.
       *     Rationale: This can be proved inductively.
       *        Assume that actors_to_run is sorted at a beginning of one round (it is at round 0: the deployment file
       *        is parsed linearly).
       *        Let's show that it is still so at the end of this round.
       *        - if an actor is added when being created, that's from maestro. It can be either at startup
       *          time (and then in PID order), or in response to a process_create simcall. Since simcalls are handled
       *          in arbitrary order (inductive hypothesis), we are fine.
       *        - If an actor is added because it's getting killed, its subsequent actions shouldn't matter
       *        - If an actor gets added to actors_to_run because one of their blocking action constituting the meat
       *          of a simcall terminates, we're still good. Proof:
       *          - You are added from ActorImpl::simcall_answer() only. When this function is called depends on the
       *            resource kind (network, cpu, disk, whatever), but the same arguments hold. Let's take communications
       *            as an example.
       *          - For communications, this function is called from SIMIX_comm_finish().
       *            This function itself don't mess with the order since simcalls are handled in FIFO order.
       *            The function is called:
       *            - before the comm starts (invalid parameters, or resource already dead or whatever).
       *              The order then trivial holds since maestro didn't interrupt its handling of the simcall yet
       *            - because the communication failed or were canceled after startup. In this case, it's called from
       *              the function we are in, by the chunk:
       *                       set = model->states.failed_action_set;
       *                       while ((synchro = extract(set)))
       *                          SIMIX_simcall_post((smx_synchro_t) synchro->data);
       *              This order is also fixed because it depends of the order in which the surf actions were
       *              added to the system, and only maestro can add stuff this way, through simcalls.
       *              We thus use the inductive hypothesis once again to conclude that the order in which synchros are
       *              popped out of the set does not depend on the user code's execution order.
       *            - because the communication terminated. In this case, synchros are served in the order given by
       *                       set = model->states.done_action_set;
       *                       while ((synchro = extract(set)))
       *                          SIMIX_simcall_post((smx_synchro_t) synchro->data);
       *              and the argument is very similar to the previous one.
       *            So, in any case, the orders of calls to CommImpl::finish() do not depend on the order in which user
       *            actors are executed.
       *          So, in any cases, the orders of actors within actors_to_run do not depend on the order in which
       *          user actors were executed previously.
       *     So, if there is no killing in the simulation, the simulation reproducibility is not jeopardized.
       *   - If there is some actor killings, the order is changed by this decision that comes from user-land
       *     But this decision may not have been motivated by a situation that were different because the simulation is
       *     not reproducible.
       *     So, even the order change induced by the actor killing is perfectly reproducible.
       *
       *   So science works, bitches [http://xkcd.com/54/].
       *
       *   We could sort the actors_that_ran array completely so that we can describe the order in which simcalls are
       *   handled (like "according to the PID of issuer"), but it's not mandatory (order is fixed already even if
       *   unfriendly).
       *   That would thus be a pure waste of time.
       */

      for (auto const& actor : actors_that_ran_) {
        if (actor->simcall_.call_ != simix::Simcall::NONE) {
          actor->simcall_handle(0);
        }
      }

      execute_tasks();
      do {
        wake_all_waiting_actors();
      } while (execute_tasks());

      /* If only daemon actors remain, cancel their actions, mark them to die and reschedule them */
      if (actor_list_.size() == daemons_.size())
        for (auto const& dmon : daemons_) {
          XBT_DEBUG("Kill %s", dmon->get_cname());
          simix_global->get_maestro()->kill(dmon);
        }
    }

    time = timer::Timer::next();
    if (time > -1.0 || not actor_list_.empty()) {
      XBT_DEBUG("Calling surf_solve");
      time = surf_solve(time);
      XBT_DEBUG("Moving time ahead : %g", time);
    }

    /* Notify all the hosts that have failed */
    /* FIXME: iterate through the list of failed host and mark each of them */
    /* as failed. On each host, signal all the running actors with host_fail */

    // Execute timers and tasks until there isn't anything to be done:
    bool again = false;
    do {
      again = timer::Timer::execute_all();
      if (execute_tasks())
        again = true;
      wake_all_waiting_actors();
    } while (again);

    /* Clean actors to destroy */
    empty_trash();

    XBT_DEBUG("### time %f, #actors %zu, #to_run %zu", time, actor_list_.size(), actors_to_run_.size());

    if (time < 0. && actors_to_run_.empty() && not actor_list_.empty()) {
      if (actor_list_.size() <= daemons_.size()) {
        XBT_CRITICAL("Oops! Daemon actors cannot do any blocking activity (communications, synchronization, etc) "
                     "once the simulation is over. Please fix your on_exit() functions.");
      } else {
        XBT_CRITICAL("Oops! Deadlock or code not perfectly clean.");
      }
      display_all_actor_status();
      simgrid::s4u::Engine::on_deadlock();
      for (auto const& kv : actor_list_) {
        XBT_DEBUG("Kill %s", kv.second->get_cname());
        simix_global->get_maestro()->kill(kv.second);
      }
    }
  } while (time > -1.0 || has_actors_to_run());

  if (not actor_list_.empty())
    THROW_IMPOSSIBLE;

  simgrid::s4u::Engine::on_simulation_end();
}
} // namespace kernel
} // namespace simgrid
