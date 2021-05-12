/* Copyright (c) 2013-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "cpu_interface.hpp"
#include "cpu_ti.hpp"
#include "src/kernel/resource/profile/Profile.hpp"
#include "src/surf/surf_interface.hpp"
#include "surf/surf.hpp"

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(res_cpu, ker_resource, "CPU resource, fueling execution activites");

namespace simgrid {
namespace kernel {
namespace resource {

/*********
 * Model *
 *********/

void CpuModel::update_actions_state_lazy(double now, double /*delta*/)
{
  while (not get_action_heap().empty() && double_equals(get_action_heap().top_date(), now, sg_surf_precision)) {
    auto* action = static_cast<CpuAction*>(get_action_heap().pop());
    XBT_DEBUG("Something happened to action %p", action);

    action->finish(kernel::resource::Action::State::FINISHED);
    XBT_DEBUG("Action %p finished", action);
  }
}

void CpuModel::update_actions_state_full(double /*now*/, double delta)
{
  for (auto it = std::begin(*get_started_action_set()); it != std::end(*get_started_action_set());) {
    auto& action = static_cast<CpuAction&>(*it);
    ++it; // increment iterator here since the following calls to action.finish() may invalidate it

    action.update_remains(action.get_variable()->get_value() * delta);
    action.update_max_duration(delta);

    if (((action.get_remains_no_update() <= 0) && (action.get_variable()->get_penalty() > 0)) ||
        ((action.get_max_duration() != NO_MAX_DURATION) && (action.get_max_duration() <= 0))) {
      action.finish(Action::State::FINISHED);
    }
  }
}

/************
 * Resource *
 ************/
Cpu::Cpu(s4u::Host* host, const std::vector<double>& speed_per_pstate)
    : Resource_T(host->get_cname()), piface_(host), speed_per_pstate_(speed_per_pstate)
{
  speed_.scale = 1;
  speed_.peak     = speed_per_pstate_.front();
  host->pimpl_cpu = this;
}

void Cpu::reset_vcpu(Cpu* that)
{
  this->pstate_ = that->pstate_;
  this->speed_  = that->speed_;
  this->speed_per_pstate_.clear();
  this->speed_per_pstate_.assign(that->speed_per_pstate_.begin(), that->speed_per_pstate_.end());
}

Cpu* Cpu::set_pstate(int pstate_index)
{
  xbt_assert(pstate_index <= static_cast<int>(speed_per_pstate_.size()),
             "Invalid parameters for CPU %s (pstate %d > length of pstates %d). Please fix your platform file, or your "
             "call to change the pstate.",
             get_cname(), pstate_index, static_cast<int>(speed_per_pstate_.size()));

  double new_peak_speed = speed_per_pstate_[pstate_index];
  pstate_ = pstate_index;
  speed_.peak = new_peak_speed;

  on_speed_change();
  return this;
}

double Cpu::get_pstate_peak_speed(int pstate_index) const
{
  xbt_assert((pstate_index <= static_cast<int>(speed_per_pstate_.size())),
             "Invalid parameters (pstate index out of bounds)");

  return speed_per_pstate_[pstate_index];
}

void Cpu::on_speed_change()
{
  s4u::Host::on_speed_change(*piface_);
}

Cpu* Cpu::set_core_count(int core_count)
{
  xbt_assert(not is_sealed(), "Core count cannot be changed once CPU has been sealed");
  xbt_assert(core_count > 0, "Host %s must have at least one core, not 0.", piface_->get_cname());
  if (dynamic_cast<CpuTiModel*>(get_model()) != nullptr)
    xbt_assert(core_count == 1, "Multi-core not handled by this model yet");

  core_count_ = core_count;
  return this;
}

int Cpu::get_core_count()
{
  return core_count_;
}

Cpu* Cpu::set_speed_profile(kernel::profile::Profile* profile)
{
  if (profile) {
    xbt_assert(speed_.event == nullptr, "Cannot set a second speed trace to Host %s", piface_->get_cname());
    speed_.event = profile->schedule(&profile::future_evt_set, this);
  }
  return this;
}

void Cpu::seal()
{
  lmm::System* lmm = get_model()->get_maxmin_system();
  if (dynamic_cast<CpuTiModel*>(get_model()) == nullptr)
    this->set_constraint(lmm->constraint_new(this, core_count_ * speed_per_pstate_.front()));
  Resource::seal();
}

/**********
 * Action *
 **********/

void CpuAction::update_remains_lazy(double now)
{
  xbt_assert(get_state_set() == get_model()->get_started_action_set(),
             "You're updating an action that is not running.");
  xbt_assert(get_sharing_penalty() > 0, "You're updating an action that seems suspended.");

  double delta = now - get_last_update();

  if (get_remains_no_update() > 0) {
    XBT_DEBUG("Updating action(%p): remains was %f, last_update was: %f", this, get_remains_no_update(),
              get_last_update());
    update_remains(get_last_value() * delta);

    XBT_DEBUG("Updating action(%p): remains is now %f", this, get_remains_no_update());
  }

  set_last_update();
  set_last_value(get_variable()->get_value());
}

xbt::signal<void(CpuAction const&, Action::State)> CpuAction::on_state_change;

void CpuAction::suspend(){
  Action::State previous = get_state();
  on_state_change(*this, previous);
  Action::suspend();
}

void CpuAction::resume(){
  Action::State previous = get_state();
  on_state_change(*this, previous);
  Action::resume();
}

void CpuAction::set_state(Action::State state)
{
  Action::State previous = get_state();
  Action::set_state(state);
  on_state_change(*this, previous);
}

/** @brief returns a list of all CPUs that this action is using */
std::list<Cpu*> CpuAction::cpus() const
{
  std::list<Cpu*> retlist;
  int llen = get_variable()->get_number_of_constraint();

  for (int i = 0; i < llen; i++) {
    /* Beware of composite actions: ptasks put links and cpus together */
    // extra pb: we cannot dynamic_cast from void*...
    Resource* resource = get_variable()->get_constraint(i)->get_id();
    auto* cpu          = dynamic_cast<Cpu*>(resource);
    if (cpu != nullptr)
      retlist.push_back(cpu);
  }

  return retlist;
}
} // namespace resource
} // namespace kernel
} // namespace simgrid
