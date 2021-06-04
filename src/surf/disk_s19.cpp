/* Copyright (c) 2013-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "disk_s19.hpp"
#include "simgrid/kernel/routing/NetPoint.hpp"
#include "simgrid/kernel/routing/NetZoneImpl.hpp"
#include "simgrid/s4u/Engine.hpp"
#include "simgrid/s4u/Host.hpp"
#include "src/kernel/EngineImpl.hpp"
#include "src/kernel/lmm/maxmin.hpp"
#include "src/kernel/resource/profile/Event.hpp"
#include "src/surf/xml/platf.hpp"
#include "surf/surf.hpp"

XBT_LOG_EXTERNAL_DEFAULT_CATEGORY(res_disk);

/*********
 * Model *
 *********/

void surf_disk_model_init_default()
{
  auto disk_model = std::make_shared<simgrid::kernel::resource::DiskS19Model>("Disk");
  simgrid::kernel::EngineImpl::get_instance()->add_model(disk_model);
  simgrid::s4u::Engine::get_instance()->get_netzone_root()->get_impl()->set_disk_model(disk_model);
}

namespace simgrid {
namespace kernel {
namespace resource {

DiskImpl* DiskS19Model::create_disk(const std::string& name, double read_bandwidth, double write_bandwidth)
{
  return (new DiskS19(name, read_bandwidth, write_bandwidth))->set_model(this);
}

void DiskS19Model::update_actions_state(double /*now*/, double delta)
{
  for (auto it = std::begin(*get_started_action_set()); it != std::end(*get_started_action_set());) {
    auto& action = *it;
    ++it; // increment iterator here since the following calls to action.finish() may invalidate it
    action.update_remains(rint(action.get_variable()->get_value() * delta));
    action.update_max_duration(delta);

    if (((action.get_remains_no_update() <= 0) && (action.get_variable()->get_penalty() > 0)) ||
        ((action.get_max_duration() != NO_MAX_DURATION) && (action.get_max_duration() <= 0))) {
      action.finish(Action::State::FINISHED);
    }
  }
}

DiskAction* DiskS19Model::io_start(const DiskImpl* disk, sg_size_t size, s4u::Io::OpType type)
{
  auto* action = new DiskS19Action(this, static_cast<double>(size), not disk->is_on());
  get_maxmin_system()->expand(disk->get_constraint(), action->get_variable(), 1.0);
  switch (type) {
    case s4u::Io::OpType::READ:
      get_maxmin_system()->expand(disk->get_read_constraint(), action->get_variable(), 1.0);
      break;
    case s4u::Io::OpType::WRITE:
      get_maxmin_system()->expand(disk->get_write_constraint(), action->get_variable(), 1.0);
      break;
    default:
      THROW_UNIMPLEMENTED;
  }
  return action;
}

/************
 * Resource *
 ************/
void DiskS19::update_penalties(double delta) const
{
  const kernel::lmm::Element* elem     = nullptr;
  const kernel::lmm::Element* nextelem = nullptr;
  int numelem                          = 0;
  while (const auto* var = get_constraint()->get_variable_safe(&elem, &nextelem, &numelem)) {
    auto* action = static_cast<DiskS19Action*>(var->get_id());
    action->sharing_penalty_ += delta;
    if (not action->is_suspended())
      get_model()->get_maxmin_system()->update_variable_penalty(action->get_variable(), action->sharing_penalty_);
  }
}

void DiskS19::set_read_bandwidth(double value)
{
  read_bw_.peak = value;

  get_model()->get_maxmin_system()->update_constraint_bound(get_constraint(), read_bw_.peak * read_bw_.scale);

  double delta = 1.0 / value - 1.0 / (read_bw_.peak * read_bw_.scale);
  update_penalties(delta);
}

void DiskS19::set_write_bandwidth(double value)
{
  write_bw_.peak = value;

  get_model()->get_maxmin_system()->update_constraint_bound(get_constraint(), write_bw_.peak * write_bw_.scale);

  double delta = 1.0 / value - 1.0 / (write_bw_.peak * write_bw_.scale);
  update_penalties(delta);
}

void DiskS19::apply_event(kernel::profile::Event* triggered, double value)
{
  /* Find out which of my iterators was triggered, and react accordingly */
  if (triggered == read_bw_.event) {
    set_read_bandwidth(value);
    tmgr_trace_event_unref(&read_bw_.event);

  } else if (triggered == write_bw_.event) {
    set_write_bandwidth(value);
    tmgr_trace_event_unref(&write_bw_.event);

  } else if (triggered == state_event_) {
    if (value > 0)
      turn_on();
    else
      turn_off();
    tmgr_trace_event_unref(&state_event_);
  } else {
    xbt_die("Unknown event!\n");
  }

  XBT_DEBUG("There was a resource state event, need to update actions related to the constraint (%p)",
            get_constraint());
}

/**********
 * Action *
 **********/

DiskS19Action::DiskS19Action(Model* model, double cost, bool failed)
    : DiskAction(model, cost, failed, model->get_maxmin_system()->variable_new(this, 1.0, -1.0, 3))
{
}

void DiskS19Action::update_remains_lazy(double /*now*/)
{
  THROW_IMPOSSIBLE;
}
} // namespace resource
} // namespace kernel
} // namespace simgrid
