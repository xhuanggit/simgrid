/* Copyright (c) 2019-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "DiskImpl.hpp"

#include "simgrid/s4u/Engine.hpp"
#include "src/kernel/EngineImpl.hpp"
#include "src/kernel/lmm/maxmin.hpp"

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(res_disk, ker_resource, "Disk resources, fuelling I/O activities");

namespace simgrid {
namespace kernel {
namespace resource {

/*********
 * Model *
 *********/

DiskModel::DiskModel(const std::string& name) : Model(name)
{
  set_maxmin_system(new lmm::System(true /* selective update */));
}

/************
 * Resource *
 ************/
DiskImpl* DiskImpl::set_host(s4u::Host* host)
{
  host_ = host;
  return this;
}

DiskImpl* DiskImpl::set_read_bandwidth(double read_bw)
{
  read_bw_ = read_bw;
  return this;
}

DiskImpl* DiskImpl::set_write_bandwidth(double write_bw)
{
  write_bw_ = write_bw;
  return this;
}

DiskImpl* DiskImpl::set_read_constraint(lmm::Constraint* constraint_read)
{
  constraint_read_ = constraint_read;
  return this;
}

DiskImpl* DiskImpl::set_write_constraint(lmm::Constraint* constraint_write)
{
  constraint_write_ = constraint_write;
  return this;
}

/** @brief Fire the required callbacks and destroy the object
 *
 * Don't delete directly a Disk, call d->destroy() instead.
 */
void DiskImpl::destroy()
{
  s4u::Disk::on_destruction(this->piface_);
  delete this;
}

bool DiskImpl::is_used() const
{
  return get_model()->get_maxmin_system()->constraint_used(get_constraint());
}

void DiskImpl::apply_event(kernel::profile::Event* /*event*/, double /*value*/)
{
  THROW_UNIMPLEMENTED;
}

void DiskImpl::turn_on()
{
  if (not is_on()) {
    Resource::turn_on();
    s4u::Disk::on_state_change(this->piface_);
  }
}
void DiskImpl::turn_off()
{
  if (is_on()) {
    Resource::turn_off();
    s4u::Disk::on_state_change(this->piface_);
  }
}

void DiskImpl::seal()
{
  xbt_assert(this->get_model(), "Cannot seal Disk (%s) without setting the model first", this->get_cname());
  lmm::System* maxmin_system = get_model()->get_maxmin_system();
  this->set_read_constraint(maxmin_system->constraint_new(this, read_bw_))
      ->set_write_constraint(maxmin_system->constraint_new(this, write_bw_))
      ->set_constraint(maxmin_system->constraint_new(this, std::max(read_bw_, write_bw_)));
  XBT_DEBUG("Create resource with read_bw '%f' write_bw '%f'", read_bw_, write_bw_);
  Resource::seal();
  turn_on();
}
xbt::signal<void(DiskAction const&, Action::State, Action::State)> DiskAction::on_state_change;

/**********
 * Action *
 **********/
void DiskAction::set_state(Action::State state)
{
  Action::State old = get_state();
  Action::set_state(state);
  on_state_change(*this, old, state);
}
} // namespace resource
} // namespace kernel
} // namespace simgrid
