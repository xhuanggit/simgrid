/* Copyright (c) 2019-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "simgrid/kernel/resource/Action.hpp"
#include "simgrid/kernel/resource/Model.hpp"
#include "simgrid/kernel/resource/Resource.hpp"
#include "simgrid/s4u/Disk.hpp"
#include "simgrid/s4u/Io.hpp"
#include "src/surf/surf_interface.hpp"
#include <xbt/PropertyHolder.hpp>

#include <map>

#ifndef DISK_INTERFACE_HPP_
#define DISK_INTERFACE_HPP_

/*********
 * Model *
 *********/

namespace simgrid {
namespace kernel {
namespace resource {
/***********
 * Classes *
 ***********/

class DiskAction;

/*********
 * Model *
 *********/
class DiskModel : public Model {
public:
  explicit DiskModel(const std::string& name);
  DiskModel(const DiskModel&) = delete;
  DiskModel& operator=(const DiskModel&) = delete;

  virtual DiskImpl* create_disk(const std::string& name, double read_bandwidth, double write_bandwidth) = 0;
};

/************
 * Resource *
 ************/
class DiskImpl : public Resource_T<DiskImpl>, public xbt::PropertyHolder {
  s4u::Host* host_           = nullptr;
  s4u::Disk piface_;
  double read_bw_ = -1.0;
  double write_bw_ = 1.0;
  lmm::Constraint* constraint_write_ = nullptr; /* Constraint for maximum write bandwidth*/
  lmm::Constraint* constraint_read_ = nullptr;  /* Constraint for maximum read bandwidth*/

protected:
  ~DiskImpl() override = default; // Disallow direct deletion. Call destroy() instead.

public:
  DiskImpl(const std::string& name, double read_bandwidth, double write_bandwidth)
      : Resource_T(name), piface_(name, this), read_bw_(read_bandwidth), write_bw_(write_bandwidth)
  {
  }
  DiskImpl(const DiskImpl&) = delete;
  DiskImpl& operator=(const DiskImpl&) = delete;

  /** @brief Public interface */
  const s4u::Disk* get_iface() const { return &piface_; }
  s4u::Disk* get_iface() { return &piface_; }
  DiskImpl* set_host(s4u::Host* host);
  s4u::Host* get_host() const { return host_; }

  DiskImpl* set_read_bandwidth(double read_bw);
  double get_read_bandwidth() const { return read_bw_; }

  DiskImpl* set_write_bandwidth(double write_bw);
  double get_write_bandwidth() const { return write_bw_; }

  DiskImpl* set_read_constraint(lmm::Constraint* constraint_read);
  lmm::Constraint* get_read_constraint() const { return constraint_read_; }

  DiskImpl* set_write_constraint(lmm::Constraint* constraint_write);
  lmm::Constraint* get_write_constraint() const { return constraint_write_; }

  /** @brief Check if the Disk is used (if an action currently uses its resources) */
  bool is_used() const override;
  void apply_event(profile::Event* event, double value) override;
  void turn_on() override;
  void turn_off() override;

  void seal() override;
  void destroy(); // Must be called instead of the destructor
  virtual DiskAction* io_start(sg_size_t size, s4u::Io::OpType type) = 0;
  virtual DiskAction* read(sg_size_t size)                           = 0;
  virtual DiskAction* write(sg_size_t size)                          = 0;
};

/**********
 * Action *
 **********/

class DiskAction : public Action {
public:
  static xbt::signal<void(DiskAction const&, Action::State, Action::State)> on_state_change;

  DiskAction(Model* model, double cost, bool failed, DiskImpl* disk, s4u::Io::OpType type)
      : Action(model, cost, failed), type_(type), disk_(disk){};

  /**
   * @brief diskAction constructor
   *
   * @param model The DiskModel associated to this DiskAction
   * @param cost The cost of this DiskAction in bytes
   * @param failed [description]
   * @param var The lmm variable associated to this DiskAction if it is part of a LMM component
   * @param disk The Disk associated to this DiskAction
   * @param type [description]
   */
  DiskAction(kernel::resource::Model* model, double cost, bool failed, kernel::lmm::Variable* var, DiskImpl* disk,
             s4u::Io::OpType type)
      : Action(model, cost, failed, var), type_(type), disk_(disk){};

  void set_state(simgrid::kernel::resource::Action::State state) override;

  s4u::Io::OpType type_;
  DiskImpl* disk_;
};

} // namespace resource
} // namespace kernel
} // namespace simgrid
#endif /* DISK_INTERFACE_HPP_ */
