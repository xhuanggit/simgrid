/* Copyright (c) 2019-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "simgrid/s4u/Disk.hpp"
#include "simgrid/s4u/Engine.hpp"
#include "simgrid/s4u/Host.hpp"
#include "simgrid/s4u/Io.hpp"
#include "simgrid/simix.hpp"
#include "src/kernel/resource/DiskImpl.hpp"

namespace simgrid {

template class xbt::Extendable<s4u::Disk>;

namespace s4u {

xbt::signal<void(Disk&)> Disk::on_creation;
xbt::signal<void(Disk const&)> Disk::on_destruction;
xbt::signal<void(Disk const&)> Disk::on_state_change;

Disk* Disk::set_name(const std::string& name)
{
  name_ = name;
  return this;
}

Disk* Disk::set_read_bandwidth(double read_bw)
{
  kernel::actor::simcall([this, read_bw] { pimpl_->set_read_bandwidth(read_bw); });
  return this;
}

Disk* Disk::set_write_bandwidth(double write_bw)
{
  kernel::actor::simcall([this, write_bw] { pimpl_->set_write_bandwidth(write_bw); });
  return this;
}

double Disk::get_read_bandwidth() const
{
  return pimpl_->get_read_bandwidth();
}

double Disk::get_write_bandwidth() const
{
  return pimpl_->get_write_bandwidth();
}

Disk* Disk::set_host(Host* host)
{
  pimpl_->set_host(host);
  return this;
}

Host* Disk::get_host() const
{
  return pimpl_->get_host();
}

const std::unordered_map<std::string, std::string>* Disk::get_properties() const
{
  return pimpl_->get_properties();
}

const char* Disk::get_property(const std::string& key) const
{
  return pimpl_->get_property(key);
}

void Disk::set_property(const std::string& key, const std::string& value)
{
  kernel::actor::simcall([this, &key, &value] { this->pimpl_->set_property(key, value); });
}

IoPtr Disk::io_init(sg_size_t size, Io::OpType type)
{
  return Io::init()->set_disk(this)->set_size(size)->set_op_type(type);
}

IoPtr Disk::read_async(sg_size_t size)
{
  return IoPtr(io_init(size, Io::OpType::READ))->vetoable_start();
}

sg_size_t Disk::read(sg_size_t size)
{
  return IoPtr(io_init(size, Io::OpType::READ))->vetoable_start()->wait()->get_performed_ioops();
}

IoPtr Disk::write_async(sg_size_t size)
{
  return IoPtr(io_init(size, Io::OpType::WRITE)->vetoable_start());
}

sg_size_t Disk::write(sg_size_t size)
{
  return IoPtr(io_init(size, Io::OpType::WRITE))->vetoable_start()->wait()->get_performed_ioops();
}

void Disk::seal()
{
  kernel::actor::simcall([this]{ pimpl_->seal(); });
  get_host()->add_disk(this);
  Disk::on_creation(*this);
}
} // namespace s4u
} // namespace simgrid

/* **************************** Public C interface *************************** */

const char* sg_disk_get_name(const_sg_disk_t disk)
{
  return disk->get_cname();
}

sg_host_t sg_disk_get_host(const_sg_disk_t disk)
{
  return disk->get_host();
}

double sg_disk_read_bandwidth(const_sg_disk_t disk)
{
  return disk->get_read_bandwidth();
}

double sg_disk_write_bandwidth(const_sg_disk_t disk)
{
  return disk->get_write_bandwidth();
}

sg_size_t sg_disk_read(sg_disk_t disk, sg_size_t size)
{
  return disk->read(size);
}
sg_size_t sg_disk_write(sg_disk_t disk, sg_size_t size)
{
  return disk->write(size);
}

void* sg_disk_get_data(const_sg_disk_t disk)
{
  return disk->get_data();
}

void sg_disk_set_data(sg_disk_t disk, void* data)
{
  disk->set_data(data);
}
