/* Copyright (c) 2009-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "simgrid/kernel/routing/NetPoint.hpp"
#include "simgrid/kernel/routing/VivaldiZone.hpp"
#include "simgrid/s4u/Engine.hpp"
#include "simgrid/s4u/Host.hpp"
#include "xbt/log.h"

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(surf_route, surf, "Routing part of surf");

namespace simgrid {

template class xbt::Extendable<kernel::routing::NetPoint>;

namespace kernel {
namespace routing {

simgrid::xbt::signal<void(NetPoint&)> NetPoint::on_creation;

NetPoint::NetPoint(const std::string& name, NetPoint::Type componentType) : name_(name), component_type_(componentType)
{
  simgrid::s4u::Engine::get_instance()->netpoint_register(this);
  simgrid::kernel::routing::NetPoint::on_creation(*this);
}

NetPoint* NetPoint::set_englobing_zone(NetZoneImpl* netzone_p)
{
  englobing_zone_ = netzone_p;
  if (netzone_p != nullptr)
    id_ = netzone_p->add_component(this);
  return this;
}

NetPoint* NetPoint::set_coordinates(const std::string& coords)
{
  if (not coords.empty())
    new vivaldi::Coords(this, coords);
  return this;
}
} // namespace routing
} // namespace kernel
} // namespace simgrid

/** @brief Retrieve a netpoint from its name
 *
 * Netpoints denote the location of host or routers in the network, to compute routes
 */
simgrid::kernel::routing::NetPoint* sg_netpoint_by_name_or_null(const char* name)
{
  return simgrid::s4u::Engine::get_instance()->netpoint_by_name_or_null(name);
}
