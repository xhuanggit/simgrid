/* Copyright (c) 2017-2021. The SimGrid Team. All rights reserved.               */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "catch.hpp"

#include "NetZone_test.hpp" // CreateHost callback
#include "simgrid/kernel/routing/FatTreeZone.hpp"
#include "simgrid/s4u/Engine.hpp"

TEST_CASE("kernel::routing::FatTreeZone: Creating Zone", "")
{
  simgrid::s4u::Engine e("test");
  simgrid::s4u::ClusterCallbacks callbacks(CreateHost{});
  REQUIRE(create_fatTree_zone("test", e.get_netzone_root(), {2, {4, 4}, {1, 2}, {1, 2}}, callbacks, 1e9, 10,
                              simgrid::s4u::Link::SharingPolicy::SHARED));
}

TEST_CASE("kernel::routing::FatTreeZone: Invalid params", "")
{
  simgrid::s4u::Engine e("test");
  simgrid::s4u::ClusterCallbacks callbacks(CreateHost{});

  SECTION("0 levels")
  {
    REQUIRE_THROWS_AS(create_fatTree_zone("test", e.get_netzone_root(), {0, {4, 4}, {1, 2}, {1, 2}}, callbacks, 1e9, 10,
                                          simgrid::s4u::Link::SharingPolicy::SHARED),
                      std::invalid_argument);
  }

  SECTION("Invalid down links")
  {
    REQUIRE_THROWS_AS(create_fatTree_zone("test", e.get_netzone_root(), {2, {4}, {1, 2}, {1, 2}}, callbacks, 1e9, 10,
                                          simgrid::s4u::Link::SharingPolicy::SHARED),
                      std::invalid_argument);
  }

  SECTION("Invalid up links")
  {
    REQUIRE_THROWS_AS(create_fatTree_zone("test", e.get_netzone_root(), {2, {4, 4}, {1}, {1, 2}}, callbacks, 1e9, 10,
                                          simgrid::s4u::Link::SharingPolicy::SHARED),
                      std::invalid_argument);
  }

  SECTION("Invalid link count")
  {
    REQUIRE_THROWS_AS(create_fatTree_zone("test", e.get_netzone_root(), {2, {4, 4}, {1, 2}, {1}}, callbacks, 1e9, 10,
                                          simgrid::s4u::Link::SharingPolicy::SHARED),
                      std::invalid_argument);
  }

  SECTION("Down links with zeroes")
  {
    REQUIRE_THROWS_AS(create_fatTree_zone("test", e.get_netzone_root(), {2, {4, 0}, {1, 2}, {1, 2}}, callbacks, 1e9, 10,
                                          simgrid::s4u::Link::SharingPolicy::SHARED),
                      std::invalid_argument);
  }

  SECTION("Up links with zeroes")
  {
    REQUIRE_THROWS_AS(create_fatTree_zone("test", e.get_netzone_root(), {2, {4, 4}, {0, 2}, {1, 2}}, callbacks, 1e9, 10,
                                          simgrid::s4u::Link::SharingPolicy::SHARED),
                      std::invalid_argument);
  }

  SECTION("Link count with zeroes")
  {
    REQUIRE_THROWS_AS(create_fatTree_zone("test", e.get_netzone_root(), {2, {4, 4}, {1, 2}, {1, 0}}, callbacks, 1e9, 10,
                                          simgrid::s4u::Link::SharingPolicy::SHARED),
                      std::invalid_argument);
  }

  SECTION("0 bandwidth")
  {
    REQUIRE_THROWS_AS(create_fatTree_zone("test", e.get_netzone_root(), {2, {4, 4}, {1, 2}, {1, 2}}, callbacks, 0, 10,
                                          simgrid::s4u::Link::SharingPolicy::SHARED),
                      std::invalid_argument);
  }

  SECTION("Negative latency")
  {
    REQUIRE_THROWS_AS(create_fatTree_zone("test", e.get_netzone_root(), {2, {4, 4}, {1, 2}, {1, 2}}, callbacks, 1e9,
                                          -10, simgrid::s4u::Link::SharingPolicy::SHARED),
                      std::invalid_argument);
  }
}
