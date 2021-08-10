/* Copyright (c) 2013-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SURF_ROUTING_FLOYD_HPP_
#define SURF_ROUTING_FLOYD_HPP_

#include <simgrid/kernel/routing/RoutedZone.hpp>

namespace simgrid {
namespace kernel {
namespace routing {

/** @ingroup ROUTING_API
 *  @brief NetZone with an explicit routing computed at initialization with Floyd-Warshal
 *
 *  The path between components is computed at creation time from every one-hop links,
 *  using the Floyd-Warshal algorithm.
 *
 *  This result in rather small platform file, slow initialization time,  and intermediate memory requirements
 *  (somewhere between the one of @{DijkstraZone} and the one of @{FullZone}).
 */
class XBT_PRIVATE FloydZone : public RoutedZone {
  /* vars to compute the Floyd algorithm. */
  std::vector<std::vector<int>> predecessor_table_;
  std::vector<std::vector<unsigned long>> cost_table_;
  std::vector<std::vector<std::unique_ptr<Route>>> link_table_;

  void init_tables(unsigned int table_size);
  void do_seal() override;

public:
  using RoutedZone::RoutedZone;
  FloydZone(const FloydZone&) = delete;
  FloydZone& operator=(const FloydZone&) = delete;

  void get_local_route(const NetPoint* src, const NetPoint* dst, Route* into, double* latency) override;
  void add_route(NetPoint* src, NetPoint* dst, NetPoint* gw_src, NetPoint* gw_dst,
                 const std::vector<s4u::LinkInRoute>& link_list, bool symmetrical) override;
};
} // namespace routing
} // namespace kernel
} // namespace simgrid

#endif /* SURF_ROUTING_FLOYD_HPP_ */
