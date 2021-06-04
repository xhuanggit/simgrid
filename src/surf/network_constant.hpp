/* Copyright (c) 2013-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef NETWORK_CONSTANT_HPP_
#define NETWORK_CONSTANT_HPP_

#include "network_interface.hpp"

namespace simgrid {
namespace kernel {
namespace resource {

class NetworkConstantModel : public NetworkModel {
public:
  using NetworkModel::NetworkModel;
  Action* communicate(s4u::Host* src, s4u::Host* dst, double size, double rate) override;
  double next_occurring_event(double now) override;
  void update_actions_state(double now, double delta) override;

  LinkImpl* create_link(const std::string& name, const std::vector<double>& bws) override;
  LinkImpl* create_wifi_link(const std::string& name, const std::vector<double>& bws) override;
};

class NetworkConstantAction final : public NetworkAction {
public:
  NetworkConstantAction(NetworkConstantModel* model_, s4u::Host& src, s4u::Host& dst, double size);
  void update_remains_lazy(double now) override;
};

} // namespace resource
} // namespace kernel
} // namespace simgrid

#endif /* NETWORK_CONSTANT_HPP_ */
