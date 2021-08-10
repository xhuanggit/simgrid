/* Copyright (c) 2007-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "src/kernel/activity/SynchroRaw.hpp"
#include "simgrid/Exception.hpp"
#include "simgrid/kernel/resource/Action.hpp"
#include "src/kernel/activity/ConditionVariableImpl.hpp"
#include "src/kernel/activity/MutexImpl.hpp"
#include "src/kernel/activity/SemaphoreImpl.hpp"
#include "src/kernel/context/Context.hpp"
#include "src/surf/cpu_interface.hpp"
#include "src/surf/surf_interface.hpp"
#include <simgrid/s4u/Host.hpp>

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(simix_synchro, simix, "SIMIX Synchronization (mutex, semaphores and conditions)");

namespace simgrid {
namespace kernel {
namespace activity {

RawImpl& RawImpl::set_host(s4u::Host* host)
{
  host_ = host;
  return *this;
}
RawImpl& RawImpl::set_timeout(double timeout)
{
  timeout_ = timeout;
  return *this;
}

RawImpl* RawImpl::start()
{
  surf_action_ = host_->get_cpu()->sleep(timeout_);
  surf_action_->set_activity(this);
  return this;
}

void RawImpl::suspend()
{
  /* The suspension of raw synchros is delayed to when the actor is rescheduled. */
}

void RawImpl::resume()
{
  /* I cannot resume raw synchros directly. This is delayed to when the actor is rescheduled at
   * the end of the synchro. */
}

void RawImpl::cancel()
{
  /* I cannot cancel raw synchros directly. */
}

void RawImpl::post()
{
  if (surf_action_->get_state() == resource::Action::State::FAILED) {
    state_ = State::FAILED;
  } else if (surf_action_->get_state() == resource::Action::State::FINISHED) {
    state_ = State::SRC_TIMEOUT;
  }

  clean_action();
  /* Answer all simcalls associated with the synchro */
  finish();
}

void RawImpl::finish()
{
  XBT_DEBUG("RawImpl::finish() in state %s", to_c_str(state_));
  xbt_assert(simcalls_.size() == 1, "Unexpected number of simcalls waiting: %zu", simcalls_.size());
  smx_simcall_t simcall = simcalls_.front();
  simcalls_.pop_front();

  if (state_ == State::FAILED) {
    simcall->issuer_->context_->set_wannadie();
    simcall->issuer_->exception_ = std::make_exception_ptr(HostFailureException(XBT_THROW_POINT, "Host failed"));
  } else {
    xbt_assert(state_ == State::SRC_TIMEOUT, "Internal error in RawImpl::finish() unexpected synchro state %s",
               to_c_str(state_));
  }

  finish_callback_();
  simcall->issuer_->waiting_synchro_ = nullptr;
  simcall->issuer_->simcall_answer();
}

} // namespace activity
} // namespace kernel
} // namespace simgrid
