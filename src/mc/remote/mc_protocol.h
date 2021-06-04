/* Copyright (c) 2015-2021. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SIMGRID_MC_PROTOCOL_H
#define SIMGRID_MC_PROTOCOL_H

// ***** Environment variables for passing context to the model-checked process

/** Environment variable name used to pass the communication socket.
 *
 * It is set by `simgrid-mc` to enable MC support in the children processes
 */
#define MC_ENV_SOCKET_FD "SIMGRID_MC_SOCKET_FD"

#ifdef __cplusplus

#include "mc/datatypes.h"
#include "simgrid/forward.h" // aid_t
#include <array>
#include <cstdint>
#include <xbt/dynar.h>
#include <xbt/mmalloc.h>
#include <xbt/utility.hpp>

// ***** Messages
namespace simgrid {
namespace mc {

XBT_DECLARE_ENUM_CLASS(MessageType, NONE, INITIAL_ADDRESSES, CONTINUE, IGNORE_HEAP, UNIGNORE_HEAP, IGNORE_MEMORY,
                       STACK_REGION, REGISTER_SYMBOL, DEADLOCK_CHECK, DEADLOCK_CHECK_REPLY, WAITING, SIMCALL_HANDLE,
                       SIMCALL_IS_VISIBLE, SIMCALL_IS_VISIBLE_ANSWER, SIMCALL_TO_STRING, SIMCALL_TO_STRING_ANSWER,
                       SIMCALL_DOT_LABEL, ASSERTION_FAILED, ACTOR_ENABLED, ACTOR_ENABLED_REPLY, FINALIZE);

} // namespace mc
} // namespace simgrid

constexpr unsigned MC_MESSAGE_LENGTH = 512;

/** Basic structure for a MC message
 *
 *  The current version of the client/server protocol sends C structures over `AF_LOCAL`
 *  `SOCK_SEQPACKET` sockets. This means that the protocol is ABI/architecture specific:
 *  we currently can't model-check a x86 process from a x86_64 process.
 *
 *  Moreover the protocol is not stable. The same version of the library should be used
 *  for the client and the server.
 */

/* Basic structure: all message start with a message type */
struct s_mc_message_t {
  simgrid::mc::MessageType type;
};

struct s_mc_message_int_t {
  simgrid::mc::MessageType type;
  uint64_t value;
};

/* Client->Server */
struct s_mc_message_initial_addresses_t {
  simgrid::mc::MessageType type;
  xbt_mheap_t mmalloc_default_mdp;
  unsigned long* maxpid;
  xbt_dynar_t actors;
  xbt_dynar_t dead_actors;
};

struct s_mc_message_ignore_heap_t {
  simgrid::mc::MessageType type;
  int block;
  int fragment;
  void* address;
  size_t size;
};

struct s_mc_message_ignore_memory_t {
  simgrid::mc::MessageType type;
  uint64_t addr;
  size_t size;
};

struct s_mc_message_stack_region_t {
  simgrid::mc::MessageType type;
  s_stack_region_t stack_region;
};

struct s_mc_message_register_symbol_t {
  simgrid::mc::MessageType type;
  std::array<char, 128> name;
  int (*callback)(void*);
  void* data;
};

/* Server -> client */
struct s_mc_message_simcall_handle_t {
  simgrid::mc::MessageType type;
  aid_t aid_;
  int times_considered_;
};

struct s_mc_message_restore_t {
  simgrid::mc::MessageType type;
  int index;
};

struct s_mc_message_actor_enabled_t {
  simgrid::mc::MessageType type;
  aid_t aid; // actor ID
};

/* RPC */
struct s_mc_message_simcall_is_visible_t { // MessageType::SIMCALL_IS_VISIBLE
  simgrid::mc::MessageType type;
  aid_t aid;
};
struct s_mc_message_simcall_is_visible_answer_t { // MessageType::SIMCALL_IS_VISIBLE_ANSWER
  simgrid::mc::MessageType type;
  bool value;
};

struct s_mc_message_simcall_to_string_t { // MessageType::SIMCALL_TO_STRING or MessageType::SIMCALL_DOT_LABEL
  simgrid::mc::MessageType type;
  aid_t aid;
  int time_considered;
};
struct s_mc_message_simcall_to_string_answer_t { // MessageType::SIMCALL_TO_STRING_ANSWER
  simgrid::mc::MessageType type;
  char value[1024];
};

#endif // __cplusplus
#endif
