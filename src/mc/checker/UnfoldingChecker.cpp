#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "UnfoldingChecker.hpp"
#include "xbt/ex.h"

namespace simgrid {
namespace mc {

class UnfoldingEvent;
bool isUnionConfig(EventSet hist1, UnfoldingEvent evt1, EventSet hist2, UnfoldingEvent evt2);

class Transition { // this class can be extended from SimGridMC
public:
  bool enabled;
  int id;
  bool isDependent(Transition* other)
  {
    // TODO
    return true;
  }
};

bool EventSet::contains(UnfoldingEvent e)
{
  return events_.find(e) != events_.end();
}

/** @brief Check if I'm dependent with another EventSet */
bool EventSet::depends(EventSet s2)
{
  for (auto e1 : events_)
    for (auto e2 : s2.events_)
      if (e1.transition->isDependent(e2.transition))
        return true;
  return false;
}

bool EventSet::isConfig()
{
  // checking conflict relation between one event and all other events in the set
  for (auto e1 : events_)
    for (auto e2 : events_) {
      if (e1 == e2)
        continue;

      if (e1.isConflict(e2))
        return false;

      // Every event of the history should be in the set
      for (auto ancestor : e2.getHistory())
        if (not this->contains(ancestor))
          return false;
    }
  return true;
}
EventSet EventSet::makeUnion(EventSet s1, EventSet s2)
{
  EventSet res;
  res.events_.insert(s1.events_.begin(), s1.events_.end());
  res.events_.insert(s2.events_.begin(), s2.events_.end());
  return res;
}
EventSet EventSet::makeIntersection(EventSet s1, EventSet s2)
{
  EventSet res;
  std::set_intersection(s1.events_.begin(), s1.events_.end(), s2.events_.begin(), s2.events_.end(),
                        std::back_inserter(res.events_));
  return res;
}
const UnfoldingEvent* EventSet::first()
{
  return &*events_.begin();
}

size_t EventSet::size() const
{
  return events_.size();
}
bool EventSet::empty() const noexcept
{
  return this->events_.empty();
}
std::set<UnfoldingEvent>::const_iterator EventSet::begin() const
{
  return events_.begin();
}
std::set<UnfoldingEvent>::const_iterator EventSet::end() const
{
  return events_.end();
}
std::set<UnfoldingEvent>::iterator EventSet::begin()
{
  return events_.begin();
}
std::set<UnfoldingEvent>::iterator EventSet::end()
{
  return events_.end();
}
bool EventSet::operator==(const EventSet& other) const
{
  return this->events_ == other.events_;
}
void EventSet::insert(UnfoldingEvent e)
{
  events_.insert(e);
}

UnfoldingEvent::UnfoldingEvent(Transition* t, EventSet causes)
{
  static int nb_events = 0; // To generate the unique identifier

  this->id         = nb_events++;
  this->causes     = causes;
  this->transition = t;
}

// Recursively compute the history of a given event by adding the causes of all ancestors
EventSet UnfoldingEvent::getHistory() const
{
  if (causes.size() == 0) // No direct ancestor => empty history
    return causes;
  else {
    EventSet res = causes;
    for (auto ancestor : causes) {
      EventSet h1 = ancestor.getHistory();
      res.events_.insert(h1.begin(), h1.end());
    }
    return res;
  }
}

/** @brief check for conflict in the history of current and provided events
 *
 * In the paper, this.isConflict(other) is written "this # other"
 */
bool UnfoldingEvent::isConflict(UnfoldingEvent otherEvent)
{

  /* if 2 event they have the same causes, just check their last transition */
  if (causes == otherEvent.causes)
    return transition->isDependent(otherEvent.transition);

  // if not, then check dependent relation on their full histories
  else
    return dependSetEvent(getHistory(), otherEvent.getHistory());
}

/** @brief Checks if current event is in immediate conflict with the provided one
 *
 * For that, there is two conditions to meet:
 *  - both events are in conflict (there is a conflict in their histories)
 *  -      Union(hist1,       hist2, evt2) is a valid configuration
 *    AND  Union(hist1, evt1, hist2)       is a valid configuration
 *
 * In the paper, e1.isImmediate(e2) will be written "e1 #ⁱ e2"
 */
bool UnfoldingEvent::isImmediateConflict(UnfoldingEvent evt2)
{
  // The first condition is easy to check
  if (!isConflict(evt2))
    return false;

  // Now, check the second condition
  EventSet hist1 = this->getHistory();
  EventSet hist2 = evt2.getHistory();

  // First compare the existing configurations
  for (auto e1 : hist1.events_)
    for (auto e2 : hist2.events_)
      if (e1.isConflict(e2))
        return false; // hist1 U hist2 is not a config => no immediate conflict
  // Compare the first config to the second new event
  for (auto e1 : hist1)
    if (e1.isConflict(evt2))
      return false;
  // Compare the second config to the first new event
  for (auto e2 : hist2)
    if (e2.isConflict(*this))
      return false;

  // Every tests passed
  return true;
}

// checking conflict relation between one event and one configuration or one history, it used when computing enC
// there is a better way by checking the event with maximal events in the configuration, (change type of enC )
bool UnfoldingEvent::conflictWithConfig(EventSet config)
{
  /* TODO: we don't really need to check the whole config. The maximal event should be enough.
   * But since that's a set and not a vector, we don't know which is the maximal event anymore...
   */
  for (auto evt : config)
    if (isConflict(evt))
      return true;
  return false;
}

void UnfoldingEvent::getEnabledTransition(std::set<Transition*>* whereto)
{
  THROW_UNIMPLEMENTED; // FIXME
}

// this operator is used for ordering in a set (need a key)
bool UnfoldingEvent::operator<(const UnfoldingEvent& other) const
{
  return id < other.id;
}

/** @brief check semantic equality (same transition, same causes) */
bool UnfoldingEvent::operator==(const UnfoldingEvent& other) const
{
  return ((transition->id == other.transition->id) && (causes == other.causes));
}

void Configuration::getEnabledTransition(std::set<Transition*>* whereto)
{
  for (auto e : maxEvent)
    e.getEnabledTransition(whereto);
}
EventSet Configuration::generateEvents(Transition* t)
{
  EventSet res;
  /* FIXME: Generate one event per transition and per subset of maxEvent */
  THROW_UNIMPLEMENTED;
  return res;
}

/* for each event in C, search all enabled transition in the state of that event
   then creating new event based on enabled transition and configuration C*/
void UnfoldingChecker::extend(Configuration C, EventSet& enC)
{

  std::set<Transition*> enabledTransitions;
  C.getEnabledTransition(&enabledTransitions);
  for (auto enabledT : C) {
    for (auto trans : enabledTransitions) {
      for (auto newEvent : C.generateEvents(trans)) {
        U.insert(newEvent);
        if (not newEvent.conflictWithConfig(C))
          enC.insert(newEvent);
      }
    }
  }
}
/** @brief Make the application reach that event state
 *  Precondition: the application is in a direct cause of e
 */
void UnfoldingEvent::execute()
{
  UnfoldingChecker::getSession().execute(this->transition);
  this->appState = std::unique_ptr<simgrid::mc::State>(new simgrid::mc::State(this->id));
}

void UnfoldingChecker::explore(Configuration C, EventSet D, EventSet A)
{
  UnfoldingEvent* e;
  EventSet enC;
  extend(C, enC);
  if (enC.empty())
    return;

  if (A.empty())
    e = enC.begin();
  else {
    // if A is not empty, chose one event in the intersection of A and enC
    EventSet intersection;
    set_intersection(enC.begin(), enC.end(), A.begin(), A.end(), std::back_inserter(intersection));
    e = *intersection.begin();
  }
  // execute event e
  e.execute(); // this function must be defined
  UnfoldingEvent* newEvent = e + e.transition;
  EventSet A1, C1;
  A1 = A.erase(e);
  C1 = C.insert(e);
  explore(C1, D, A1);
  EventSet J, tempJ, dif; // dif = J-C

  J = computeAlt(J, C, D, tempJ, U);
  std::set_difference(J.begin(), J.end(), C.begin(), C.end(), std::inserter(dif, dif.end()));
  if (!J.empty())
    explore(C, D.insert(e), dif);
  remove(e, C, D); // this function must be defined
}

EventSet generateEvent(EventSet C, Transition t, UnfoldingEvent evtFather, EventSet& enC)
{
  /* we only create new event if its transition depend on the transition of its ancestor,
   and it has not been created before (existing in U)  */

  EventSet result, dif;
  EventSet causesCandidate;
  causesCandidate.insert(evtFather);
  UnfoldingEvent e;
  std::EventSet::iterator rit;

  /* calculating all history candidates for transition t
   * putting all event which is independent with ancestor into causesCandidates
   */

  if (t.isDependent(evtFather.transition)) {
    for (auto e : C)
      if (!evtFather.transition.isDependent(e.transition)) {
        // checking 2 events have the same ancestor?
        if (evtFather.causes == e.causes)
          causesCandidate.insert(e);
      }
    // generating events from causesCandidate by iterating subset of causesCandidates to have new events
    EventSet causes;
    causes.insert(evtFather);
    UnfoldingEvent e = new UnfoldingEvent(nbe, t, causes);
    nbe++;
    EventSet temp;
    causesCandidate.erase(evtFather);
    result = generateFromCandidate(result, t, evtFather, temp, causesCandidate);
    for (auto it : result)
      // checking whether e already created? if not add to e to U
      if (eventInSet(it, U)) {
        it.id = nbe;
        nbe++;
        U.insert(it);

        // this following condition need add one more condition causes of e belong to C

        if (!e.conflictWithConfig(enC))
          enC.insert(e);
      }
  } else
    break;

  return result;
}

/* this function is used to generate events from a candidate Causes
 */
void generateFromCandidate(EventSet& result, Transition t, UnfoldingEvent evtFather, EventSet& temp,
                           EventSet& candidateCauses)
{
  UnfoldingEvent e;

  if (candidateCauses.empty()) {
    EventSet& temp1 = temp;
    temp1.insert(evtFather);
    e = new UnfoldingEvent(0, t, temp1); // need update id
    result.insert(e);
    return;
  } else {
    EventSet::iterator it = candidateCauses.begin();
    UnfoldingEvent a      = *it;
    EventSet t1, t2, t3, t4;
    t1 = temp;
    t3 = candidateCauses;
    t1.insert(a);
    t3.erase(a);
    generateFromCandidate(result, t, evtFather, t1, t3);
    t2 = temp;
    t2 = t4;
    t4.erase(a);
    generateFromCandidate(result, t, evtFather, t2, t4);
  }
}

void computeAlt(EventSet J, EventSet C, EventSet D, EventSet tempJ, EventSet U)
{
  if (!J.empty())
    return;
  EventSet t = C;
  t.insert(tempJ.begin(), tempJ.end());
  // checking C v U is a configuration?
  if (isConfig(t)) {
    int count = 0;
    for (auto it : D)
      for (auto it1 : t)
        if (it.isImmediateConflict(it1) && eventInSet(it1, U)) {
          count++;
          break;
        }
    if (D.size() == count) {
      J = tempJ;
      return;
    }
  }

  if (U.empty()) {

    return;
  } else {
    EventSet::iterator it = U.begin();
    UnfoldingEvent a      = *it;
    EventSet t1, t2, t3, t4;
    t1 = tempJ;
    t3 = U;
    t1.insert(a);
    t3.erase(a);
    computeAlt(J, C, D, t1, t3);
    t2 = tempJ;
    t4 = U;
    t4.erase(a);
    computeAlt(J, C, D, t2, t4);
  }
}

void remove(UnfoldingEvent e, EventSet C, EventSet& D)
{
  EventSet unionSet, res, res1;
  unionSet = C;
  unionSet.insert(D.begin(), D.end());

  for (auto e1 : U)
    for (auto e2 : unionSet) {
      if (e1.isImmediateConflict(e2))
        res.insert(e1);
      break;
    }

  for (auto e1 : res) {
    EventSet h = e1.getHistory();
    res        = res.insert(h.begin(), h.end());
  }

  res = res.insert(unionSet.begin(), unionSet.end());

  // move e from U to G if the condition is satisfied
  if (not eventInSet(e, res)) {
    U.erase(e);
    G.insert(e);
  }

  // move history of ê from U to G
  for (auto e1 : U)
    if (e1.isImmediateConflict(e)) {
      EventSet h = e1.getHistory();
      h.insert(e1);
      for (auto e2 : h)
        if (not eventInSet(e2, res)) {
          U.erase(e2);
          G.insert(e2);
        }
    }
}

Checker* createUnfoldingChecker(Session& session)
{
  return new UnfoldingChecker(session);
}
}
}
