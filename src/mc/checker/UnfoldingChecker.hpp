#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <set>

namespace simgrid {
namespace mc {
class UnfoldingEvent;
class Transition;

class EventSet {
public:
  bool contains(UnfoldingEvent e);
  bool depends(EventSet s2);
  bool isConfig();
  static EventSet makeUnion(EventSet s1, EventSet s2);
  static EventSet makeIntersection(EventSet s1, EventSet s2);

  const UnfoldingEvent* first();

  size_t size() const;
  bool empty() const noexcept;
  std::set<UnfoldingEvent>::const_iterator begin() const;
  std::set<UnfoldingEvent>::const_iterator end() const;
  std::set<UnfoldingEvent>::iterator begin();
  std::set<UnfoldingEvent>::iterator end();
  void insert(UnfoldingEvent);

  bool operator==(const EventSet& other) const;

  std::set<UnfoldingEvent> events_;
};

class Configuration : public EventSet {
public:
  std::set<UnfoldingEvent> maxEvent; // Events recently added to events_

  void getEnabledTransition(std::set<Transition*>* whereto);
  EventSet generateEvents(Transition* t);
};

class UnfoldingEvent {
public:
  int id                       = -1;
  simgrid::mc::State* appState = nullptr;
  Transition* transition; // The last transition made to reach that state
  EventSet causes;        // used to store directed ancestors of event e

  UnfoldingEvent(Transition* t, EventSet causes);

  bool dependSetEvent(EventSet s1, EventSet s2);
  EventSet getHistory() const;
  bool isConflict(UnfoldingEvent otherEvent);
  bool isImmediateConflict(UnfoldingEvent other);
  bool conflictWithConfig(EventSet config);

  void getEnabledTransition(std::set<Transition*>* whereto);
  void UnfoldingEvent::execute();

  bool operator<(const UnfoldingEvent& other) const;
  bool operator==(const UnfoldingEvent& other) const;
};

class UnfoldingChecker {
  EventSet A, C, D, G, U;

  static Session& getSession();

private:
  void explore(EventSet C, EventSet D, EventSet A);
  void extend(Configuration C, EventSet& enC);

  static Session& session;
};
}
}
