/**
 * @class EDFStreamScheduler EDFStreamScheduler.h
 *
 * This scheduler assigns streams to resource blocks by their deadline.
 *
 * The scheduler uses the EDFStreamSched algorithm described in S. Auroux, 
 * D. Parruca and H. Karl, "Joint real-time scheduling and interference 
 * coordination for wireless factory automation.". This algorithm tries to 
 * assigns streams to resource blocks in such a way that all deadlines are 
 * met and the minimum number of resource blocks is used. 
 *
 * The resource block assignment is recomputed every T seconds, where T is 
 * the least common multiple of the periods of all streams in the local cell.
 *
 * The actual scheduling of packets is computed on a per-timeslot basis by 
 * EDF scheduling for all streams assigned to each resource block. The 
 * assignment of resource blocks via the EDFStreamSched ensures that all 
 * RBEDFSchedulers can meet all deadlines for all streams assigned to them, 
 * as long as there are enough resource blocks to service all local 
 * communication streams.
 */

#pragma once

#include "StreamScheduler.h"
#include "MessageTypes.h"

#include <unordered_map>
#include <utility>
#include <vector>

struct EDFRb;
struct EDFStream{
  EDFStream(MessageDirection direction,
      double period,
      double deadline,
      unsigned long id,
      bool d2d):
    direction(direction),
    period(period),
    deadline(deadline),
    id(id),
    d2d(d2d){
  }
  inline bool operator==(const EDFStream& rhs){
    return (this->direction==rhs.direction) && (this->id==rhs.id);
  }
  MessageDirection direction;
  std::vector<EDFRb*> viableRbs;
  double period;
  double deadline;
  unsigned long id;
  bool d2d;
};

struct EDFRb{
  EDFRb(MessageDirection direction,
      int id):
    direction(direction),
    utilization(0.0),
    id(id){}
  inline bool operator==(const EDFRb& rhs){
    return (this->direction==rhs.direction) && (this->id==rhs.id);
  }
  MessageDirection direction;
  std::vector<EDFStream*> viableStreams;
  double utilization;
  int id;
};

class EDFStreamScheduler: public StreamScheduler{
  private:
    virtual void getViability(std::vector<EDFStream>& streams,
        std::vector<EDFRb>& blocks);
    virtual bool schedulabilityTest(const EDFRb& rb, const EDFStream& newStream);
    virtual double computeSchedulingInterval(const std::vector<EDFStream>& streams);

  protected:
    virtual void initialize();
    virtual void scheduleDynStreams();
};
