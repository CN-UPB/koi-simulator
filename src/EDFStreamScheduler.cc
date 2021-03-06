/**
 * @file EDFStreamScheduler.cc 
 * Implementation of the EDFStreamScheduler class.
 *
 * For details on the scheduling algorithm, see S. Auroux, D. Parruca and 
 * H. Karl, "Joint real-time scheduling and interference coordination for 
 * wireless factory automation."
 */

#include "EDFStreamScheduler.h"
#include "includes.h"
#include "MessageTypes.h"
#include "StreamInfo_m.h"
#include "util.h"

#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

using std::vector;
using std::unordered_map;
using std::pair;

Define_Module(EDFStreamScheduler);

void EDFStreamScheduler::initialize(){
	StreamScheduler::initialize();
}


void EDFStreamScheduler::getViability(vector<EDFStream>& streams,
    vector<EDFRb>& blocks){
  // For now, the viability of a specific stream on a given resource block 
  // is generated randomly, as we do not yet have a process to quantify 
  // transmission quality per channel/resource block.
  for(EDFStream& s:streams){
    for(EDFRb& rb:blocks){
      if(((s.direction==rb.direction) || s.direction==MessageDirection::d2d) 
          && dblrand()>0.5){
        rb.viableStreams.push_back(&s);
        s.viableRbs.push_back(&rb);
      }
    }
  }
}

bool EDFStreamScheduler::schedulabilityTest(const EDFRb& rb,
    const EDFStream& newStream){
  // Currently, any packet needs only a single TTI to transmit, so all 
  // streams have an execution time e equal to the length of a tti.
	omnetpp::simtime_t e(tti);
  double newStreamLoad = e.dbl()/std::min(newStream.period,
      newStream.deadline);
  return rb.utilization+newStreamLoad <= 1.0;
}

double EDFStreamScheduler::computeSchedulingInterval(
    const std::vector<EDFStream>& streams){
  auto getPeriod = [](const EDFStream& s) -> double {return s.period;};
  vector<double> periods(streams.size());
  std::transform(streams.begin(),streams.end(),periods.begin(),getPeriod);
  return lcmSequence(periods);  
}

void EDFStreamScheduler::scheduleDynStreams(){
  if(!this->infos.empty()){
    // First, clear the current assignment
    this->rbAssignments.clear();
    // Convert StreamInfo messages to EDFStreams
    vector<EDFStream> streams;
    for(StreamInfo* inf:this->infos){
      if(inf->getD2d()){
        streams.emplace_back(MessageDirection::d2d,inf->getInterarrival(),
            inf->getDeadline(),inf->getStreamId(),inf->getD2d());
      }
      else{
        streams.emplace_back(MessageDirection::up,inf->getInterarrival(),
            inf->getDeadline(),inf->getStreamId(),inf->getD2d());
        streams.emplace_back(MessageDirection::down,inf->getInterarrival(),
            inf->getDeadline(),inf->getStreamId(),inf->getD2d());
      }
    }
    // Generate EDFRb instances for all UP and DOWN resource blocks
    vector<EDFRb> blocks;
    for(int i:assignedUpRB){
      blocks.emplace_back(MessageDirection::up,i);
    }
    for(int i:assignedDownRB){
      blocks.emplace_back(MessageDirection::down,i);
    }

    // Determine which streams are viable on which resource blocks
    getViability(streams,blocks);
    // Sort streams by the number of resource blocks they can be scheduled on
    auto compareStreamsRbNum = [](const EDFStream& first,
        const EDFStream& second) 
      -> bool {return first.viableRbs.size() < second.viableRbs.size();};
    std::sort(streams.begin(),streams.end(),compareStreamsRbNum);
    // Sort resource blocks by the number of streams which can be run on them
    auto compareRbsStreamNum = [](const EDFRb& first,
        const EDFRb& second) 
      -> bool {return first.viableStreams.size() > second.viableStreams.size();};
    std::sort(blocks.begin(),blocks.end(),compareRbsStreamNum);

		omnetpp::simtime_t e(tti);
    // Schedule the streams on the resource blocks
    bool scheduled;
    int unscheduledStreams = 0;
    int schedStreams = 0;
    for(EDFStream s:streams){
      scheduled = false;
      auto compareRbUtil = [](const EDFRb* first, const EDFRb* second)
        -> bool {return first->utilization > second->utilization;};
      std::sort(s.viableRbs.begin(),s.viableRbs.end(),compareRbUtil);
      for(EDFRb* rb:s.viableRbs){
        if(schedulabilityTest(*rb,s)){
          rb->utilization += e.dbl()/std::min(s.period,s.deadline);
          // Schedule streams s on resource block rb
          this->rbAssignments[s.id][s.d2d ? MessageDirection::d2d : s.direction]
            = std::make_pair<MessageDirection,int>((MessageDirection&&)rb->direction,(int&&)rb->id);
          scheduled = true;
          schedStreams++;
          break;
        }
      }
      if(!scheduled){
        ++unscheduledStreams;
      }
    }
    // Set validity period for schedule to least common multiple of the 
    // periods of all streams.
    this->streamSchedPeriod = computeSchedulingInterval(streams);
  }
}

