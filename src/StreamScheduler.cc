/**
 * \file StreamScheduler.cc
 *
 * This file contains the default implementation of the StreamScheduler interface.
 *
 * This implementation of the StreamScheduler implements simple round robin 
 * assignment of streams to resource blocks.
 */

#include "MessageTypes.h"
#include "ScheduleInfo_m.h"
#include "StreamScheduler.h"
#include "TransReqList_m.h"
#include <iostream>
#include <memory>
#include <stdexcept>

using std::vector;
using std::unordered_map;

Define_Module(StreamScheduler);

void StreamScheduler::initialize(){
	this->initOffset = par("initOffset");
	this->downRB = par("downResourceBlocks");
	this->upRB = par("upResourceBlocks");
	this->upStatic = par("upStatic");
	this->downStatic = par("downStatic");
	this->numberOfMs = par("numberOfMobileStations");
	this->infos = vector<StreamInfo*>();
	this->streamSchedPeriod = par("streamSchedPeriod");
	this->tti = par("tti");
	this->epsilon = par("epsilon");
	this->sinrEstimate.resize(numberOfMs,nullptr);
	this->longtermSinrEstimate.resize(numberOfMs,nullptr);
	this->estimateBS = nullptr;
	this->longtermEstimateBS = nullptr;
	this->staticSchedLength = par("staticSchedLength");

	// Set the default packet sort criterion to creation time
	this->defaultPacketSorter = [](const KoiData *left,const KoiData *right) 
		-> bool{ return left->getCreationTime()<right->getCreationTime(); };

	if(!upStatic || !downStatic){
		// Produce the first schedule right at the init offset
		// We will need to make certain that all mobile stations have reported 
		// their streams at the time the first schedule is computed.
		scheduleAt(simTime()+initOffset,
				new cMessage("",MessageType::scheduleStreams));
		scheduleAt(simTime()+initOffset,
				new cMessage("",MessageType::scheduleRBs));
	}
	if(upStatic || downStatic){
		scheduleAt(simTime()+initOffset-epsilon,
				new cMessage("",MessageType::genStaticSchedule));
	}
	scheduleAt(simTime()+epsilon,
			new cMessage("",MessageType::scheduleInfo));
}

void StreamScheduler::scheduleDynStreams(){
	if(!this->infos.empty()){
		// First, clear the current assignment
		this->rbAssignments.clear();
		// This simple algorithm assigns streams to resource blocks 
		// by a round robin.
		int rbUp = 0;
		int rbDown = 0;
		for(auto info:this->infos){
			if(info->getD2d()){
				// If the stream is a D2D link, the scheduler 
				// does not only need to decide which ressource 
				// block the stream should use, but also 
				// which frequency band half, the one assigned 
				// to down or up traffic.
				//
				// We don't need assignments for up/down here 
				// as D2D streams only have one direction, 
				// directly from the sending MS to the receiving
				// MS.
				if(rbUp/(double)upRB<=rbDown/(double)downRB){
					this->rbAssignments[info->getStreamId()][MessageDirection::d2d] 
						= std::make_pair<MessageDirection,int>(MessageDirection::up,rbUp%upRB);
					rbUp++;
				}
				else{
					this->rbAssignments[info->getStreamId()][MessageDirection::d2d] 
						= std::make_pair<MessageDirection,int>(MessageDirection::down,rbDown%downRB);
					rbDown++;
				}
			}
			else{
				this->rbAssignments[info->getStreamId()][MessageDirection::up] 
					= std::make_pair<MessageDirection,int>(MessageDirection::up,rbUp%upRB);
				this->rbAssignments[info->getStreamId()][MessageDirection::down] 
					= std::make_pair<MessageDirection,int>(MessageDirection::down,rbDown%downRB);
				rbUp++;
				rbDown++;
			}
			delete info;
		}
		// Remove all stream infos
		this->infos.clear();
	}
}

std::unordered_map<int,ScheduleList> StreamScheduler::scheduleStatStreams(){
	throw std::runtime_error("Static Scheduling not implemented for this Stream Scheduler.");
}

void StreamScheduler::distributeStaticSchedules(){
	std::unordered_map<int,ScheduleList> schedules(scheduleStatStreams());
	for(auto& origin:schedules){
		StaticSchedule* sched = new StaticSchedule();
		sched->setScheduleLength(staticSchedLength);
		sched->setOrigin(origin.first);
		sched->setSchedule(origin.second);
		if(origin.first == -1){
			// Send schedule to local BS
			send(sched,"toBs");
		}
		else{
			// Send to MS
			send(sched,"toMs",origin.first);
		}
	}
}

void StreamScheduler::handleMessage(cMessage *msg){
  switch(msg->getKind()){
    case MessageType::streamInfo:{
        StreamInfo *tmp = dynamic_cast<StreamInfo*>(msg);
        this->infos.push_back(tmp);
        } break;
    case MessageType::streamTransReq:{
        StreamTransReq *req = dynamic_cast<StreamTransReq*>(msg);
        ResAssign& assignment = rbAssignments[req->getStreamId()][req->getMessageDirection()];
        this->requests[assignment.first][assignment.second].push_back(req);
        } break;
    case MessageType::scheduleStreams:{
				this->scheduleDynStreams();
        scheduleAt(simTime()+this->streamSchedPeriod,msg);
        } break;
    case MessageType::genStaticSchedule:{
				this->distributeStaticSchedules();
        } break;
		case MessageType::scheduleRBs:
				scheduledStations.clear();
        // Iterate over all message directions (up/down)
        for(auto iterDir=this->requests.begin();
            iterDir!=requests.end();
            ++iterDir){
          // Iterate over all Resource blocks in the current 
          // transmission direction
          for(auto iterRb=iterDir->second.begin();
              iterRb!=iterDir->second.end();
              ++iterRb){
            TransReqList *lst = new TransReqList();
            lst->setRequests(iterRb->second);
            // Gather SINR estimates for all 
            // MS with streams in this request
						std::shared_ptr<unordered_map<int,SINR*>> estimates = std::make_shared<unordered_map<int,SINR*>>();
            for(auto& req:iterRb->second){
              if(req->getMessageDirection()==MessageDirection::down){
                // Request from BS
                (*estimates)[-1] = estimateBS;
              }
              else{
                // Request from MS
                (*estimates)[req->getSrc()] = sinrEstimate[req->getSrc()];
              }
            }
            lst->setEstimates(estimates);
            // Clear all transmission requests for this 
            // stream from the stream schedulers queue.
						lst->setMessageDirection(iterDir->first);
            switch(iterDir->first){
              case MessageDirection::up:
                send(lst,"upRB$o",iterRb->first);
                break;
              case MessageDirection::down:
                send(lst,"downRB$o",iterRb->first);
                break;
            }
          }
        }
        scheduleAt(simTime()+tti,msg);
				scheduleAt(simTime()+epsilon,new cMessage("",MessageType::sendSchedules));
        break;
    case MessageType::streamSched:{
        StreamTransSched *sched = dynamic_cast<StreamTransSched*>(msg);
				if(scheduledStations.find(sched->getSrc())==scheduledStations.end()){
					scheduledStations.insert(sched->getSrc());
				}
				delete sched;
        } break;
		case MessageType::longTermSinrEst:
    case MessageType::sinrEst:{
        SINR *sinrEst = dynamic_cast<SINR*>(msg);
        this->handleSINREstimate(sinrEst);        
        } break;
		case MessageType::sendSchedules:{
				StreamTransSched* sched(nullptr);
				for(const int& id:scheduledStations){
					sched = new StreamTransSched();
					sched->setSrc(id);
					if(id==-1){
						// -1 is the local base station
						send(sched,"toBs");
					}
					else{
						// Send to the indicated mobile station
						send(sched,"toMs",id);
					}
				}
				// At this point, all requests for the current TTI have been handled, 
				// and the messages can be deleted.
        // Iterate over all message directions (up/down)
        for(auto iterDir=this->requests.begin();
            iterDir!=requests.end();
            ++iterDir){
          // Iterate over all Resource blocks in the current 
          // transmission direction
          for(auto iterRb=iterDir->second.begin();
              iterRb!=iterDir->second.end();
              ++iterRb){
            for(auto& req:iterRb->second){
							delete req;
            }
            iterRb->second.clear();
          }
        }
			} break;
		case MessageType::scheduleInfo:{
				if(msg->isSelfMessage()){
					// Self message generated in initialize, at this point any 
					// comparator from RB Schedulers will have arrived and the comparator
					// function will be set correctly.
					ScheduleInfo *inf = new ScheduleInfo();
					inf->setSortfn(defaultPacketSorter);
					inf->setUpStatic(upStatic);
					inf->setDownStatic(downStatic);
					// Forward the scheduling info to all local 
					// MS and the local BS
					send(inf->dup(),"toBs");
					for(int i=0; i<gateSize("toMs"); ++i){
						send(inf->dup(),"toMs",i);
					}
					delete inf;
				}
				else{
					// If this is not a self message, it must have arrived from the 
					// RB schedulers with a specific packet sorting criterion.
					ScheduleInfo *inf = dynamic_cast<ScheduleInfo*>(msg);
					defaultPacketSorter = inf->getSortfn();
				}
				delete msg;
			} break;
    default:
        std::cerr << "Received invalid Message in "
          << "StreamScheduler handleMessage"
          << std::endl;
  }
}

void StreamScheduler::printAssignment(){
  std::cout << "Stream" << "\t" << "Band" << "\t" << "RB" << std::endl;
  for(auto iterStream=this->rbAssignments.begin();
      iterStream!=this->rbAssignments.end();
      ++iterStream){
    for(auto iterDir=iterStream->second.begin();
        iterDir!=iterStream->second.end();
        ++iterDir){
      std::cout << iterStream->first << "\t"
        << iterDir->second.first << "\t"
        << iterDir->second.second << std::endl;
    }
  }
}

void StreamScheduler::handleSINREstimate(SINR *msg){
  if(msg->getMsId()==-1){
    // SINR estimate for the local BS
		if(msg->getKind()==MessageType::longTermSinrEst){
			if(longtermEstimateBS!=nullptr){
				delete longtermEstimateBS;
			}
			longtermEstimateBS = msg;
		}
		else{
			if(estimateBS!=nullptr){
				delete estimateBS;
			}
			estimateBS = msg;
		}
  }
  else{
    // SINR estimate for a local mobile station
		if(msg->getKind()==MessageType::longTermSinrEst){
			if(longtermSinrEstimate[msg->getMsId()]!=nullptr){
				delete longtermSinrEstimate[msg->getMsId()];
			}
			longtermSinrEstimate[msg->getMsId()] = msg;
		}
		else{
			if(sinrEstimate[msg->getMsId()]!=nullptr){
				delete sinrEstimate[msg->getMsId()];
			}
			sinrEstimate[msg->getMsId()] = msg;
		}
  }
}
