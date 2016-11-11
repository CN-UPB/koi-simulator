/*
 * MsMac.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#include "MsMac.h"
#include "PositionExchange_m.h"
#include "KoiData_m.h"
#include "QueueSort_m.h"
#include "TransmitRequest_m.h"
#include "Schedule_m.h"
#include "StreamInfo_m.h"
#include "StreamTransReq_m.h"
#include "StreamTransSched_m.h"
#include "SINR_m.h"
#include "TransInfo_m.h"
#include "util.h"
#include <algorithm>
#include <stdlib.h>
#include <cmath>
#include <random>
#include <set>
#include <utility>

using std::set;

Define_Module(MsMac);

inline simtime_t MsMac::positionResendTime()  {
    if(simTime() == initOffset)  {
        return simTime() + (positionResendInterval * tti) - (2 * epsilon); //send with the schedules; one sync point
    }
    else  {
        return simTime() + positionResendInterval * tti;
    }
}

/* calc the init pos of the ms in the cell */
Position MsMac::initMsPosition(int quadrant, double alpha, double beta, double gamma)  { //for random movement of MS
	Position msPos;
	
	double angle1;
	double angle2;
	Position tri1;
	Position tri2;
	
	velocity.push_back(0);
	velocity.push_back(0);

	//triangles of the cell
	/*
	if(quadrant == 0)  {
	    angle1 = 330;
	    angle2 = 30;
	}
	else if(quadrant == 1)  {
		angle1 = 30;
		angle2 = 90;
	}
	else if(quadrant == 2)  {
		angle1 = 90;
		angle2 = 150;
	}
	else if(quadrant == 3)  {
		angle1 = 150;
		angle2 = 210;
	}
	else if(quadrant == 4)  {
		angle1 = 210;
		angle2 = 270;
	}
	else if(quadrant == 5)  {
		angle1 = 270;
		angle2 = 330;
	}

	//get the point of the triangle
	tri1.x = cos(angle1/180*M_PI) * radius + initBsPos.x;
	tri1.y = sin(angle1/180*M_PI) * radius + initBsPos.y;

	tri2.x = cos(angle2/180*M_PI) * radius + initBsPos.x;
	tri2.y = sin(angle2/180*M_PI) * radius + initBsPos.y;


	//calc the points on the triangle (barycentric coordinates)
	msPos.x = (alpha * tri1.x + beta * tri2.x + gamma * initBsPos.x) / (alpha + beta + gamma);
	msPos.y = (alpha * tri1.y + beta * tri2.y + gamma * initBsPos.y) / (alpha + beta + gamma);
	
	return msPos;
	*/
	
	double cell0x = 500.0;
	double cell0y = 1000.0;
	double cell1x = 500.0;
	double cell1y = 0.0;
	
	double cell2x = 933.0;
	double cell2y = 750.0;
	double cell3x = 500.0;
	double cell3y = 500.0;
	double cell4x = 67.0;
	double cell4y = 750.0;

	double cell5x = 67.0;
	double cell5y = 250.0;
	double cell6x = 933.0;
	double cell6y = 250.0;
	
	ofstream PL;
	PL.open ("BS_" + std::to_string( (long long) bsId ) + "_Pathloss.txt", fstream::app);

	// Use c++ random number to allow different position seeds easily
	static int seedinc = 0;
	int PosSeed = par("PosSeed");
	std::mt19937 gen(PosSeed + seedinc);
	do{
		std::uniform_real_distribution<double> dis(-radius, radius);
		
		msPos.x = initBsPos.x + dis(gen);
		msPos.y = initBsPos.y + dis(gen);
		
		if(bsId == 3){
			double dist1 = sqrt(pow((msPos.x - initBsPos.x),2) + pow((msPos.y - initBsPos.y),2));
			
			double dist2 = sqrt(pow((msPos.x - cell0x),2) + pow((msPos.y - cell0y),2));
			double dist3 = sqrt(pow((msPos.x - cell1x),2) + pow((msPos.y - cell1y),2));
			double dist4 = sqrt(pow((msPos.x - cell2x),2) + pow((msPos.y - cell2y),2));
			double dist5 = sqrt(pow((msPos.x - cell4x),2) + pow((msPos.y - cell4y),2));
			double dist6 = sqrt(pow((msPos.x - cell5x),2) + pow((msPos.y - cell5y),2));
			double dist7 = sqrt(pow((msPos.x - cell6x),2) + pow((msPos.y - cell6y),2));
			if(dist2 < dist1 || dist3 < dist1 || dist4 < dist1 || dist5 < dist1 || dist6 < dist1 || dist7 < dist1){
				msPos.x = initBsPos.x;
				msPos.y = initBsPos.y;
				continue;
			}
		}
		
	}while(sqrt(pow((msPos.x - initBsPos.x),2) + pow((msPos.y - initBsPos.y),2)) > radius || sqrt(pow((msPos.x - initBsPos.x),2) + pow((msPos.y - initBsPos.y),2)) < 20.0);
	
	seedinc++;
	
	PL << msPos.x << " " << msPos.y << " " << 1/pow(10, (22*log10(sqrt(pow(msPos.x - 500.0,2) + pow(msPos.y - 500.0,2))) + 28 + 20*log10(2.5)) /10) << " " << 1/pow(10, (22*log10(sqrt(pow(msPos.x - 500.0,2) + pow(msPos.y - 1000.0,2))) + 28 + 20*log10(2.5)) /10) << " " << 1/pow(10, (22*log10(sqrt(pow(msPos.x - 500.0,2) + pow(msPos.y - 0.0,2))) + 28 + 20*log10(2.5)) /10) << " " << 1/pow(10, (22*log10(sqrt(pow(msPos.x - 933.0,2) + pow(msPos.y - 750.0,2))) + 28 + 20*log10(2.5)) /10) << " " << 1/pow(10, (22*log10(sqrt(pow(msPos.x - 67.0,2) + pow(msPos.y - 750.0,2))) + 28 + 20*log10(2.5)) /10) << " " << 1/pow(10, (22*log10(sqrt(pow(msPos.x - 67.0,2) + pow(msPos.y - 250.0,2))) + 28 + 20*log10(2.5)) /10) << " " << 1/pow(10, (22*log10(sqrt(pow(msPos.x - 933.0,2) + pow(msPos.y - 250.0,2))) + 28 + 20*log10(2.5)) /10) << std::endl;
	
	PL.close();
	return msPos;
} //end of MsMac::initMsPosition 

Position MsMac::initMsPositionLinear()  { //for MS Position along a straight road
	Position msPos;
	ofstream P_MS;
	P_MS.open ("MS_Positions.txt", fstream::app);

        velocity.push_back(30); //velocity in x-direction in m/s
	velocity.push_back(0); //velocity in y-direction in m/s

	msPos.x = 6;
        msPos.y = radius - 9; //for movement of the MS on a straight road
	P_MS << msPos.x << " " << msPos.y << "\n" << std::endl;
	
        P_MS.close();
	return msPos;
}

Position MsMac::initMsPositionRand(){
	double angle = uniform(0,360);
	// Choose radius between 10m and the cell radius
	// 10 meters is the minimum dist for the METIS model
	double distToBs = uniform(10,radius);
	// Convert angle to radians for math library functions
	angle = angle * (M_PI/180);
	Position initPos;
	initPos.x = initBsPos.x+distToBs*std::cos(angle);
	initPos.y = initBsPos.y+distToBs*std::sin(angle);
	return initPos;
}

Position MsMac::initMsPositionLine(){
	Position initPos;
	double msDists = (radius-1.01)/numberOfMobileStations;
	initPos.x = initBsPos.x+1.01+(msId*msDists);
	initPos.y = initBsPos.y;
	return initPos;
}

void MsMac::initialize()  {
	positionResendInterval = par("positionResendInterval");
	msId = par("msId");
	bsId = par("bsId");
	numberOfMobileStations = par("numberOfMobileStations");
	epsilon = par("epsilon");
	radius = par("radius");
	initBsPos.x = par("initBsXPos");
	initBsPos.y = par("initBsYPos");
	double initPosAlpha = par("initPosAlpha");
	double initPosBeta = par("initPosBeta");
	double initPosGamma = par("initPosGamma");
	int initQuadrant = par("initQuadrant");
	initOffset = par("initOffset");
	tti = par("tti");
	packetLength = par("packetLength");
	transmissionPower = par("transmissionPower");

	switch((int)par("positioning")){
		case MsMac::Placement::params:
			msPosition.x = par("initMsXPos");
			msPosition.y = par("initMsYPos");
			break;
		case MsMac::Placement::bySector:
			msPosition = initMsPosition(initQuadrant, 
					initPosAlpha, 
					initPosBeta, 
					initPosGamma); //for random MS positions in hexagonal cell
			break;
		case MsMac::Placement::linear:
			msPosition = initMsPositionLinear(); //for MS position along a straight road
			break;
		case MsMac::Placement::uniformRand:
			msPosition = initMsPositionRand();
			break;
		case MsMac::Placement::line:
			msPosition = initMsPositionLine();
			break;
		default:
			// TODO Notify the user that the value for MS positioning
			// is invalid.
			std::cout << "Invalid Ms placement algorithm " << (int) par("positioning") << std::endl;
	}

	std::string fname("rate-ms-"+std::to_string(bsId)+"-"+std::to_string(msId));
	rateFile = std::move(getResultFile(fname));
        
	//every tti send transmit requests to stream scheduler
	scheduleAt(simTime() + initOffset-epsilon, new cMessage("GEN_TRANSMIT_REQUEST"));

	// Send MS Position once at the very beginning for cluster generation
	PositionExchange *posEx = new PositionExchange("MS_POS_UPDATE");
	posEx->setId(msId);
	posEx->setPosition(msPosition);
	send(posEx, "toBsMac");
}

void MsMac::finish(){
	rateFile.close();
}

void MsMac::handleMessage(cMessage *msg)  {
	if(msg->getKind()==MessageType::streamSched)  {
		StreamTransSched *schedule = dynamic_cast<StreamTransSched*>(msg);
		std::pair<set<int>,set<int>> infos = std::make_pair(set<int>(),set<int>());

		if(schedule->getSrc() == msId)  {
			KoiData *currPacket = nullptr;
			int rate = 0;
			for(auto streamIter = streamQueues.begin();
					streamIter!=streamQueues.end(); ++streamIter){
				list<KoiData*>& currList = streamIter->second;
				for(auto packetIter = currList.begin(); 
						packetIter!=currList.end();){
					currPacket = *packetIter;
					if(currPacket->getScheduled()){
						packetIter = currList.erase(packetIter);
						currPacket->setTransPower(transmissionPower);
						// Set CQI for a fixed value until we decide on how to 
						// compute it
						currPacket->setCqi(15);
						this->take(currPacket);
						sendDelayed(currPacket, epsilon, "toPhy");
						// Store RB in frequency half dependent sets, so that we send 
						// at most 1 TransInfo for any given resource block, even 
						// if we transmit multiple packets using that block.
						if(currPacket->getMessageDirection()==MessageDirection::up
								|| currPacket->getMessageDirection()==MessageDirection::d2dUp){
							infos.first.insert(currPacket->getResourceBlock());
						}
						else{
							infos.second.insert(currPacket->getResourceBlock());
						}
						rate += currPacket->getBitLength();
					}
					else{
						++packetIter;
					}
				}
			}
		rateFile << rate << std::endl;
		}
		for(auto& rb:infos.first){
			TransInfo *info = new TransInfo();
			info->setBsId(bsId);
			info->setPower(transmissionPower);
			info->setRb(rb);
			info->setMsId(msId);
			info->setMessageDirection(MessageDirection::up);
			send(info,"toBsMac");
		}
		for(auto& rb:infos.second){
			TransInfo *info = new TransInfo();
			info->setBsId(bsId);
			info->setPower(transmissionPower);
			info->setRb(rb);
			info->setMsId(msId);
			info->setMessageDirection(MessageDirection::d2dDown);
			send(info,"toBsMac");
		}
		delete schedule;
	}
	else if(msg->getKind()==MessageType::transInfo){
		send(msg,"toPhy");
	}
	else if(msg->getKind()==MessageType::sortOrder){
		QueueSort *s = dynamic_cast<QueueSort*>(msg);
		comparator = s->getSortfn();
		delete msg;
	}
	else if(msg->isName("GEN_TRANSMIT_REQUEST"))  {
		// Send requests for each stream originating from this MS to the 
		// scheduler if that stream has packets.
		for(auto iter=this->streamQueues.begin(); iter!=streamQueues.end();
				++iter){
			if(!iter->second.empty()){
				StreamTransReq *req = new StreamTransReq();
				KoiData *queueHead = dynamic_cast<KoiData*>(iter->second.front());
				req->setSrc(this->msId);
				req->setDest(queueHead->getDest());
				req->setStreamId(iter->first);
				req->setPeriod(queueHead->getInterarrival());
				req->setPackets(&(iter->second));
				req->setRequestOrigin(msId);
				if(queueHead->getD2d()){
					req->setMessageDirection(MessageDirection::d2d);
				}
				else{
					req->setMessageDirection(MessageDirection::up);
				}
				send(req,"toScheduler");
			}
		}
		scheduleAt(simTime() + tti-epsilon, msg);
	}
	else if(msg->getKind()==MessageType::sinrEst)  {
		// Forward estimates to BS Mac
		send(msg,"toBsMac");
	}
	else if(msg->arrivedOn("fromApp"))  {
		// Packet arrived for sending from traffic generator
		switch(msg->getKind()){
			case MessageType::streamInfo:{
				// Add queue for the new stream
				StreamInfo *tmp = dynamic_cast<StreamInfo*>(msg);
				this->streamQueues[tmp->getStreamId()];
				send(tmp->dup(),"toScheduler");
				send(tmp->dup(),"toBsMac");
				delete msg;
			} break;
			case MessageType::koidata:{
				KoiData *data = dynamic_cast<KoiData*>(msg);
				list<KoiData*>& squeue(streamQueues[data->getStreamId()]);
				auto p = std::lower_bound(squeue.begin(),squeue.end(),data,comparator);
				this->streamQueues[data->getStreamId()].insert(p,data);
			} break;
		}
	}
	else if(msg->arrivedOn("fromPhy"))  {
		// Unpack the data bundle and forward data packets to app
		if(msg->getKind()==MessageType::koidata){
			send(msg,"toApp");
		}
	}
}

MsMac::~MsMac()  {
}
