/*
 * MsMac.cc
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#include "MsMac.h"
#include "DataPacket_m.h"
#include "PositionExchange_m.h"
#include "DataPacketBundle_m.h"
#include "TransmitRequest_m.h"
#include "Schedule_m.h"
#include "SINR_m.h"
#include <stdlib.h>
#include <cmath>

Define_Module(MsMac);

inline simtime_t MsMac::positionResendTime()  {
    if(simTime() == initOffset)  {
        return simTime() + (positionResendInterval * tti) - (2 * epsilon); //send with the schedules; one sync point
    }
    else  {
        return simTime() + positionResendInterval * tti;
    }
}

void MsMac::updateDisplayString()  {
    std::stringstream converter;
    //correction the center of the cell to (radius,radius)
    converter << "p=" << msPosition.x - initBsPos.x + radius << "," << msPosition.y - initBsPos.y + radius << ";bgb=200,300;i=device/cellphone;";
    std::string displayString = converter.str();

    // Setting a module's position, icon and status icon:
    cDisplayString& dispStr = getParentModule()->getDisplayString();
    dispStr.parse(displayString.c_str());
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

void MsMac::initialize()  {
    	positionResendInterval = par("positionResendInterval");
    	msId = par("msId");
    	bsId = par("bsId");
    	epsilon = par("epsilon");
   	radius = par("radius");
	initBsPos.x = par("initBsXPos");
	initBsPos.y = par("initBsYPos");
    	currentChannel = par("currentChannel");
	double initPosAlpha = par("initPosAlpha");
	double initPosBeta = par("initPosBeta");
	double initPosGamma = par("initPosGamma");
	int initQuadrant = par("initQuadrant");
   	initOffset = par("initOffset");
   	tti = par("tti");
    	downResourceBlocks = par("downResourceBlocks");
    	packetLength = par("packetLength");

        //msPosition = initMsPosition(initQuadrant, initPosAlpha, initPosBeta, initPosGamma); //for random MS positions in hexagonal cell
    	msPosition = initMsPositionLinear(); //for MS position along a straight road
	//for Tkenv
	updateDisplayString();
	
	//TODO: get from ini
	SINR_ = zeros(100);

    //resend the ms position every x times to the BsMac layer
    scheduleAt(simTime() + initOffset + tti - epsilon, new cMessage("RESEND_POS")); //originally set to simTime() + initOffset 

    //every tti send the bs the transmit requests
    scheduleAt(simTime() + initOffset + tti, new cMessage("GEN_TRANSMIT_REQUEST"));
    
    // Send MS Position once at the very beginning for cluster generation
    PositionExchange *posEx = new PositionExchange("MS_POS_UPDATE");
    posEx->setId(msId);
    posEx->setPosition(msPosition);
    send(posEx, "toBsMac");
}

void MsMac::handleMessage(cMessage *msg)  {
    if(msg->isName("SCHEDULE"))  {
        Schedule *schedule = (Schedule *) msg;
        
        //drop schedule which not corresponds to the current channel of this MS
	//TODO execute only if this is a DOWN schedule
        if(schedule->getChannel() == currentChannel)  {
            DataPacketBundle *packetBundle = new DataPacketBundle("DATA_BUNDLE"); //only send out one bundle of packets
            packetBundle->setMsId(msId);
            packetBundle->setBsId(bsId);
            int numberOfRB = 0;
            
            vector<double> sinr_values;
            vector<double> RBs;
            
            // Get the number of Ressource Blocks, that are assigned to this MS
            for(int i = 0; i < (int) schedule->getUpScheduleArraySize(); ++i)  {
                int currentMsId = schedule->getUpSchedule(i);
		//std::cout << "Up Schedule MS Id: " << currentMsId << std::endl;
                if(currentMsId == msId){
					sinr_values.push_back(SINR_(i));
					RBs.push_back(i);
					numberOfRB++;
				}
            }
            
            double channel_capacity = getChannelCapacity(sinr_values);
            int cqi;
            if(sinr_values.size() > 0){
				cqi = SINR_to_CQI(*(std::min_element(sinr_values.begin(), sinr_values.end())));
			}else{
				cqi = 1;
			}
            
            if(channel_capacity > 0)  {
		int packetsToSend;
		if(packetQueue.length()>channel_capacity/packetLength){
			packetsToSend = channel_capacity/packetLength;
		}
		else{
			packetsToSend = packetQueue.length();
		}
                packetBundle->setPacketsArraySize(packetsToSend);
                //send packets to the bs if the own ms id is in the up schedule
                ASSERT(schedule->getBsId() == bsId);
                for(uint i = 0; i < packetBundle->getPacketsArraySize(); i++){
			KoiData *packet = dynamic_cast<KoiData*>(packetQueue.pop());
			packetBundle->setPackets(i, *packet);
			delete packet;
		}
		packetBundle->setRBsArraySize(RBs.size());
		packetBundle->setCqi(cqi);
		for(uint i = 0; i < packetBundle->getRBsArraySize(); i++){
			packetBundle->setRBs(i,RBs[i]);
		}
                //ev << "MS " << msId << " sending " << packetBundle->getPacketsArraySize() << "packets to the BS " << schedule->getBsId() << " channel" << endl;
                sendDelayed(packetBundle, epsilon, "toPhy");
            }
        
		/*
        //drop schedule which not corresponds to the current channel of this MS
        if(schedule->getChannel() == currentChannel)  {
            DataPacketBundle *packetBundle = new DataPacketBundle("DATA_BUNDLE"); //only send out one bundle of packets
            packetBundle->setMsId(msId);
            packetBundle->setBsId(bsId);
            int neededSize = 0;
            for(int i = 0; i < (int) schedule->getUpScheduleArraySize(); ++i)  { //TODO could be done more efficient
                int currentMsId = schedule->getUpSchedule(i);
                if(currentMsId == msId)
                    neededSize++;
            }
            if(neededSize > 0)  {
                packetBundle->setPacketsArraySize(neededSize);
                //send packets to the bs if the own ms id is in the up schedule
                ASSERT(schedule->getBsId() == bsId);
                int curBundlePos = 0;
                for(int i = 0; i < (int) schedule->getUpScheduleArraySize(); ++i)  {
                    int currentMsId = schedule->getUpSchedule(i);
                    //ev << "Slot " << i << " - " << currentMsId << endl;
                    if(currentMsId == msId)  {
                        DataPacket *packet = (DataPacket *) packetQueue.pop();
                        packet->setResourceBlock(i);
                        packetBundle->setPackets(curBundlePos, *packet);
                        curBundlePos++;
                    }
                }
                //ev << "MS " << msId << " sending " << packetBundle->getPacketsArraySize() << "packets to the BS " << schedule->getBsId() << " channel" << endl;
                sendDelayed(packetBundle, epsilon, "toPhy");
            }
            * */
        }
        delete schedule;
    }
    else if(msg->isName("GEN_TRANSMIT_REQUEST"))  {
        //send the bs the sending wishes
        int length = packetQueue.length();
        if(length != 0)  { //only send the transmit req when the queue is not empty
            TransmitRequest *transReq = new TransmitRequest("TRANSMIT_REQUEST");
            transReq->setId(msId);
            transReq->setPacketCount(packetQueue.length());
            send(transReq, "toBsMac");
        }
        scheduleAt(simTime() + tti, msg);
    }
    else if(msg->isName("RESEND_POS"))  {
	
        PositionExchange *posEx = new PositionExchange("MS_POS_UPDATE");
		msPosition.x = msPosition.x + (positionResendInterval/1000.0) * velocity.at(0);
		msPosition.y = msPosition.y + (positionResendInterval/1000.0) * velocity.at(1);
		//cout << "MS Pos: " << msPosition.x << " " << msPosition.y << endl;
		//cout << "MS Vel: " << velocity.at(0) << " " << velocity.at(1) << endl;
        posEx->setId(msId);
        posEx->setPosition(msPosition);
        send(posEx->dup(), "toPhy"); //send position to the own channel module
        send(posEx, "toBsMac");
        scheduleAt(positionResendTime(), msg);
    }
    else if(msg->isName("SINR_ESTIMATION"))  {
		SINR *sinrMessage = (SINR *) msg;
		// Add MS ID Information, MS Channel doesnt have that information.
		sinrMessage->setMsId(msId);
		vec sinr_new(downResourceBlocks);
		for(int i = 0;i < downResourceBlocks; i++){
			sinr_new.set(i,sinrMessage->getSINR(i));
		}
		SINR_ = sinr_new;
		send(sinrMessage, "toBsMac");
    }
    //DataPacket
    else if(msg->arrivedOn("fromApp"))  {
	packetQueue.insert((DataPacket *) msg);
    }
    //DataPacket
    else if(msg->arrivedOn("fromPhy"))  {
	// Unpack the data bundle and forward data packets to app
	if(msg->isName("DATA_BUNDLE")){
		DataPacketBundle *bundle = dynamic_cast<DataPacketBundle*>(msg);
		for(unsigned int i=0; i<bundle->getPacketsArraySize(); i++){
			send(bundle->getPackets(i).dup(),"toApp");
		} 
	}
	delete msg;
    }
}

MsMac::~MsMac()  {
    //ev << "BS " << bsId << "; MS " << msId << ": Queue " << packetQueue.length() << endl;
}
