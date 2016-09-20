/*
 * Channel.cc
 *
 *  Created on: Jul 15, 2014
 *      Author: Thomas Prinz
 * 
 */
 
#include "Channel.h"
#include "MessageTypes.h"

const double Channel::speedOfLightVac = 299792458;

void Channel::addTransInfo(TransInfo* trans){
	int dir = trans->getMessageDirection();
	if(dir==MessageDirection::up || dir==MessageDirection::d2dUp){
		transInfo.first[trans->getRb()].push_front(trans);
	}
	else if(dir==MessageDirection::down || dir==MessageDirection::d2dDown){
		transInfo.second[trans->getRb()].push_front(trans);
	}
}
