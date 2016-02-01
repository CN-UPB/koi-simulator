/*
 * SimpleChannelCalc.cc
 *
 *  Created on: Jul 9, 2013
 *      Author: Sascha Schmerling
 */

#include "SimpleChannelCalc.h"

/* randomNumber: Interval [0,1] for packet loss */
bool SimpleChannelCalc::calc(int numberOfNops, double randomNumber, double packetLoss)  {
    //simple delay variant
    for(int i = 0; i < numberOfNops; i++)
        asm("nop;");
    if(randomNumber <= packetLoss)
        return false;
    return true;
}
