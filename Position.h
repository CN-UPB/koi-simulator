/*
 * Position.h
 *
 *  Created on: Jul 4, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include <omnetpp.h>

struct Position  {
    double x;
    double y;
};

void doPacking(cCommBuffer *buffer, Position &pos);

void doUnpacking(cCommBuffer *buffer, Position &pos);
