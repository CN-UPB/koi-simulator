/*
 * Position.cc
 *
 *  Created on: Jul 4, 2013
 *      Author: Sascha Schmerling
 */

#include "Position.h"

void doPacking(cCommBuffer *buffer, Position &pos) {
    buffer->pack(pos.x);
    buffer->pack(pos.y);
}

void doUnpacking(cCommBuffer *buffer, Position &pos) {
    buffer->unpack(pos.x);
    buffer->unpack(pos.y);
}

