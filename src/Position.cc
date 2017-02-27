/*
 * Position.cc
 *
 *  Created on: Jul 4, 2013
 *      Author: Sascha Schmerling
 */

#include "Position.h"

using omnetpp::cCommBuffer;

double Position::dist2d(const Position& other) const{
	double distX = pow(x-other.x,2);
	double distY = pow(y-other.y,2);
	return sqrt(distX+distY);
}

double Position::dist3d(const Position& other) const{
	double distX = pow(x-other.x,2);
	double distY = pow(y-other.y,2);
	double distZ = pow(z-other.z,2);
	return sqrt(distX+distY+distZ);
}

double Position::distance(const Position& other) const{
	return y==other.y ? dist2d(other) : dist3d(other);
}

void doPacking(cCommBuffer *buffer, Position &pos) {
    buffer->pack(pos.x);
    buffer->pack(pos.y);
    buffer->pack(pos.z);
}

void doUnpacking(cCommBuffer *buffer, Position &pos) {
    buffer->unpack(pos.x);
    buffer->unpack(pos.y);
    buffer->unpack(pos.z);
}

