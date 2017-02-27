/*
 * Position.h
 *
 *  Created on: Jul 4, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include "includes.h"

struct Position  {
    double x;
    double y;
		double z;
		double distance(const Position& other) const;

	private:
		double dist2d(const Position& other) const;
		double dist3d(const Position& other) const;
};

void doPacking(omnetpp::cCommBuffer *buffer, Position &pos);

void doUnpacking(omnetpp::cCommBuffer *buffer, Position &pos);
