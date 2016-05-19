/*
 * NeighbourIdMatching.h
 *
 *  Created on: Jul 17, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include "includes.h"
#include <map>

struct Matching {
    int dataStrId;
    int gateId;
};

typedef std::pair<int, int> DataGateIdPair;
typedef std::map<int, DataGateIdPair> NeighbourMap;

class NeighbourIdMatching  {
    private:
        NeighbourMap matching;
	cModule *cell;
	int ownBsId;

        void setupNeighbours(int ownBsId, int maxNumberOfNeighbours, cModule *cell);

    public:
        NeighbourIdMatching(int ownBsId, int maxNumberOfNeighbours, cModule *cell);
        void insert(int bsId, int dataStrId, int gateId);
        int getDataStrId(int bsId);
        int getGateId(int bsId);
        NeighbourMap *getNeighbourMap();
	int getNumberOfMS(int bsId);
        int numberOfNeighbours();
        void showMatching();
};
