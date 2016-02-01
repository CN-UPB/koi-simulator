/*
 * Physical.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include <omnetpp.h>

class BsPhy : public cSimpleModule  {
    private:
        int numberOfMobileStations;

    protected:
        virtual void initialize();
        virtual void handleMessage(cMessage *msg);
};

