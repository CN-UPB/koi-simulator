/*
 * Physical.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include "includes.h"

class BsPhy : public omnetpp::cSimpleModule  {
    private:
        int numberOfMobileStations;

    protected:
        void initialize() override;
        void handleMessage(omnetpp::cMessage *msg) override;
};

