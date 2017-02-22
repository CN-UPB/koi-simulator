/*
 * MsPhy.h
 *
 *  Created on: Jul 1, 2013
 *      Author: Sascha Schmerling
 */

#pragma once

#include "includes.h"

class MsPhy : public omnetpp::cSimpleModule  {
    private:

    protected:
        virtual void initialize();
        virtual void handleMessage(omnetpp::cMessage *msg);
};
