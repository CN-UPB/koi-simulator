/*
 * VisibilityRegion.h
 *
 *  Created on: May 7, 2014
 *      Author: Thomas Prinz
 */

#pragma once

#include "includes.h"
#include <vector>

struct VisibilityRegion  {
    int origin;
    double x;
    double y;
    std::vector<int> bsIds;
    int clusterId;
    bool last;
    int assigner;
};

void doPacking(cCommBuffer *buffer, VisibilityRegion &vr);

void doUnpacking(cCommBuffer *buffer, VisibilityRegion &vr);
