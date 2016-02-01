/*
 * VisibilityRegion.h
 *
 *  Created on: May 7, 2014
 *      Author: Thomas Prinz
 */

#pragma once

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
