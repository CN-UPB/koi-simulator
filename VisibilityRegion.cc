/*
 * VisibilityRegion.cc
 *
 *  Created on: May 7, 2014
 *      Author: Thomas Prinz
 */

#include "VisibilityRegion.h"

void doPacking(cCommBuffer *buffer, VisibilityRegion &vr) {
    buffer->pack(vr.origin);
    buffer->pack(vr.x);
    buffer->pack(vr.y);
    buffer->pack(vr.bsIds.size());
    for(uint i = 0;i < vr.bsIds.size(); i++){
       buffer->pack(vr.bsIds[i]);
    }
    buffer->pack(vr.clusterId);
    buffer->pack(vr.last);
    buffer->pack(vr.assigner);
}

void doUnpacking(cCommBuffer *buffer, VisibilityRegion &vr) {
    buffer->unpack(vr.origin);
    buffer->unpack(vr.x);
    buffer->unpack(vr.y);
    std::vector<int>::size_type length;
    buffer->unpack(length);
    vr.bsIds.reserve(length);
    for(uint i = 0;i < length; i++){
        int temp;
        buffer->unpack(temp);
        vr.bsIds.push_back(temp);
    }
    buffer->unpack(vr.clusterId);
    buffer->unpack(vr.last);
    buffer->unpack(vr.assigner);
}

