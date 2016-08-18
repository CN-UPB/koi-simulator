/**
 * \file EDFRBSched.cc
 *
 * This file contains the implementation of the EDF RB Scheduler.
 *
 * This resource block scheduler schedules packets depending on their deadlines.
 * The packet with the earliest deadline will be scheduled first.
 */

#include "EDFRBSched.h"

Define_Module(EDFRBSched);

bool EDFRBSched::comparator(const KoiData *left, const KoiData *right) const{
	return left->getDeadline()<right->getDeadline();
}
