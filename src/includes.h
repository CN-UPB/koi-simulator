/**
 * \file includes.h
 *
 * Include file for the HORIZON Omnet++ header.
 *
 * This special include file is necessary because including the omnetpp.h 
 * produces a large number of warnings during compilation, which make finding 
 * warnings from actual model code harder. Thus, the Omnet++/HORIZON 
 * header is included here with suppressed warnings.
 */

#pragma once
#pragma GCC system_header
#include <omnetpp.h>
