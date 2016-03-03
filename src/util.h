/*
 * util.h
 *
 *  Created on: Jul 28, 2014
 *      Author: Thomas Prinz
 */

#pragma once

#include "includes.h"
#include <itpp/itbase.h>
#include <algorithm>

using namespace itpp;
using namespace std;

int SINR_to_CQI(double minSinr);

// Includes the bit/symbol rate (for MCS and Codingrate) for a given CQI
double getSpectralEfficiency(int CQI);

double getChannelCapacity(vector<double> sinrValues);

double getBler(int cqi, double sinr, cSimpleModule* module);

double getEffectiveSINR(vector<double> sinrValues, vec eesm_beta_values);

double getPer(vec bler);
