/*
 * util.h
 *
 *  Created on: Jul 28, 2014
 *      Author: Thomas Prinz
 */

#pragma once

#include "includes.h"
#include <itpp/itbase.h>
#include <fstream>
#include <algorithm>

using namespace itpp;
using namespace std;

int SINR_to_CQI(double minSinr);

// Includes the bit/symbol rate (for MCS and Codingrate) for a given CQI
double getSpectralEfficiency(int CQI);

double getChannelCapacity(const vector<double>& sinrValues, int subcarriers);

double getBler(int cqi, double sinr, cSimpleModule* module);

double getEffectiveSINR(vector<double> sinrValues, vec eesm_beta_values);

double getPer(vec bler);

int lcm(int a, int b);

double lcmSequence(const vector<double>& elems);

int lcmSequence(const vector<int>& elems);

int integerOoM(double val);

std::ofstream getResultFile(std::string& fname);
