/*
 * util.cc
 *
 *  Created on: Jul 28, 2014
 *      Author: Thomas Prinz
 */

#include "util.h"
#include <fstream>
#include <numeric>
#include <sys/stat.h>
#include <sys/types.h>

using namespace itpp;
using namespace std;

#define BLER_equals_PER 1

int SINR_to_CQI(double minSinr){
	if(minSinr < -6.75){
		return 1;
	}else if(minSinr < -4.75){
		return 2;
	}else if(minSinr < -2.75){
		return 3;
	}else if(minSinr < -1){
		return 4;
	}else if(minSinr < 1){
		return 5;
	}else if(minSinr < 3){
		return 6;
	}else if(minSinr < 5){
		return 7;
	}else if(minSinr < 6.75){
		return 8;
	}else if(minSinr < 8.75){
		return 9;
	}else if(minSinr < 10.75){
		return 10;
	}else if(minSinr < 12.5){
		return 11;
	}else if(minSinr < 14.5){
		return 12;
	}else if(minSinr < 16.25){
		return 13;
	}else if(minSinr < 18){
		return 14;
	}else{
		return 15;
	}
}

// Includes the bit/symbol rate (for MCS and Codingrate) for a given CQI
double getSpectralEfficiency(int CQI){
	switch(CQI){
		case 1:
			return 0.15;
		case 2:
			return 0.23;
		case 3:
			return 0.38;
		case 4:
			return 0.6;
		case 5:
			return 0.88;
		case 6:
			return 1.18;
		case 7:
			return 1.48;
		case 8:
			return 1.91;
		case 9:
			return 2.41;
		case 10:
			return 2.73;
		case 11:
			return 3.32;
		case 12:
			return 3.90;
		case 13:
			return 4.52;
		case 14:
			return 5.12;
		case 15:
			return 5.55;
		default:
			cout << "Error: Invalid CQI" << CQI << endl;
			return -1;
	}
}

double getChannelCapacity(const vector<double>& sinrValues,int subcarriers){
	int numberRB = sinrValues.size();
	double minSinr;
	if(numberRB > 0){
		minSinr = *(std::min_element(sinrValues.begin(), sinrValues.end()));
	}else{
		return 0;
	}
	int CQI = SINR_to_CQI(minSinr);
	int OFDMA_Symbols = 7;	// Fix for LTE
	
	double capacityPerRB = getSpectralEfficiency(CQI) * subcarriers * OFDMA_Symbols;
	double result = capacityPerRB * numberRB;
	// Round down result, we cannot send half bits.
	return floor(result);
}

double getEffectiveSINR(vector<double> sinrValues, vec eesm_beta_values){
	double minSinr = *(std::min_element(sinrValues.begin(), sinrValues.end()));
	int CQI = SINR_to_CQI(minSinr);
	double beta;
	switch(CQI){
		case 1:
			beta = eesm_beta_values(0);
			break;
		case 2:
			beta = eesm_beta_values(1);
			break;
		case 3:
			beta = eesm_beta_values(2);
			break;
		case 4:
			beta = eesm_beta_values(3);
			break;
		case 5:
			beta = eesm_beta_values(4);
			break;
		case 6:
			beta = eesm_beta_values(5);
			break;
		case 7:
			beta = eesm_beta_values(6);
			break;
		case 8:
			beta = eesm_beta_values(7);
			break;
		case 9:
			beta = eesm_beta_values(8);
			break;
		case 10:
			beta = eesm_beta_values(9);
			break;
		case 11:
			beta = eesm_beta_values(10);
			break;
		case 12:
			beta = eesm_beta_values(11);
			break;
		case 13:
			beta = eesm_beta_values(12);
			break;
		case 14:
			beta = eesm_beta_values(13);
			break;
		case 15:
			beta = eesm_beta_values(14);
			break;
		default:
			cout << "ERROR. CQI value out of range." << endl;
	}
	double temp = 0;
	for(uint i = 0; i < sinrValues.size(); i++){
		temp += exp(-1.0 * (sinrValues.at(i) / beta));
	}
	temp = (temp / sinrValues.size());
	double effective_sinr = -1.0 * beta * log(temp);
	return effective_sinr;
}

double getBler(int cqi, double sinr, omnetpp::cSimpleModule* module){
	
	static string blerIn = module->par("bler_table");
	static mat blerTable = mat(blerIn);
	
	/**
	cout << blerTable.rows() << " " << blerTable.cols() << endl;
	cout << "CQI: " << cqi << endl;
	cout << "SINR: " << sinr << endl;
	**/
	
	double bler=0;
	int line=0;

	// Select the mutual information value according to the snr2MI lookup table
	// w.r.t. the instRbSnr: For all avaialble modulation types...
	for (line = 0; line < blerTable.rows(); line++) {
		// Check if the instRbSnr is not greater than the MI threshold.
		//cout<<"sinr: "<<sinr<<" table Sinr: "<<blerTable(line,2*cqi)<<endl;
		if (!(sinr > blerTable(line,2*(cqi-1)))) {
			// Check if the instSNR isn't to low too reach the first step in order to
			// avoid a NULL pointer.
			if (line > 0) {
				bler = blerTable(line,2*cqi-1);
				//cout << "BLER: " << bler << endl;
			}
			else {
				bler = 1;
			}

			break;
		}
	}

	if(line==0)
		bler = 1;
	else if(line == blerTable.rows())
		bler = 0;
		
	return bler;
}

double getPer(vec bler){
	if(BLER_equals_PER){
		return bler(0);
	}else{
		double result = 1;
		for(int i = 0; i < bler.size(); i++){
			result = result * bler(i);
		}
		return result;
	}
}

double toDB(double val){
	return 10*std::log10(val);
}

int lcm(int a, int b){
  int temp = itpp::gcd(a,b);
  return temp ? (a/temp)*b : 0;
}

int integerOoM(double val){
  // Find out which order of magnitude moves the decimal point out of the 
  // number, leaving only a natural number behind.
  // If the value has no fractional component, return 0 (10^0 = 1)
  double tmp;
  if(std::modf(val,&tmp)==0.0){
    return 0;
  }
  // Kickstart the process with the tenth logarithm of the value
  int oom = std::floor(std::log10(val));
  double curr(val*std::pow(10,oom));
  while(std::modf(curr,&tmp)!=0.0){
    curr *= 10;
    ++oom;
  }
  return oom;
}

double lcmSequence(const vector<double>& elems){
  // First, we need all the input numbers as integers, because LCM
  // (and GCD, which the algorithm is based on) are only defined for natural
  // numbers. Rounding is out of the question. 
  // Instead, represent all numbers in the format x*10^y with such y that 
  // x is an integer.
  int largestOoM(0);
  int curr;
  for(const double& val:elems){
    curr = integerOoM(val);
    if(curr>largestOoM){
      largestOoM = curr;
    }
  }
  auto makeInt = [&largestOoM](double v) 
    -> int {return v*std::pow(10,largestOoM);};
  vector<int> intVals(elems.size());
  std::transform(elems.begin(),elems.end(),intVals.begin(),makeInt);
  return lcmSequence(intVals)*std::pow(10,-largestOoM);
}

int lcmSequence(const vector<int>& elems){
  return std::accumulate(++(elems.begin()),elems.end(),elems[0],lcm);
}

std::ofstream getResultFile(std::string& fname){
	int run = std::stoi(omnetpp::getEnvir()->getConfig()->substituteVariables("${runnumber}"));
	fname = "./results/"+fname+"_run-"+std::to_string(run)+".dat";
	struct stat info;
	if(stat("./results",&info)!=0){
		// The results file does not exist, create it
		mkdir("./results",0770);
	}
	std::ofstream resfile;
	resfile.open(fname,std::ofstream::trunc);
	return resfile;
}
