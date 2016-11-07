/**
 * @file Coding.cpp
 * @brief Implementation for the Coding class
 */

#include "Coding.h"
#include "util.h"

#include <fstream>
#include <map>
#include <mutex>
#include <string>

using std::map;
using std::string;

int Coding::refBits = 0;
double Coding::reftti = 0;
std::once_flag Coding::tFlag;
Coding::MCSTable Coding::tMCS;

void Coding::loadTable(const string& tpath, int rbBW){
	std::ifstream tableFile;
	tableFile.open(tpath,std::ifstream::in);
	// Read in the reference bitrate
	std::string line;
	std::getline(tableFile,line);
	refBits = std::stoi(line,nullptr);	
	// Read in the reference TTI
	std::getline(tableFile,line);
	reftti = std::stod(line,nullptr);	
	// Discard next line, it only contains the table header
	std::getline(tableFile,line);
	// Read in table lines until eof
	double bwDB = toDB(rbBW);
	int rx = 0;
	int tx = 0;
	int mcs = 0;
	double normSnr = 0.0;
	double refbw = 0.0;
	size_t currPos = 0;
	for(;tableFile.good(); std::getline(tableFile,line)){
		tx = std::stoi(line,&currPos);
		rx = std::stoi(line.substr(currPos),&currPos);
		mcs = std::stoi(line.substr(currPos),&currPos);
		normSnr = std::stod(line.substr(currPos),&currPos) - bwDB;
		refbw = std::stod(line.substr(currPos),&currPos);
		tMCS[tx][rx][mcs] = std::make_tuple(normSnr,refbw);
	}
	tableFile.close();
}

void Coding::init(const string& tpath,double tti, double rbBW){
	tti = tti;
	rbBW = rbBW;
	// Load the table only once per simulation
	std::call_once(tFlag,loadTable,tpath,rbBW);
}

unsigned Coding::getNumBits(double bwMCS){
	return (rbBW/bwMCS)*(tti/reftti)*refBits;
}

unsigned Coding::getRBCapacity(double sinr,int numTx,int numRx){
	// Determine which MCS to use
	map<int,std::tuple<double,double>>& schemes = tMCS[numTx][numRx];
	double mcsBw = 0.0;
	for(auto it = schemes.rbegin(); it!=schemes.rend(); ++it){
		if(sinr>=std::get<0>(it->second)){
			mcsBw = std::get<1>(it->second);
			break;
		}
	}
	return getNumBits(mcsBw);
}
