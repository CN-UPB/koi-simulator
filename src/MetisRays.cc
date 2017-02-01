/**
 * @file MetisRays.cc
 * @brief Implmentation of the MetisRay header
 */

#include "includes.h"
#include "MetisRays.h"
#include "VecNd.h"

#include <array>
#include <complex>
#include <vector>

using std::array;
using std::complex;
using std::vector;

double getMSGain(double AoA, double EoA){ // Hertzian dipole; TODO: Change the return value to a vector for getting both theta and phi components
	//return (-1.0 * sin(EoA)); 
	return -1.0;				// For Z-oriented Hertzian dipole, F_phi is always zero and F_theta is -sin(theta), theta is in radians
}

double getBSGain(double AoD, double EoD){ // Hertzian dipole; TODO: Change the return value to a vector for getting both theta and phi components, EoD and AoD are in radians
	//return (-1.0 * sin(EoD));
	return -1.0;
}

Ray Ray::initialize(
		double azimuthASA,
		double azimuthASD,
		double zenithASA,
		double zenithASD,
		double k_0,
		double moveAngle,
		double velocity,
		const array<double,3>& senderAntennaPos,
		const array<double,3>& receiverAntennaPos,
		const double randomPhase
		){
	double AoA[3];
	double AoD[3];
	double receiverGain;
	double senderGain;
	AoA[0] = sin(zenithASA) * cos(azimuthASA);
	AoA[1] = sin(zenithASA) * sin(azimuthASA);
	AoA[2] = cos(zenithASA);

	AoD[0] = sin(zenithASD) * cos(azimuthASD);
	AoD[1] = sin(zenithASD) * sin(azimuthASD);
	AoD[2] = cos(zenithASD);

	// Ray arrival component
	complex<double> expArrival = exp( complex<double>(0.0,k_0 * (AoA[0] * receiverAntennaPos[0] + AoA[1] * receiverAntennaPos[1] + AoA[2] * receiverAntennaPos[2])) );
	// Ray Departure component
	complex<double> expDeparture = exp( complex<double>(0.0,k_0 * (AoD[0] * senderAntennaPos[0] + AoD[1] * senderAntennaPos[1] + AoD[2] * senderAntennaPos[2])) );

	// Ray Polarization
	receiverGain = getMSGain(azimuthASA*PI/180, zenithASA*PI/180);
	senderGain = getBSGain(azimuthASD*PI/180, zenithASD*PI/180);
	complex<double> pol = receiverGain * senderGain * exp(complex<double>(0, randomPhase));

	// Precomputable part of the doppler component
	array<double,3> velVec = {cos(moveAngle),sin(moveAngle),0};
	double velComp = k_0*(AoA[0]*velVec[0]+AoA[1]*velVec[1]+AoA[2]*velVec[2])*velocity;
	return Ray(velComp,expArrival*expDeparture*pol);
}

inline std::complex<double> Ray::value(const double t) const {
	return precompVal*exp( complex<double>(0,velComp * t));
}

LOSRay LOSRay::initialize(
		double dirAoA,
		double dirAoD,
		double dirZoA,
		double dirZoD,
		double k_0,
		double moveAngle,
		double velocity,
		const array<double,3>& senderAntennaPos,
		const array<double,3>& receiverAntennaPos,
		double randomPhase
		){
	double AoA[3];
	double AoD[3];
	double receiverGain;
	double senderGain;
	AoA[0] = sin(dirZoA*(180/PI)) * cos(dirAoA*(180/PI));
	AoA[1] = sin(dirZoA*(180/PI)) * sin(dirAoA*(180/PI));
	AoA[2] = cos(dirZoA*(180/PI));

	AoD[0] = sin(dirZoD*(180/PI)) * cos(dirAoD*(180/PI));
	AoD[1] = sin(dirZoD*(180/PI)) * sin(dirAoD*(180/PI));
	AoD[2] = cos(dirZoD*(180/PI));

	complex<double> expArrival = exp( complex<double>(0.0,k_0 * (AoA[0] * receiverAntennaPos[0] + AoA[1] * receiverAntennaPos[1] + AoA[2] * receiverAntennaPos[2])) );
	complex<double> expDeparture = exp( complex<double>(0.0,k_0 * (AoD[0] * senderAntennaPos[0] + AoD[1] * senderAntennaPos[1] + AoD[2] * senderAntennaPos[2])) );

	receiverGain = getMSGain(dirAoA, dirZoA);
	senderGain = getBSGain(dirAoD, dirZoD);
	complex<double> pol = receiverGain * senderGain * exp(complex<double>(0, randomPhase));

	// Precomputable part of the doppler component
	array<double,3> velVec = {cos(moveAngle),sin(moveAngle),0};
	double velComp = k_0*(AoA[0]*velVec[0]+AoA[1]*velVec[1]+AoA[2]*velVec[2])*velocity;
	
	return LOSRay(velComp,expArrival*expDeparture*pol);
}

vector<Ray> RayCluster::genNLOSRays(
		size_t numRays,
		double wavenumber,
		const vector<double>& zenithASA,
		const vector<double>& zenithASD,
		const vector<double>& azimuthASA,
		const vector<double>& azimuthASD,
		double moveAngle,
		double velocity,
		const array<double,3>& senderAntennaPos,
		const array<double,3>& receiverAntennaPos,
		const vector<double>& randomPhase,
		vector<int> *subcluster
		){
	vector<Ray> res(numRays);
	size_t rayIdx;
	for(size_t m = 0; m < numRays; m++){
		rayIdx = (subcluster!=nullptr) ? (*subcluster)[m] : m;
		res[m] = Ray::initialize(
				azimuthASA[rayIdx],
				azimuthASD[rayIdx],
				zenithASA[rayIdx],
				zenithASD[rayIdx],
				wavenumber,
				moveAngle,
				velocity,
				senderAntennaPos,
				receiverAntennaPos,
				randomPhase[rayIdx]
				);
	}
	return res;
}

RayCluster RayCluster::initialize(
		size_t numRays,
		double prefactor,
		double wavenumber,
		const vector<double>& zenithASA,
		const vector<double>& zenithASD,
		const vector<double>& azimuthASA,
		const vector<double>& azimuthASD,
		double moveAngle,
		double velocity,
		const array<double,3>& senderAntennaPos,
		const array<double,3>& receiverAntennaPos,
		const vector<double>& randomPhase,
		vector<int> *subcluster
		){
	// Precompute NLOS ray components
	vector<Ray> rays(genNLOSRays(
				numRays,
				wavenumber,
				zenithASA,
				zenithASD,
				azimuthASA,
				azimuthASD,
				moveAngle,
				velocity,
				senderAntennaPos,
				receiverAntennaPos,
				randomPhase,
				subcluster
				));
	return RayCluster(std::move(rays),prefactor);
}

RayCluster RayCluster::initialize(
		double k,
		size_t numRays,
		double prefactor,
		double wavenumber,
		const vector<double>& zenithASA,
		const vector<double>& zenithASD,
		const vector<double>& azimuthASA,
		const vector<double>& azimuthASD,
		double moveAngle,
		double velocity,
		const array<double,3>& senderAntennaPos,
		const array<double,3>& receiverAntennaPos,
		const vector<double>& randomPhase,
		double randomPhaseLOS,
		vector<int> *subcluster,
		double dirAoA,
		double dirAoD,
		double dirZoA,
		double dirZoD
		){
	// Precompute NLOS ray components
	vector<Ray> rays(genNLOSRays(
				numRays,
				wavenumber,
				zenithASA,
				zenithASD,
				azimuthASA,
				azimuthASD,
				moveAngle,
				velocity,
				senderAntennaPos,
				receiverAntennaPos,
				randomPhase,
				subcluster
				));
	// Precompute LOS ray
	LOSRay losRay(LOSRay::initialize(
				dirAoA,dirAoD,dirZoA,dirZoD,
				wavenumber,
				moveAngle,
				velocity,
				senderAntennaPos,receiverAntennaPos,
				randomPhaseLOS));
	return RayCluster(std::move(rays),prefactor,std::move(losRay),k);
}

complex<double> RayCluster::clusterValue(const double t) const {
	complex<double> val(0.0,0.0);
	// Add up all NLOS rays in the cluster
	for(const Ray& r:rays){
		val += r.value(t);
	}
	// Multiply with prefactor
	val *= sq_P_over_M;
	if(los){
		// Cmmunicating stations have Line of Sight, add a LOS ray
		val *= std::sqrt(1.0/(k+1.0));
		val += std::sqrt(k/(k+1.0))*losRay.value(t);
	}
	return val;
}

