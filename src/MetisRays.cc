/**
 * @file MetisRays.cc
 * @brief Implmentation of the MetisRay header
 */

#pragma once

#include "includes.h"
#include "MetisRays.h"
#include "VecNd.h"

#include <array>
#include <complex>
#include <vector>

using std::array;
using std::complex;
using std::vector;

Ray Ray::initialize(
		double azimuthASA,
		double azimuthASD,
		double zenithASA,
		double zenithASD,
		double k_0,
		const array<double,3>& senderAntennaPos,
		const array<double,3>& receiverAntennaPos,
		const vector<double>& randomPhase
		){
	double AoA[3];
	double AoD[3];
	double receiverGain;
	double senderGain;
	AoA[0] = sin(zenithASA*pi/180) * cos(azimuthASA*pi/180);
	AoA[1] = sin(zenithASA*pi/180) * sin(azimuthASA*pi/180);
	AoA[2] = cos(zenithASA*pi/180);

	AoD[0] = sin(zenithASD*pi/180) * cos(azimuthASD*pi/180);
	AoD[1] = sin(zenithASD*pi/180) * sin(azimuthASD*pi/180);
	AoD[2] = cos(zenithASD*pi/180);

	// Ray arrival component
	complex<double> expArrival = exp( complex<double>(0.0,k_0 * (AoA[0] * receiverAntennaPos[0] + AoA[1] * receiverAntennaPos[1] + AoA[2] * receiverAntennaPos[2])) );
	// Ray Departure component
	complex<double> expDeparture = exp( complex<double>(0.0,k_0 * (AoD[0] * senderAntennaPos[0] + AoD[1] * senderAntennaPos[1] + AoD[2] * senderAntennaPos[2])) );

	// Ray Polarization
	receiverGain = getMSGain(azimuthASA*pi/180, zenithASA*pi/180);
	senderGain = getBSGain(azimuthASD*pi/180, zenithASD*pi/180);
	complex<double> pol = receiverGain * senderGain * exp(complex<double>(0, randomPhase[0]));
	return Ray(azimuthASA*pi/180,expArrival*expDeparture*pol);
}

std::complex<double> Ray::value(double t, double moveAngle,
		double velocity, double k_0){
	std::complex<double> doppler = exp( complex<double>(0,k_0 * velocity * cos(dirAoA - moveAngle) * t ) );
	return doppler*precompVal;
}

LOSRay LOSRay::initialize(
		double dirAoA,
		double dirAoD,
		double dirZoA,
		double dirZoD,
		double k_0,
		const array<double,3>& senderAntennaPos,
		const array<double,3>& receiverAntennaPos,
		double randomPhase
		){
	double AoA[3];
	double AoD[3];
	double receiverGain;
	double senderGain;
	AoA[0] = sin(dirZoA) * cos(dirAoA);
	AoA[1] = sin(dirZoA) * sin(dirAoA);
	AoA[2] = cos(dirZoA);

	AoD[0] = sin(dirZoD) * cos(dirAoD);
	AoD[1] = sin(dirZoD) * sin(dirAoD);
	AoD[2] = cos(dirZoD);

	complex<double> expArrival = exp( complex<double>(0.0,k_0 * (AoA[0] * receiverAntennaPos[0] + AoA[1] * receiverAntennaPos[1] + AoA[2] * receiverAntennaPos[2])) );
	complex<double> expDeparture = exp( complex<double>(0.0,k_0 * (AoD[0] * senderAntennaPos[0] + AoD[1] * senderAntennaPos[1] + AoD[2] * senderAntennaPos[2])) );

	receiverGain = getMSGain(dirAoA, dirZoA);
	senderGain = getBSGain(dirAoD, dirZoD);
	complex<double> pol = receiverGain * senderGain * exp(complex<double>(0, randomPhase));
	
	return LOSRay(dirAoA,expArrival*expDeparture*pol);
}

vector<Ray> RayCluster::genNLOSRays(
		size_t numRays,
		double wavenumber,
		const vector<double>& zenithASA,
		const vector<double>& zenithASD,
		const vector<double>& azimuthASA,
		const vector<double>& azimuthASD,
		const array<double,3>& senderAntennaPos,
		const array<double,3>& receiverAntennaPos,
		const VectorNd<double,2>& randomPhase,
		vector<int> *subcluster
		){
	vector<Ray> res(numRays);
	size_t rayIdx;
	for(size_t m = 0; m < numRays; m++){
		rayIdx = (subcluster!=nullptr) ? (*subcluster)[m] : m;
		res[rayIdx] = Ray::initialize(
				azimuthASA[rayIdx],
				azimuthASD[rayIdx],
				zenithASA[rayIdx],
				zenithASD[rayIdx],
				wavenumber,
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
		const array<double,3>& senderAntennaPos,
		const array<double,3>& receiverAntennaPos,
		const VectorNd<double,2>& randomPhase,
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
		const array<double,3>& senderAntennaPos,
		const array<double,3>& receiverAntennaPos,
		const VectorNd<double,2>& randomPhase,
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
				senderAntennaPos,
				receiverAntennaPos,
				randomPhase,
				subcluster
				));
	// Precompute LOS ray
	LOSRay losRay(LOSRay::initialize(
				dirAoA,dirAoD,dirZoA,dirZoD,
				wavenumber,
				senderAntennaPos,receiverAntennaPos,
				randomPhaseLOS));
	return RayCluster(std::move(rays),prefactor,std::move(losRay),k);
}

complex<double> RayCluster::clusterValue(double t, double moveAngle,
		double velocity, double k_0){
	complex<double> val(0.0,0.0);
	// Add up all NLOS rays in the cluster
	for(Ray& r:rays){
		val += r.value(t,moveAngle,velocity,k_0);
	}
	// Multiply with prefactor
	val *= sq_P_over_M;
	if(los){
		// Cmmunicating stations have Line of Sight, add a LOS ray
		val *= std::sqrt(1.0/(k+1.0));
		val += std::sqrt(k/(k+1.0))*losRay.value(t,moveAngle,velocity,k_0);
	}
	return val;
}

