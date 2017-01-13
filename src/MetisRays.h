/**
 * @file MetisRays.h
 * @author Michael Meier <mmeier86@mail.upb.de>
 * @brief Classes to hold precomputed values for METIS ray computations
 */

#pragma once

#include "VecNd.h"

#include <array>
#include <complex>
#include <vector>

class Ray{
	private:
		double dirAoA;
		std::complex<double> precompVal;
	
	public:
		Ray(const double dirAoA, std::complex<double> precompVal)
				:dirAoA(dirAoA),precompVal(precompVal){}
		Ray(){}
		static Ray initialize(
				double azimuthASA,
				double azimuthASD,
				double zenithASA,
				double zenithASD,
				double k_0,
				const std::array<double,3>& senderAntennaPos,
				const std::array<double,3>& receiverAntennaPos,
				const std::vector<double>& randomPhase
				);
		virtual std::complex<double> value(double t, double moveAngle,
				double velocity, double k_0);

};

class LOSRay: public Ray{
	public:
		LOSRay(const double dirAoA, std::complex<double> precompVal)
				:Ray(dirAoA,precompVal){}
		LOSRay(){}
		static LOSRay initialize(
				double dirAoA,
				double dirAoD,
				double dirZoA,
				double dirZoD,
				double k_0,
				const std::array<double,3>& senderAntennaPos,
				const std::array<double,3>& receiverAntennaPos,
				double randomPhase
				);
};

class RayCluster{
	private:
		bool los;
		double k;
		std::vector<Ray> rays;
		double sq_P_over_M;
		LOSRay losRay;
		static std::vector<Ray> genNLOSRays(
				size_t numRays,
				double wavenumber,
				const std::vector<double>& zenithASA,
				const std::vector<double>& zenithASD,
				const std::vector<double>& azimuthASA,
				const std::vector<double>& azimuthASD,
				const std::array<double,3>& senderAntennaPos,
				const std::array<double,3>& receiverAntennaPos,
				const VectorNd<double,2>& randomPhase,
				std::vector<int> *subcluster
				);
	public:
		RayCluster(std::vector<Ray>&& rays, double sq_P_over_M,LOSRay losRay,double k)
				: los(true),
				k(k),
				rays(rays),
				sq_P_over_M(sq_P_over_M),
				losRay(losRay){}
		RayCluster(std::vector<Ray>&& rays, double sq_P_over_M)
				: los(false),rays(rays),sq_P_over_M(sq_P_over_M){}
		RayCluster():los(false){}
		static RayCluster initialize(
				size_t numRays,
				double prefactor,
				double wavenumber,
				const std::vector<double>& zenithASA,
				const std::vector<double>& zenithASD,
				const std::vector<double>& azimuthASA,
				const std::vector<double>& azimuthASD,
				const std::array<double,3>& senderAntennaPos,
				const std::array<double,3>& receiverAntennaPos,
				const VectorNd<double,2>& randomPhase,
				std::vector<int> *subcluster
				);
		static RayCluster initialize(
				double k,
				size_t numRays,
				double prefactor,
				double wavenumber,
				const std::vector<double>& zenithASA,
				const std::vector<double>& zenithASD,
				const std::vector<double>& azimuthASA,
				const std::vector<double>& azimuthASD,
				const std::array<double,3>& senderAntennaPos,
				const std::array<double,3>& receiverAntennaPos,
				const VectorNd<double,2>& randomPhase,
				double randomPhaseLOS,
				std::vector<int> *subcluster,
				double dirAoA,
				double dirAoD,
				double dirZoA,
				double dirZoD
				);
		std::complex<double> clusterValue(double t, double moveAngle,
				double velocity, double k_0);
};

