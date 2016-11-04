/**
 * @file Coding.h 
 * @brief A class for computing a number of bits transferable over a given RB
 */

#include "includes.h"

#include <map>
#include <mutex>
#include <string>
#include <tuple>

class Coding{
	public:
		using MCSTable = std::map<int,std::map<int,
					std::map<int,std::tuple<double,double>>>>;

	private:
		/**
		 * Bandwith per Resource block in Hz
		 */
		double rbBW;
		/**
		 * Rate bits/tti used in MCS table calculation
		 */
		static int refBits;
		/**
		 * TTI in seconds used in MCS table calculation
		 */
		static simtime_t reftti;
		/**
		 * Flag signaling that table loading should only be conducted once
		 */
		static std::once_flag tFlag;
		/**
		 * Table with normalized SNR values and RB bandwidths per antenna and MCS
		 */
		static Coding::MCSTable tMCS;
		/**
		 * The TTI length for the current system in seconds
		 */
		simtime_t tti;

	private:
		unsigned getNumBits(double bwMCS);
		static void loadTable(const std::string& tpath, int rbBW);

	public:
		Coding(){}
		void init(const std::string& tpath,simtime_t tti, double rbBW);
		unsigned getRBCapacity(double sinr,int numTx,int numRx);

};