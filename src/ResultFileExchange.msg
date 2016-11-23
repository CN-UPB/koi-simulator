// This message is used to send pointers to the result file streams 
// from the BS to the local MS.
cplusplus {{
	#include <fstream>
	using OFPtr = std::ofstream*;
}}

class noncobject OFPtr;

message ResultFileExchange  {
	OFPtr ptr;
}
