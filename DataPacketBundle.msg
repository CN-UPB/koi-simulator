cplusplus {{
	#include "DataPacket_m.h"
}}

message DataPacket;

message DataPacketBundle  {
    int msId;
    int bsId;
	DataPacket packets[];
	int RBs[];
	int cqi;
}
