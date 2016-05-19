cplusplus {{
	#include "KoiData_m.h"
}}

message KoiData;

message DataPacketBundle  {
    int msId;
    int bsId;
    KoiData packets[];
    int RBs[];
    int cqi;
    double transPower;
}
