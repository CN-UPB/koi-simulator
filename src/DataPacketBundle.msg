cplusplus {{
	#include "KoiData_m.h"
	#include "MessageTypes.h"
}}

message KoiData;
enum MessageDirection;
enum MessageType;

message DataPacketBundle  {
    int msId;
    int bsId;
    KoiData packets[];
    int RBs[];
    int cqi;
    double transPower;
    int messageDirection @enum(MessageDirection);
    kind = MessageType::bundle;
}
