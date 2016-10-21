cplusplus {{
	#include "Channel.h"
	using ChPtr = Channel*;
}}

class noncobject ChPtr;

message PointerExchange  {
	ChPtr ptr;
}
