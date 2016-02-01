cplusplus {{
	#include "Position.h"
}}

struct Position;

message BsMsPositions  {
	int bsId;
	Position positions[];
}