cplusplus {{
	#include "cluster.h"
}}

// Extends void is needed, so that omnet knows, that it is not
// derived from some interesting class ( for code generation )
class Cluster extends void;

message ClusterMessage  {
    Cluster c;
}
