// Representation of a wireless communication cell
//
// This module represents a communication cell, holding a single BaseStation 
// as well as a number of MobileStation instances. In addition, each cell has 
// a Scheduler module. Each cell also has a connection to all other cells, but 
// at the moment, those are only used for information traffic. There is no 
// simulated traffic between cells.
//
// <h3>Submodules</h3>
// The modules combined in a LteCell fullfill the following tasks:
//
// <p><strong>bs:</strong> Each LteCell has execatly one <em>BaseStation</em> module.
// This module serves as a transceiver for simulated traffic packets between 
// MobileStation module in the cell. Additionally, the BaseStation is connected 
// to all other LteCell modules in the simulation. At the moment, this
// connection is only used for information traffic, for example the TransInfo 
// messages from neigbouring cells.
//
// <p><strong>ms:</strong> The <em>MobileStation</em> modules in a cell are traffic 
// generators of the simulation. They generate and consume simulated traffic 
// in the form of KoiData packets. Each MobileStation may communicate with any 
// other MobileStation in the same cell. They communicate either via the local 
// BaseStation or directly via D2D. Their count in each cell is determined 
// at the beginning of each simulation and does not change during a simulation 
// run.
//
// <p><strong>sched:</strong> The <em>Scheduler</em> module is responsible for 
// scheduling simulated traffic packets inside a single cell. They are informed 
// about all communication streams in their respective cell at the start 
// of a simulation run. They ensure that any resource block is only used 
// by a single (local!) transmission in any given TTI. 
// <h3>Connections</h3>
// Each LteCell has a bidirectional connection to all other cells in the 
// simulation. Internally, that connection is served by the BaseStation module 
// as the local organizer. In addition, all local MobileStation modules 
// are connected to each other via MsChannel channels. They are also connected 
// to the local BaseStation via BsChannel modules for the uplink to the 
// BaseStation and MsChannel modules for the downlink from BaseStation to 
// MobileStation. Additionaly, all modules, the BaseStation as well as the 
// MobileStations, are connected to the Scheduler module.
module LteCell
{
	parameters:
		// The number of neighbouring LteCells
		int maxNumberOfNeighbours;
		// The number of MobileStations in the local cell
		int numberOfMobileStations;
		// The radius defining this cell's reach/size
		double radius;
		int currentChannel;
		// X coordinate of the LteCell's center position
		double xPos;
		// Y coordinate of the LteCell's center position
		double yPos;
		// Size of the simulation are in X direction
		double playgroundSizeX;
		// Size of the simulation are in Y direction
		double playgroundSizeY;
		// The period after which assignments of streams to resource blocks 
		// are recomputed.
		double streamSchedPeriod @unit(s);
		// The id of the BaseStation in this LteCell. This id is unique.
		int bsId;
		@display("p=$xPos,$yPos;r=$radius,,black;i=block/network2_l;bgb=330,202");

	gates:
		// For information exchange with neighbouring cells. 
		// <strong>Note:</strong> Gate indices do <em>NOT</em> correspond to 
		// BaseStation ids!
		input fromCell[];
		// For information exchange with neighbouring cells. 
		// <strong>Note:</strong> Gate indices do <em>NOT</em> correspond to 
		// BaseStation ids!
		output toCell[];

    submodules:
        bs: BaseStation {
            parameters:
                xPos = xPos;
                yPos = yPos;
                bsId = bsId;
                currentChannel = currentChannel;
                numberOfMobileStations = numberOfMobileStations;

            gates:
                fromMs[numberOfMobileStations];
                toMs[numberOfMobileStations];

                fromCell[maxNumberOfNeighbours];
                toCell[maxNumberOfNeighbours];

                fromMsMac[numberOfMobileStations];
                toMsMac[numberOfMobileStations];
        }
        ms[numberOfMobileStations]: MobileStation {
            parameters:
                bsId = bsId;
                msId = index;
                initBsXPos = xPos;
                initBsYPos = yPos;
                radius = radius;
                currentChannel = currentChannel;
            gates:
                toMs[numberOfMobileStations];
                fromMs[numberOfMobileStations];
        }
        sched: Scheduler {
            parameters:
                numberOfMobileStations = numberOfMobileStations;
                streamSchedPeriod = streamSchedPeriod;
        }

    connections allowunconnected:
        //connect BaseStation with MobileStations
        for i=0..numberOfMobileStations-1 {
            bs.toMs[i] --> ms[i].fromBs;
            ms[i].toBs --> bs.fromMs[i];
        }

        // Connect MS with each other for D2D communication
        for i=0..numberOfMobileStations-1, for j=0..numberOfMobileStations-1, if i!=j {
            ms[i].toMs[j] --> ms[j].fromMs[i];
        }

        //connect BaseStation with other Cells
        for i=0..(maxNumberOfNeighbours-1) {
            bs.toCell[i] --> toCell[i];
            fromCell++ --> bs.fromCell[i];
        }

        //for direct schedule exchange of MS/BS
        for i=0..numberOfMobileStations-1 {
            bs.toMsMac[i] --> ms[i].fromBsMac;
            ms[i].toBsMac --> bs.fromMsMac[i];
        }
        // Access to the Streaming scheduler module for MS
        for i=0..numberOfMobileStations-1 {
            ms[i].toScheduler --> sched.fromMs[i];
            sched.toMs[i] --> ms[i].fromScheduler;
        }
        // Access to the Streaming scheduler for the BS
        sched.toBs --> bs.fromScheduler;
        bs.toScheduler --> sched.fromBs;
}
