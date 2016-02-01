//Namespace(c=4, chnIni=True, lp=4, n='donald_7cells', ned=True, parIni=True, r=250.0, x=6000.0, y=3464.0)
network donald_7cells  {
    parameters:
	double playgroundSizeX = 1000.0;
	double playgroundSizeY = 1000.0;
	double radius = 250.0;

	//network parameters
	int numberOfCells = 7;

	int maxNumberOfNeighbours;
	int numberOfMobileStations;

	cell[*].ms[*].initQuadrant = intuniform(0, 5);
	cell[*].ms[*].initPosAlpha = uniform(0, 1);
	cell[*].ms[*].initPosBeta = uniform(0, 1);
	cell[*].ms[*].initPosGamma = uniform(0, 1);
	cell[0].xPos = 500.0;
	cell[0].yPos = 1000.0;
	cell[1].xPos = 500.0;
	cell[1].yPos = 0.0;
	
	cell[2].xPos = 933.0;
	cell[2].yPos = 750.0;
	cell[3].xPos = 500.0;
	cell[3].yPos = 500.0;
	cell[4].xPos = 67.0;
	cell[4].yPos = 750.0;

	cell[5].xPos = 67.0;
	cell[5].yPos = 250.0;
	cell[6].xPos = 933.0;
	cell[6].yPos = 250.0;
        
	@display("bgb=$playgroundSizeX,$playgroundSizeY");

    submodules:
   watch: StopWatch;
	cell[numberOfCells]: LteCell {
		parameters:
			numberOfMobileStations = numberOfMobileStations;
			bsId = index;
			radius = radius;
		gates:
			fromCell[maxNumberOfNeighbours];
			toCell[maxNumberOfNeighbours];
	}

    connections allowunconnected:
    
        cell[0].toCell[0] --> cell[1].fromCell[0];
        cell[0].toCell[1] --> cell[2].fromCell[0];
        cell[0].toCell[2] --> cell[3].fromCell[0];
        cell[0].toCell[3] --> cell[4].fromCell[0];
        cell[0].toCell[4] --> cell[5].fromCell[0];
        cell[0].toCell[5] --> cell[6].fromCell[0];
        cell[1].toCell[0] --> cell[0].fromCell[1];
        cell[1].toCell[1] --> cell[2].fromCell[1];
        cell[1].toCell[2] --> cell[3].fromCell[1];
        cell[1].toCell[3] --> cell[4].fromCell[1];
        cell[1].toCell[4] --> cell[5].fromCell[1];
        cell[1].toCell[5] --> cell[6].fromCell[1];
        cell[2].toCell[0] --> cell[0].fromCell[2];
        cell[2].toCell[1] --> cell[1].fromCell[2];
        cell[2].toCell[2] --> cell[3].fromCell[2];
        cell[2].toCell[3] --> cell[4].fromCell[2];
        cell[2].toCell[4] --> cell[5].fromCell[2];
        cell[2].toCell[5] --> cell[6].fromCell[2];
        cell[3].toCell[0] --> cell[0].fromCell[3];
        cell[3].toCell[1] --> cell[1].fromCell[3];
        cell[3].toCell[2] --> cell[2].fromCell[3];
        cell[3].toCell[3] --> cell[4].fromCell[3];
        cell[3].toCell[4] --> cell[5].fromCell[3];
        cell[3].toCell[5] --> cell[6].fromCell[3];
        cell[4].toCell[0] --> cell[0].fromCell[4];
        cell[4].toCell[1] --> cell[1].fromCell[4];
        cell[4].toCell[2] --> cell[2].fromCell[4];
        cell[4].toCell[3] --> cell[3].fromCell[4];
        cell[4].toCell[4] --> cell[5].fromCell[4];
        cell[4].toCell[5] --> cell[6].fromCell[4];
        cell[5].toCell[0] --> cell[0].fromCell[5];
        cell[5].toCell[1] --> cell[1].fromCell[5];
        cell[5].toCell[2] --> cell[2].fromCell[5];
        cell[5].toCell[3] --> cell[3].fromCell[5];
        cell[5].toCell[4] --> cell[4].fromCell[5];
        cell[5].toCell[5] --> cell[6].fromCell[5];
        cell[6].toCell[0] --> cell[0].fromCell[6];
        cell[6].toCell[1] --> cell[1].fromCell[6];
        cell[6].toCell[2] --> cell[2].fromCell[6];
        cell[6].toCell[3] --> cell[3].fromCell[6];
        cell[6].toCell[4] --> cell[4].fromCell[6];
        cell[6].toCell[5] --> cell[5].fromCell[6];
}
