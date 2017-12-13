//Namespace(c=1, chnIni=True, lp=1, n='Debug', ned=True, parIni=True, r=30.0, x=100.0, y=100.0)
network Verification  {
    parameters:
	double playgroundSizeX = 100.0;
	double playgroundSizeY = 100.0;
	double radius = 20.0;

	//network parameters
	int numberOfCells = 2;

	int maxNumberOfNeighbours;
	int numberOfMobileStations;

        cell[*].ms[*].initQuadrant = intuniform(0, 5);
        cell[*].ms[*].initPosAlpha = uniform(0, 1);
        cell[*].ms[*].initPosBeta = uniform(0, 1);
        cell[*].ms[*].initPosGamma = uniform(0, 1);
	@display("bgb=$playgroundSizeX,$playgroundSizeY");

    submodules:
	cell[numberOfCells]: LteCell  {
		parameters:
			numberOfMobileStations = numberOfMobileStations;
			bsId = index;
			radius = radius;
			playgroundSizeX = playgroundSizeX;
			playgroundSizeY = playgroundSizeY;
		gates:
			fromCell[maxNumberOfNeighbours];
			toCell[maxNumberOfNeighbours];

	}	

    connections allowunconnected:
        cell[0].toCell[0] --> cell[1].fromCell[0];
        cell[1].toCell[0] --> cell[0].fromCell[1];
}
