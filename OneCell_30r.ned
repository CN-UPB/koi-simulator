//Namespace(c=1, chnIni=True, lp=1, n='OneCell_30r', ned=True, parIni=True, r=30.0, x=60.0, y=60.0)
network OneCell_30r  {
    parameters:
	double playgroundSizeX = 60.0;
	double playgroundSizeY = 60.0;
	double radius = 30.0;

	//network parameters
	int numberOfCells = 1;

	int maxNumberOfNeighbours;
	int numberOfMobileStations;

        cell[*].ms[*].initQuadrant = intuniform(0, 5);
        cell[*].ms[*].initPosAlpha = uniform(0, 1);
        cell[*].ms[*].initPosBeta = uniform(0, 1);
        cell[*].ms[*].initPosGamma = uniform(0, 1);
	cell[0].xPos = 31;
	cell[0].yPos = 31;
	@display("bgb=$playgroundSizeX,$playgroundSizeY");

    submodules:
   watch: StopWatch;
	cell[numberOfCells]: LteCell  {
		parameters:
			numberOfMobileStations = numberOfMobileStations;
			bsId = index;
			radius = radius;
		gates:
			fromCell[maxNumberOfNeighbours];
			toCell[maxNumberOfNeighbours];

	}	

    connections allowunconnected:
}
