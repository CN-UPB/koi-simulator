network ShopFloor  {
	parameters:
		int numCellsX;
		int numCellsY;
		int cellX;
		int cellY;

		int maxNumberOfNeighbours;
		int numberOfMobileStations;
		int numberOfCells = numCellsX*numCellsY;
		double lteFactor;
		double chnBandwidth @unit(Hz);
		double bandwidthPerRB @unit(Hz);
		double realBw = bandwidthPerRB * lteFactor @unit(Hz); 
		int upResourceBlocks = floor(chnBandwidth/realBw/2);
		int downResourceBlocks = floor(chnBandwidth/realBw/2);
		double tti @unit(s);
		double realtti = tti/lteFactor @unit(s);

submodules:
		cell[numberOfCells]: LteCell  {
			parameters:
				numberOfMobileStations = numberOfMobileStations;
				upResourceBlocks = upResourceBlocks;
				downResourceBlocks = downResourceBlocks;
				bandwidthPerRB = realBw;
				bsId = index;
				radius = cellX/2;
				playgroundSizeX = numCellsX*cellX;
				playgroundSizeY = numCellsY*cellY;
				xPos = (index % numCellsX)*cellX+(cellX/2);
				yPos = int(index/numCellsY)*cellY+(cellY/2);
				tti = realtti;
			gates:
				fromCell[maxNumberOfNeighbours];
				toCell[maxNumberOfNeighbours];
			}	
		connections allowunconnected:
				for i=0..maxNumberOfNeighbours-1, for j=0..maxNumberOfNeighbours-1 {
					cell[i].toCell[j] --> cell[j].fromCell[i] if i!=j;
				}
}
