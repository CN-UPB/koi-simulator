import sys
import argparse

parser = argparse.ArgumentParser(description='Creates a ned file for a specific network of the abstractLTE model.')
parser.add_argument('-x', required=True, type=float,
                   help='the x value of the playground size for the hole network')
parser.add_argument('-y', required=True, type=float,
                   help='the y value of the playground size for the hole network')
parser.add_argument('-r', required=True, type=float,
                   help='the radius of one LTE cell')
parser.add_argument('-lp', required=True, type=int,
                   help='the number of logical partitions for the cell distribution')
parser.add_argument('-n', required=True, type=str,
                   help='name of the network also for the output files')
parser.add_argument('-c', required=True, type=int,
                   help='number of channels')
parser.add_argument('-ned', required=False, type=bool, default=True,
                   help='generate the .ned file for the network')
parser.add_argument('-parIni', required=False, type=bool, default=True,
                   help='generate the partition .ini file for the network')
parser.add_argument('-chnIni', required=False, type=bool, default=True,
                   help='generate the channel .ini file for the network')

args = parser.parse_args()

for x in range(100000, 160000, 1000):
    for y in range(100000, 160000, 1000):
        networkName = args.n
        playgroundSizeX = x
        playgroundSizeY = y
        channels = args.c
        radius = args.r
        lp = args.lp

        bsPos = []
        bsChn = []
        #rows without offsets
        dX = 3 * radius
        dY = (3**0.5) * radius
        xOffset = radius
        row = 0
        column = 0
        cellCount = 0
        while(1):
            xPos = xOffset + dX * column
            if(xPos > (playgroundSizeX - radius)):
                column = 0
                xPos = xOffset + dX * column
                row += 1
            yPos = radius + dY * row
            if(yPos >= playgroundSizeY - radius):
                break
            bsPos.append((xPos, yPos))
            if(cellCount%2 == 0):
                tmpChn = row % channels
                bsChn.append(tmpChn)
            else:
                tmpChn = (row+1) % channels
                bsChn.append(tmpChn)
            column += 1
            cellCount += 1

        #rows with offsets
        dX = 3 * radius
        dY = (3**0.5) * radius
        xOffset = 2.5 * radius
        yOffset = (3**0.5)/2 * radius + radius
        row = 0
        column = 0
        while(1):
            xPos = xOffset + dX * column
            if(xPos > (playgroundSizeX - radius)):
                column = 0
                xPos = xOffset + dX * column
                row += 1
            yPos = yOffset + dY * row
            if(yPos >= playgroundSizeY - radius):
                break
            column += 1
            cellCount += 1
            bsPos.append((xPos, yPos))
            tmpChn = (row+3) % channels
            bsChn.append(tmpChn)

        numberOfCells = cellCount
        if(numberOfCells == 2*3072):
            print "Number of cells: " + str(int(numberOfCells))
            print "x: " + str(int(x))
            print "y: " + str(int(y))
