The KoI Simulator
=================

This package contains the source code and simulation parameters for the KOI 
scheduling simulations.

For installation/build instructions, please consult the INSTALL file.

To execute an example simulation, build the simulation code as described in 
the INSTALL file. Then, have a look into the 
simulations folder for 
configurations.

Generating the Documentation
----------------------------

To generate the documentation, you need the software Doxygen, to generate 
the C++ API docs. Additionally, you will need a working version of the 
Omnet++ Eclipse IDE. First, change to the directory where you downloaded 
the simulation software and run
	
	$ doxygen Doxyfile

This command will place the C++ API docs in the directory `doc/doxy`.
To generate the documentation for the NED files, open the Omnet++ IDE 
with a valid workspace. Then, import the KOI Simulator project by 
choosing *Import...* from the *File* menu, then 
*Existing projects into Workspace* as the import source. 

Then, on the *Import Projects* mask, click *Browse* and navigate to the 
directory where you stored the simulator. Afterwards, a project named 
*abstractLte* should show up under *Projects*. Select it and click *finish*.

After you've imported the project successfully, right click it in your project 
explorer and choose *Generate NED Documentation*. Now, you will find 
the full documentation under the `doc/neddoc` directory. To read the 
freshly generated documentation, point your browser to `doc/index.html`.

The documentation can be also found in [Read the Docs.](https://koi-documentation.readthedocs.io/en/latest/)

History
--------

This simulation was developed in the [KoI project](http://www.koi-projekt.de). 
It is partially based on older simulation software, coming out of RWTH Aachen, TU Berlin, 
KTH Stockholm, and Paderborn University. A grateful hattip to all the work that preceeded this one! 

This simulator was used in and developed in the context of the followin papers (among others): 

M. Bohge, J. Gross, and A. Wolisz
“Optimal Soft Frequency Reuse and Dynamic Sub-carrier Assignments in Cellular OFDMA Networks” 
Transactions on Emerging Telecommunications Technologies, vol. 21, no. 8, pp. 704–713, 2010.

D. Parruca and J. Gross
“Throughput Analysis of Proportional Fair Scheduling For Sparse and Ultra-Dense Interference-Limited OFDMA/LTE Networks”
IEEE Transactions on Wireless Communications, vol. 15, no. 10, pp. 6857-6870, 2016.

If you use this simulator in other project, we'd appreciate if you'd let us know. 
