Small scale parameters.
^^^^^^^^^^^^^^^^^^^^^^^
The calculation of the small scale parameters is done taking into consieration the scenario (LOS/NLOS), in which the communication takes place. The following steps and corresponding formulas are given in [NKRR14]_.

- Firstly, the delays for each cluster are generated. 

- Secondly, the power for each cluster is calculated taking into consideration the cluster delays.

- Calculate the azimuth and zenith angles of arrival and departures.

- For each pair of sender/receiver, random phases, uniformly distributed are added to each ray of each cluster. The number of clusters varies based on the scenario.

- For each pair of sender/receiver, cross-Polarization values are added to each ray of each cluster. The number of clusters varies based on the scenario.

