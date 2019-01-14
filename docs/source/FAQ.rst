
FAQ
---

**1. What other kinds of channel model can be used for the simulation?** 

In addition to the METIS channel model, there are another two channel models (Exponential and Factory channel) which can be used in the simulation environment, to model the channel of communication between mobile stations or mobile station-base station.

**2. How are the channel coefficients calculated in the Exponential channel model**

The channel coefficients are calculated by multiplying the value of the path gain with a random value obtained from an exponential distribution. Moreover, the path gain is calculated taking into consideration the distance between the sender-receiver and the value of the path loss exponent.

**3. How are the channel coefficients calculated in the Factory channel model?**

For the factory channel model, the channel coefficients are calculated as a multiplication of the path gain, fading, and shadowing. 
The path gain is calculated taking into consideration the distance between the sender-receiver, the value of the pathloss at the reference distance and the value of the path loss exponent.
Moreover, fading and shadowing variables have random values obtained from exponential and normal distributions respectively. 

**4. How can the values be assigned to parameters?** 

The parameters can get their values in different ways
- from the configuration (omnet.ini ) file.
- from the NED file.
- from the use.
