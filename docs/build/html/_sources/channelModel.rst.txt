Channel Model
-------------

The channel model implemented in the simulator is the METIS channel model [NKRR14]_. The geometry-based stochastic model have the following approach:

.. figure:: channelModel.png
    :width: 800px
    :align: center
    :height: 500px
    :alt: alternate text
    :figclass: align-center 
    
    METIS procedure for the generation of the channel coefficient [NKRR14]_. 

The following parameters are calculated for the downstream and upstream directions taking into the consideration the type of communication (D2D or communication through the base station) as given in [NKRR14]_.

.. toctree::
   :maxdepth: 2

   GeneralParameters
   smallScale
   coeffGen



