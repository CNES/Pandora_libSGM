Python SGM, Points management
==================================

LrManager manages the points used with respect to a direction, at a given state.
This class is used in the python version of SGM


Background
----------

* Input of the libSGM : a 3-dimension cost volume, size: W, H, D

For the SGM computation, for the direction d, we can see the sgm problem as a propagation by planes.


Solution
------------------

The first solution is to initialize the LrManager for each direction.
For a direction d, one or several planes are initialized, corresponding to the "origin" of the propagation.
The cost on each point of the planes are computed with the cost values from the previous planes.

The second solution is to initialize all the LrManagers, split all the points associated to its direction, and propagate the cost from every starting points in parallel.


The LrManager is used for the two methods. In both cases, the .__init__() is called, creating the initial planes.
The .next() method is only called by the parallel method. This function update the planes, "moving" each plane  with the wanted direction.


The following diagram represents a LrManager for the direction (1, 1):

    .. image:: ../images/lr_manager.png
