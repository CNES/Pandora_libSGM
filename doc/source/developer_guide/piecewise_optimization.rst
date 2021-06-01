Piecewise optimization optimization
===================================

If the user provides a segmentation composed of more than one segment, a piecewise optimization is applied.

Background
----------

* Input of the libSGM : a 3-dimension cost volume, size: W, H, D, a 2-dimension segmentation, size: W, H
* Ouput of the libSGM : a 3-dimension cost volume, size : W, H, D

For the SGM computation, temporary data are needed. For each point, we must know the value of previous computed point, depending on the direction.
As a reminder, SGM is performed along 8 directions in our implementation, with 2 passes through the cost volume : a top-down one and a bottom-up one.


Solution
--------

For each direction d, when we compute the optimized cost of a pixel p, if p belongs to a different segment than the previous point, we do not use the history of the direction.

A reset factor is computed for the first disp, and is applied to history as follows :

:math:`Lr(p, d) = C(p, d)  + ResetFactor(p, p-r) * [\min_{d'}{Lr(p-r, d') + V(d, d'))} - \min_{d'}{Lr(p-r, d'))}]`
with :math:`p` the current pixel, :math:`d` the current disparity, :math:`r` the current direction, and :math:`V()` the penalty function.


Buffers of previous segments are updated at each pixel.

The following diagram explains the concept:

    .. image:: ../images/piecewise_optimization.png
   
