Memory management and optimization
==================================

It aims at reducing the memory consumption.

Background
----------

* Input of the libSGM : a 3-dimension cost volume, size: W, H, D
* Ouput of the libSGM : a 3-dimension cost volume, size : W, H, D

For the SGM computation, temporary data are needed. For each point, we must know the value of previous computed point, depending on the direction.
As a reminder, SGM is performed along 8 directions in our implementation, with 2 passes through the cost volume : a top-down one and a bottom-up one.

Naive solution
--------------

The most intuitive solution, but the least effective one, is to store 3-dimension temporary cost volumes, one per direction.
The size of each one is W, H, D

So, with this solution, libSGM requires :math:`2 \times W \times H \times D + 8 \times W \times H \times D`


Improved solution
------------------

There is no need to store as much information as in the previous solution. It is necessary and sufficient.
to keep information from the previous point (along the depth axis).

Let's take the example of aggregation along to the horizontal left to right path.
To compute a point at (x,y,d) position, we need to access the values of the points calculated in position (x-1,y) along the depth axis.
So, it is necessary to have stored them in a buffer of size D.
Moreover, this buffer is updated at each new point (x,y,d) calculated. That's why we use a temporary value.

The following diagram explains the concept:

    .. image:: ../images/sgm_memory_optim.png
   
We want to compute the point (x,y,d), in orange on the diagram. So, we need to read values previously computed at the point (x-1,y,d-1), (x-1,y,d) and (x-1,y,d+1).
These are stored in the buffer, except the value (x-1,y,d-1) which is stored separately (in purple). Why ? Because in the buffer, at (d-1) position, it is not the value (x-1,y,d-1) anymore which is stored. But this is the value representing (x,y,d-1) that we just computed. Because the buffer is updated in real-time.
This is the explanation for temporary value being (the purple one on the diagram).

At this step, all required previous values are available: the two in red and the one in purple. So, we can compute value for the orange point.
Once this value is computed the information must be stored in the buffer. The value of the orange point is copied in position (d) on the buffer.
Just before, we took care to copy the old value of the buffer into the temporary value (purple).

The total size of temporary information is equal to D+1 for this direction.

For other directiions, memory management also works with buffers. However, it is not necessary to store the information of a point (along the D axis) as for the left->right direction. 
But it is necessary to store the information of a whole line. We then use storage buffers of size W*D. As with the purple point in the previous example, some useful information 
are erased by the buffer updating process.

* For up to down vertical direction: storing an extra point, size = 1
* For diagonal direction from the left corner: storing an extra point along depth axis and another point, size = D +1
* For diagonal direction from right angle: no additional storage

In conclusion, the size of temporary stored data, for the 4 directions of the top-down pass is :math:`3 \times W \times D +2 \times D +2`.
With the other pass, the bottom-up one, the total size is equal to :math:`6 \times W \times D + 4 \times D + 4`.