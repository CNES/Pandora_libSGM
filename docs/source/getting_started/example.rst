Examples
========

SGM applied on random values
----------------------------

.. sourcecode:: python

	import numpy as np
	from libSGM import sgm_wrapper

	if __name__ == '__main__':

	    direction = np.array([[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]], dtype=np.int32)

	    # Create a random cost volume
	    cost_volume_in = 1000*np.random.rand(1000, 1000, 100).astype(np.float32)

	    # Create two penalty arrays corresponding to P1, P2 penalties of sgm algorithm
	    penalty_p1 = np.zeros(np.shape(cost_volume_in), dtype=cost_volume_in.dtype.type) + 8
	    penalty_p2 = np.zeros(np.shape(cost_volume_in), dtype=cost_volume_in.dtype.type) + 32

	    # Sgm computation 
	    cost_volumes_out = sgm_wrapper.sgm_api(cost_volume_in, penalty_p1, penalty_p2, direction, invalid_value=100.0,
	                                           cost_paths=False, overcounting=False)
