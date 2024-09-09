Examples
========

SGM applied on random values
----------------------------

.. sourcecode:: python

  import numpy as np 
  from libSGM import sgm_wrapper

  direction = np.array([[0, 1], [1, 0], [1, 1], [1, -1], [0, -1], [-1, 0], [-1, -1], [-1, 1]], dtype=np.int32)

  # Create a random cost volume
  cost_volume_in = 1000*np.random.rand(1000, 1000, 100).astype(np.float32)

  # Create two penalty arrays corresponding to P1, P2 penalties of sgm algorithm
  penalty_p1 = np.zeros(np.shape(cost_volume_in), dtype=cost_volume_in.dtype.type) + 8
  penalty_p2 = np.zeros(np.shape(cost_volume_in), dtype=cost_volume_in.dtype.type) + 32

  # Default optimization layer in pandora plugin, for a piecewise optimization layer array use 3sgm method
  # use image shape compatible with cost_volume_in
  optimization_layer = np.ones((1000,1000), dtype=np.float32)

  # Sgm computation 
  cost_volume_out = sgm_wrapper.sgm_api(cost_volume_in, penalty_p1, penalty_p2, direction, invalid_value=100.0, segmentation=optimization_layer, cost_paths=False, overcounting=False)

  # Show Cost Volume 
  print(cost_volume_out["cv"])
