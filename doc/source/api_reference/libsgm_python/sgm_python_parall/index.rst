:py:mod:`libsgm.libsgm_python.sgm_python_parall`
================================================

.. py:module:: libsgm.libsgm_python.sgm_python_parall`

This module contains functions to execute python SGM library.



Module Contents
---------------

Functions
~~~~~~~~~

.. py:function:: run_sgm_parall(cv_in: np.ndarray,p1_in: np.ndarray,p2_in: np.ndarray,directions: np.ndarray,
    cost_paths: bool = False,overcounting: bool = False) -> Dict:

   Run Python LibSGM

   :param cv_in: cost volume
   :type cv_in:  np.ndarray
   :param p1_in: p1 penalties
   :type p1_in:  3D np.ndarray
   :param p2_in: p2 penalties
   :type p2_in:  3D np.ndarray
   :param directions: directions used in SGM
   :type directions: 2D np.ndarray
   :param cost_paths: True if Cost Volumes along direction are to be returned
   :type cost_paths: bool
   :param overcounting: over-counting correction option
   :type overcounting: bool
   :return: the aggregated cost volume and the minimum cost along 8 directions
   :rtype: dict of 3D np.ndarray


.. py:function:: compute_costs(starting_points: List,cv_in: np.ndarray,p1_in: np.ndarray,p2_in: np.ndarray,
    dir_str: np.ndarray,cost_paths: bool = False) -> Dict:

   Compute cost volume, starting from given points / directions

   :param starting_points: List of starting points [i0, j0, dir_i, dir_j, num_dir]
   :type starting_points:  List
   :param cv_in: cost volume
   :type cv_in:  3D np.ndarray
   :param p1_in: p1 penalties
   :type p1_in:  3D np.ndarray
   :param p2_in: p2 penalties
   :type p2_in:  3D np.ndarray
   :param dir_str: list of direction names
   :type dir_str:  List
   :param cost_paths: True if Cost Volumes along direction are to be returned
   :type cost_paths: bool
   :return: the aggregated cost volume and the minimum cost along directions
   :rtype: dict of 3D np.ndarray