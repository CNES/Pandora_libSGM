:py:mod:`libsgm.libsgm_python.sgm_python`
=========================================

.. py:module:: libsgm.libsgm_python.sgm_python`

This module contains functions to execute python SGM library.



Module Contents
---------------

Functions
~~~~~~~~~


.. py:function:: run_sgm(cv_in: np.ndarray,p1_in: np.ndarray,p2_in: np.ndarray,directions: np.ndarray,
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
   :return: the aggregated cost volume and the minimum cost along directions
   :rtype: dict of 3D np.ndarray


.. py:function:: compute_lr(cv_in_2d_front: np.ndarray,lr_2d_previous: np.ndarray,disp: int,
    p1_in_1d: np.ndarray,p2_in_1d: np.ndarray) -> np.ndarray:

   Compute Lr of current plane, at a given disparity

   :param cv_in_2d_front: cost volume
   :type cv_in_2d_front:  np.ndarray
   :param lr_2d_previous: previous lr computed
   :type lr_2d_previous: np.ndarray
   :param d: current disparity
   :type d:  int
   :param p1_in_1d: p1 penalties, dim=1
   :type p1_in_1d: np.ndarray
   :param p2_in_1d: p2 penalties, dim=1
   :type p2_in_1d: np.ndarray
   :return: partial lrs
   :rtype: np.ndarray
