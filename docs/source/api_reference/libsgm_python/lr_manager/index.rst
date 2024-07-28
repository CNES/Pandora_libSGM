:py:mod:`libsgm.libsgm_python.lr_manager`
=========================================

.. py:module:: libsgm.libsgm_python.lr_manager

This module contains classes and functions to manage parallel planes indices moving along given direction.



Module Contents
---------------

Classes
~~~~~~~


.. py:class:: LrManager

   LrManager class

   .. py:attribute:: cv_shape
      :type: List



   .. py:attribute:: direction
      :type: Tuple



   .. py:function:: next(self) -> None:

      Update the state of the planes, moves forward

      :return: None


   .. :py:function:: set_current_lr(self, lr_s: np.ndarray) -> None:

      Set current lr

      :param lr_s: partial cost lr
      :type lr_s:  np.ndarray
      :return: None


   .. py:function:: get_previous_lr(self, num_plane: int) -> np.ndarray:

      Get previous lr

      :param num_plane: numero of plane
      :type num_plane:  int
      :return: previous partial cost lr
      :rtype: np.ndarray
