#!/usr/bin/env python
# coding: utf8
#
# Copyright (c) 2024 Centre National d'Etudes Spatiales (CNES).
#
# This file is part of LIBSGM
#
#     https://github.com/CNES/Pandora_libsgm
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""
This module contains functions to execute SGM library.
"""

cimport cython
from cpython cimport Py_INCREF, PyObject
from libc.stdint cimport *
from libc.stdlib cimport *
from libcpp cimport bool
from libcpp.string cimport string

import numpy as np

cimport numpy as np
cimport openmp

#Numpy initialization to avoid segfaults
np.import_array()

#Fused type to allow multiple types definitions
#The following list represent Cost volume types avalaible
#You can add another one
ctypedef fused my_cv_type:
    uint8_t
    float

# ArrayWrapper class  to deallocate our array
# when Python Object is deleted

cdef class Array3DWrapper:

    cdef void* data_ptr
    cdef int nb_row , nb_col , nb_disp

    cdef set_data(self, int nb_row, int nb_col , int nb_disp, void * data_ptr):
        self.data_ptr = data_ptr
        self.nb_row = nb_row
        self.nb_col = nb_col
        self.nb_disp = nb_disp


    def __array__(self):
        """
        __Array__ methode called when numpy tries to
            get an array from the object
        """

        cdef np.npy_intp shape[3]
        shape[0] = <np.npy_intp> self.nb_row
        shape[1] = <np.npy_intp> self.nb_col
        shape[2] = <np.npy_intp> self.nb_disp

        #Create a 3D array pf length (size_row,size_col,size_depth)
        ndarray = np.PyArray_SimpleNewFromData(3,shape,np.NPY_FLOAT32, self.data_ptr)
        return ndarray


    def __dealloc__(self):
        """
        Free the array. This is called by Python when all the references
        to the object are gone
        """

        free(<void*>self.data_ptr)


cdef class Array3DWrapperINT:


    cdef void* data_ptr
    cdef int nb_row , nb_col , nb_disp

    cdef set_data(self, int nb_row, int nb_col , int nb_disp, void * data_ptr):
        self.data_ptr = data_ptr
        self.nb_row = nb_row
        self.nb_col = nb_col
        self.nb_disp = nb_disp


    def __array__(self):
        """
        __Array__ methode called when numpy tries to
        get an array from the object
        """

        cdef np.npy_intp shape[3]
        shape[0] = <np.npy_intp> self.nb_row
        shape[1] = <np.npy_intp> self.nb_col
        shape[2] = <np.npy_intp> self.nb_disp

        #Create a 3D array pf length (size_row,size_col,size_depth)
        ndarray = np.PyArray_SimpleNewFromData(3,shape,np.NPY_INT32, self.data_ptr)
        return ndarray


    def __dealloc__(self):
        """
        Free the array. This is called by Python when all the references
        to the object are gone
        """

        free(<void*>self.data_ptr)

# Prototype of C function
cdef extern from "lib/sgm.hpp":

    ctypedef struct CostVolumes:
        float * cost_volume
        int * cost_volume_min



    CostVolumes sgm[T](T * cv_in, T * p1_in,T * p2_in , int* directions, unsigned long int nb_row,
        unsigned long int nb_col, unsigned int nb_disp, T invalid_value, float * segmentation, bool cost_paths, bool overcounting)

#Create this function is necessary. The goal is to make fused type and memory view work together
#Avoid error as " Cannot coerce to a type that is not specialized"
def sgm_call(my_cv_type[:,:,::1] cv_in_memview, my_cv_type[:,:,::1] p1_in_memview, my_cv_type[:,:,::1] p2_in_memview, \
         int[:,::1] directions, my_cv_type invalid_value, float[:,::1] segmentation_memview, bool cost_paths, \
                                                                                           bool overcounting):
    """
    Internal function. Wrap the call of sgm C function.

    :param cv_in_memview: cost volume memory view
    :type cv:  3D array contiguous memoryview
    :param p1_in_memview: p1 penalties
    :type p1_in_memview:  3D array contiguous memoryview
    :param p2_in_memview: p2 penalties
    :type p2_in_memview:  3D array contiguous memoryview
    :param directions: directions used in SGM
    :type directions: 2D array
    :param invalid value: value representing invalid cost
    :type invalid_value: uint8 or float
    :param segmentation_memview: segmentation
    :type segmentation_memview:  2D array contiguous memoryview
    :param cost_paths: True if Cost Volumes along direction are to be returned
    :type cost_paths: bool
    :param overcounting: over-counting correction option
    :type overcounting: bool
    :return: the aggregated cost volume and the minimum cost along 8 directions
    :rtype: dict of 3D arrays
    """

    nb_row, nb_col, nb_disp = cv_in_memview.shape[0], cv_in_memview.shape[1], cv_in_memview.shape[2]
    nb_dir = 8

    #Float pointer for
    cdef CostVolumes cv_out

    #Calling C function
    cv_out = sgm(&cv_in_memview[0,0,0], &p1_in_memview[0,0,0], &p2_in_memview[0,0,0], &directions[0,0], nb_row, \
                                    nb_col,  nb_disp, invalid_value, &segmentation_memview[0,0], cost_paths, \
                                                                      overcounting)

    # Build Python object which is an numpy array, data coming from the Cn pointer
    cdef np.ndarray ndarray_out_cv
    array_wrapper_cv = Array3DWrapper()
    array_wrapper_cv.set_data(nb_row, nb_col, nb_disp, < void * > cv_out.cost_volume)
    ndarray_out_cv = np.array(array_wrapper_cv, copy=False)
    # Assign our object to the base of the ndarray object
    np.PyArray_SetBaseObject(ndarray_out_cv, array_wrapper_cv)
    # Increment the reference count, as the above assignement was done in C
    # and Python does not know that there is this additional reference
    Py_INCREF(array_wrapper_cv)


    # Build Python object which is an numpy array, data coming from the Cn pointer
    # following variables are defined outside 'if' because of cython constraints.
    cdef np.ndarray ndarray_out_cv_min
    array_wrapper_cv_min = Array3DWrapperINT()
    if (cost_paths):
        array_wrapper_cv_min.set_data(nb_row, nb_col, nb_dir, < void * > cv_out.cost_volume_min)
        ndarray_out_cv_min = np.array(array_wrapper_cv_min, copy=False)
        # Assign our object to the base of the ndarray object
        np.PyArray_SetBaseObject(ndarray_out_cv_min, array_wrapper_cv_min)
        # Increment the reference count, as the above assignement was done in C
        # and Python does not know that there is this additional reference
        Py_INCREF(array_wrapper_cv_min)
        cv = {'cv' : ndarray_out_cv, 'cv_min' : ndarray_out_cv_min}
    else :
        cv = {'cv' : ndarray_out_cv, 'cv_min' : None}

    return cv


def sgm_api( cv_in not None, p1_in not None, p2_in not None, direction not None, invalid_value, segmentation not None, \
                                                                       cost_paths = False, overcounting = False):
    """
    Compute aggregated cost volume using C++ library where sgm method is implemented

    :param cv: cost volume without aggregation
    :type cv:  3D array contiguous
    :param p1_in: p1 penalties
    :type p1_in:  3D array contiguous
    :param p2_in: p2 penalties
    :type p2_in:  3D array contiguous
    :param directions: directions used in SGM
    :type directions: 2D array
    :param invalid value: value representing invalid cost
    :type invalid value: uint8 or float
    :param segmentation: segmentation
    :type segmentation:  2D array contiguous
    :param cost_paths: True if Cost Volumes along direction are to be returned
    :type cost_paths: bool
    :param overcounting: over-counting correction option
    :type overcounting: bool
    :return: the aggregated cost volume ('cv') and the minimum cost along 8 directions ('cv_min')
    :rtype: dict of 3D arrays
    """

    #Calling wrapper
    ndarray_out = sgm_call(cv_in, p1_in, p2_in, direction,  invalid_value, segmentation, cost_paths, overcounting)

    return ndarray_out
