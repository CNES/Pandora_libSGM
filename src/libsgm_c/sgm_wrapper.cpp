// Copyright (C) 2024 Centre National d'Etudes Spatiales (CNES).
//
// This file is part of cars-filter
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "sgm.cpp"




namespace py = pybind11;

template<typename T, typename Tout>
py::dict pySgmApi(py::array_t<T, py::array::c_style> cv_in,
                 py::array_t<T, py::array::c_style> p1_in,
                 py::array_t<T, py::array::c_style> p2_in,
                 py::array_t<int, py::array::c_style> directions,
                 float invalid_value,
                 py::array_t<float, py::array::c_style> segmentation,
                 bool cost_paths,
                 bool overcounting)
{

    auto cv_in_shape = cv_in.shape();
    auto p1_in_shape = p1_in.shape();
    auto p2_in_shape = p2_in.shape();
    auto segmentation_shape = segmentation.shape();
    auto directions_shape = directions.shape();

    unsigned long int nb_rows = cv_in_shape[0];
    unsigned long int nb_cols = cv_in_shape[1];
    unsigned int nb_disps = cv_in_shape[2];
    unsigned int nb_directions = 8;

    //Check dimensions
    if (cv_in.ndim() != 3) {
        throw std::invalid_argument("cv_in must be a 3D array.");
    }
    if (p1_in.ndim() != 3) {
        throw std::invalid_argument("p1_in must be a 3D array.");
    }
    if (p2_in.ndim() != 3) {
        throw std::invalid_argument("p2_in must be a 3D array.");
    }
    if (directions.ndim() != 2) {
        throw std::invalid_argument("direction must be a 2D array.");
    }
    if (segmentation.ndim() != 2) {
        throw std::invalid_argument("segmentation must be a 2D array.");
    }


    if (p1_in_shape[0] != cv_in_shape[0] || p1_in_shape[1] != cv_in_shape[1]) {
        throw std::invalid_argument("p1_in dimensions must match cv_in dimensions.");
    }
    if (p2_in_shape[0] != cv_in_shape[0] || p2_in_shape[1] != cv_in_shape[1]){
        throw std::invalid_argument("p2_in dimensions must match cv_in dimensions.");
    }
    if (segmentation_shape[0] != cv_in_shape[0] || segmentation_shape[1] != cv_in_shape[1]) {
        throw std::invalid_argument("segmentation dimensions must match the height and width of cv_in.");
    }
    if (directions_shape[0] != nb_directions) {
        throw std::invalid_argument("SGM only support 8 dimensions");
    }

    /* Request buffers descriptor from Python */
    T* cv_in_buf = const_cast<T*>(cv_in.data());
    T* p1_in_buf = const_cast<T*>(p1_in.data());
    T* p2_in_buf = const_cast<T*>(p2_in.data());
    int* directions_buf = const_cast<int*>(directions.data());
    float* segmentation_buf = const_cast<float*>(segmentation.data());

    CostVolumes<Tout> cv_out = sgm<T, Tout>(
        cv_in_buf,
        p1_in_buf,
        p2_in_buf,
        directions_buf,
        nb_rows,
        nb_cols,
        nb_disps,
        invalid_value,
        segmentation_buf,
        cost_paths,
        overcounting
    );

    py::dict result;
    result["cv"] = py::array_t<Tout>(std::vector<size_t>{nb_rows, nb_cols, nb_disps}, cv_out.cost_volume);
    if (cost_paths) {
        result["cv_min"] = py::array_t<int>(std::vector<size_t>{nb_rows, nb_cols, 8}, cv_out.cost_volume_min);
    } else {
        result["cv_min"] = py::array_t<int>();  // Return an empty array if cost_paths is false
    }
    delete cv_out.cost_volume;
    delete cv_out.cost_volume_min;
    return result;
}

//wrap as Python module
PYBIND11_MODULE(c_libsgm, m)
{
  m.doc() = "Wrapper module of sgm.cpp, a c++ module for sgm";
  m.def("sgm_api",
        &pySgmApi<uint8_t, uint16_t>,
        "Compute aggregated cost volume following Semi-Global algorithm by Hirschmuller",
        py::arg("cv_in").noconvert(), // not allow py:array to convert this arg
        py::arg("p1_in"), // next argument can be casted
        py::arg("p2_in"),
        py::arg("directions"),
        py::arg("invalid_value"),
        py::arg("segmentation"),
        py::arg("cost_paths") = false,
        py::arg("overcounting") = false,
        R"pbdoc(
            Python SGM wrapper

            :param cv_in: Input cost volume
            :type cv_in: uint8_t numpy ndarray
            :param p1_in: p1 matrix
            :type p1_in: uint8_t numpy ndarray
            :param p2_in: p2 matrix
            :type p2_in: uint8_t numpy ndarray
            :param directions: directions to explore
            :type directions: uint8_t numpy ndarray
            :param invalid_value: invalid value to use
            :type invalid_value: uint8_t
            :param segmentation: segmentation matrix
            :type segmentation: uint8_t numpy ndarray
            :param cost_paths: activate cost paths
            :type cost_paths: bool
            :param overcounting: activate overcounting
            :type overcounting: bool
            :return: ("cv": optimize cost volume, "cv_min": cost paths)
            :rtype: dict
        )pbdoc"
  );
  m.def("sgm_api",
        &pySgmApi<float, float>,
        "Compute aggregated cost volume following Semi-Global algorithm by Hirschmuller",
        py::arg("cv_in").noconvert(),  // not allow py:array to convert this arg
        py::arg("p1_in"), // next argument can be casted
        py::arg("p2_in"),
        py::arg("directions"),
        py::arg("invalid_value"),
        py::arg("segmentation"),
        py::arg("cost_paths") = false,
        py::arg("overcounting") = false,
        R"pbdoc(
            Python SGM wrapper

            :param cv_in: Input cost volume
            :type cv_in: float32 numpy ndarray
            :param p1_in: p1 matrix
            :type p1_in: float32 numpy ndarray
            :param p2_in: p2 matrix
            :type p2_in: float32 numpy ndarray
            :param directions: directions to explore
            :type directions: uint8_t numpy ndarray
            :param invalid_value: invalid value to use
            :type invalid_value: float32
            :param segmentation: segmentation matrix
            :type segmentation: uint8_t numpy ndarray
            :param cost_paths: activate cost paths
            :type cost_paths: bool
            :param overcounting: activate overcounting
            :type overcounting: bool
            :return: ("cv": optimize cost volume, "cv_min": cost paths)
            :rtype: dict
        )pbdoc"
  );
}