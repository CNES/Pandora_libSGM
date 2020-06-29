## libSGM

C++ implementation of Semi-Global Matching (SGM) algorithm that is wrapped with cython to provide a `libsgm` python module. 

SGM references this work is based upon: 
 - [H. Hirschmuller, "Stereo Processing by Semiglobal Matching and Mutual Information," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 30, no. 2, pp. 328-341, Feb. 2008. doi: 10.1109/TPAMI.2007.1166] 
 - [Ernst, Ines & Hirschmüller, Heiko. (2008). Mutual Information Based Semi-Global Stereo Matching on the GPU. Proceedings of the International Symposium on Visual Computing. 5358. 10.1007/978-3-540-89639-5_22.]
- [Hirschmüller, Heiko & Buder, Maximilian & Ernst, Ines. (2012). Memory Efficient Semi-Global Matching. ISPRS Annals of Photogrammetry, Remote Sensing and Spatial Information Sciences. I-3. 10.5194/isprsannals-I-3-371-2012.]

## Dependencies

Dependencies are 

    python, numpy, gcc
    
For developers, there is another one

    cython

## Installation
First of all, make sure numpy is already installed.

    pip install numpy 

This package can be install through the following commands:

    cd libsgm
    pip install .

For developers,before pip command, do not forget to delete sources/sgm_wrapper.cpp file
If not, pip won't use cython to generate a new .cpp file. It will use the last updated .cpp file on git repository. This may cause conflicts. 

**Build documentation**

Make sure Doxygen, latex and dvipng is already installed

```
pip install exhale
pip install sphinx-rtd-theme
python setup.py build_sphinx
```

Documentation is built in libsgm/doc/build/html

## Usage

libSGM must be used as a package : 

    from libSGM import sgm_wrapper

## Things to know

* **Input Cost volume type**


Input Cost volume type for sgm_wrapper function can be *float* or *uint8*. 
Actually, libSGM would be able to receive any type thanks to "template" methods. You just need to add the type you want in sgm_wrapper.pyx on `ctypedef fused my_cv_type` and instantiate function at the end of lib/sgm.cpp . 

## References

If you use this CNES software, please cite the following paper: 

Cournet, M., Sarrazin, E., Dumas, L., Michel, J., Guinet, J., Youssefi, D., Defonte, V., Fardet, Q., 2020. 
Ground-truth generation and disparity estimation for optical satellite imagery.
ISPRS - International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences.
