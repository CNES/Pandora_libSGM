Before coding
=============

Thinks and tips to know before coding

About cython
------------

We have decided to distribute the generated .c files as well as our Cython sources.
So, without Cython, one can install libsgm package. 

This means that developers have to maintain the generated .c files. It can be done automatically, by Cython compilation.
To generate them, make sur Cython is installed and use the following command

.. sourcecode:: text

    $ python setup.py build_ext


Unit testing
------------

`Google test framework <https://github.com/google/googletest>`__ is used on this project.
Don't forget to write and maintain unit tests. 

