Installation
==============

Requirements
-------------

Aiida-phonopy has been developped under the following conditions.

* VASP
* alamode (alm, anphon, ALM >= 1.4.0)
* phonopy >= 2.15.1
* custodian >= 2022.5.26
* pymatgen v2022.3.24
* ase v3.22.1
* spglib v1.16.5
* seekpath

Preparation
-------------

* Download all the data of phonondb in .../phonondb-20180417/mp-***.

Installation of python libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    $ conda create -n alm python=3.7
    $ pip install pymatgen 
    $ conda install -c conda-forge phonopy
    $ pip install ase
    $ pip install seekpath
    $ pip isntall custodian
    $ conda install -c conda-forge eigen

.. Installation of Eigen
.. ^^^^^^^^^^^^^^^^^^^^^^^
.. 
.. .. code-block:: bash
..     
..     $ cd .../eigen-3.4.0
..     $ mkdir build
..     $ cd ./build
..     $ cmake3 ..
..     $ cmake3 . -DCMAKE_INSTALL_PREFIX=/home/*****/usr/local
..     $ make install
.. 
.. * Check /home/*****/usr/local/include/eigen3


Installation of ALM
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash
    
    $ git clone https://github.com/ttadano/ALM.git
    $ cd ./ALM/python
    $ python setup.py install
    

