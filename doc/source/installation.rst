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

    $ conda create -n alm python=3.8.13
    $ conda activate alm
    $ pip install pymatgen 
    $ conda install -c conda-forge phonopy
    $ pip install ase
    $ pip install seekpath
    $ pip install custodian
    $ conda install -c conda-forge eigen
    $ conda install -c conda-forge gcc
    $ pip install xmltodict
    $
    $ conda install -c conda-forge mkl
    $
    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib


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


Setting for VASP with pymatgen
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add the following line. In the directory, potpaw_PBE exists.

.. code-block:: bash
    
    $ cat ~/.bash_profile
    
    ...
    export VASP_PP_PATH=$HOME/.../.../Potentials
    ...


Installation of ALM
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash
    
    $ git clone https://github.com/ttadano/ALM.git
    $ cd ./ALM
    $ git pull
    $ cd ./python
    $ python setup.py install

.. For Grand-Chariot, the following line may need to be added in setup.py.
.. 
.. .. code-block:: bash
.. 
..     os.environ["CC"] = /usr/bin/gcc


Modification of .bash_profile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    
You may also need to add .../alm/lib in LD_LIBRARY_PATH.

.. conda-block:: bash
    


