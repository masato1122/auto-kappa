==============
Installation
==============

Requirements
=============

Auto-kappa requires the following packages. 
While the python packages can be automatically downloaded during the installation process,
users need to separately install VASP and Alamode.

* VASP
* alamode (alm, anphon, ALM >= 1.4.0)
* phonopy 
* custodian 
* pymatgen
* ase
* spglib
* seekpath


Preparation
============

VASP and Alamode
-------------------

Auto-kappa needs ``vasp``, ``alm``, and ``anphon`` commands.
Please install VASP and 
`ALAMODE <https://alamode.readthedocs.io/en/latest/index.html>`_
in advance.

* To allow ASE to call POTCAR files of VASP, 
  set ``VASP_PP_PATH`` variable in your shell configuration file, 
  which may be ``~/.bash_profile``. 
  See the `ASE page <https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html>`_ for detail.

.. code-block:: bash
    
    export VASP_PP_PATH=$HOME/(directory in which potpaw_PBE is located)
    
``potpaw_PBE.54.tar.gz`` supposed to be expanded in ``potpaw_PBE``.


Phonondb
---------

* Download all the data of phonondb with wget.

* You can also find all the data (phonondb-20180417.tar.gz) in Box: 
  `.../Box/ppdb1/Others <https://app.box.com/s/69nioqnpu6xxis5q4f4ua3sqxwwvla36>`_.

* Expand all the tar.gz file in a directory. The name of expanded directories 
  in which ``mp-***`` directories are located.

.. code-block:: bash    
    
    $ cd .../phonondb (arbitrary directory name)
    $ wget http://phonondb.mtl.kyoto-u.ac.jp/_downloads/mp-25-20180417.tar.lzma 
    $ tar xvf mp-25-20180417.tar.lzma
    $ ls
    $ mp-25-20180417

    
Installation of auto-kappa
============================

Because auto-kappa is currently located in a private repogitory in Github,
if you'd like to use it, please contact M. Ohnishi (ohnishi@t.u-tokyo.ac.jp).

1. Send your ssh public key to the developper (M. Ohnishi).

2. Once your key was registered, you can download auto-kappa with git command.

.. code-block:: bash
    
    $ git clone -b develop git@github.com:masato1122/auto-kappa.git


If you cannot download. Please add the following contents in ~/.ssh/config.
If you newly create the config file, you also need to modify the permission.

.. code-block:: bash
    
    Host github.com
        HostName ssh.github.com
        Port 443
        IdentityFile ~/.ssh/id_rsa  # If you changed the directory, modify this part.
        User git

    Host ssh.github.com
        Port 443
        IdentityFile ~/.ssh/id_rsa  # Same as the above
        User git
    

You can change the permission with ``chmod``.

.. code-block:: bash
    
    $ chmod 600 ~/.ssh/config


3. Create a virtual environment, ``kappa``, with conda.

.. code-block:: bash

    $ conda create -n kappa python==3.9
    $ conda init
    $ exit (You once need to logout and login to the server.)
    
    
    Login the server again and confirm the virtual environment was created.
    $ conda env list
    ...
    kappa       /home/***/***/envs/kappa
    ...
    
    
    Activate the virtual environment.
    $ conda activate kappa


To set ``kappa`` as the default, add the following line in ``.bash_profile``.

.. code-block:: bash

    source activate kappa


4. Continue to install auto-kappa.

.. code-block:: bash
    
    $ cd (arbitrary directory)/auto-kappa
    $ git config pull.rebase false
    $ git pull     ## update the package
    $ sh install.sh
     
    Check if auto_kappa is installed or not.
    $ python
    >>> import auto_kappa
    >>> exit()
    
    $ akrun -h


5. Run test examples.

.. code-block:: bash
    
    $ cd (move to an arbitrary directory outside auto-kappa directory)
    $ cp -r (auto-kappa directory)/auto-kappa/examples ./
    $ cd examples
    $ ls
    alm anphon database vasp_ase vasp_custodian
    
    $ cd alm
    $ sh run.sh
    
    $ cd ../anphon
    $ sh run.sh
    
    $ sh ../vasp_ase
    $ sh run.sh
    
    $ sh ../vasp_custodian
    $ sh run.sh    
    # This job takes time. You can stop after checking output files such as OUTCAR
    # OSZICAR, etc. were created.
    
    $ sh ../database
    $ sh run.sh
    # This job also takes time. You can stop a few minutes after starting the job.


6. ``database`` example

The calculation in ``database`` calculates thermal conductivity of Silicon automatically.
Because every process is included in this job, it takes a few hours.
It is recommended to use a job scheduler to submit this job.
An example of job script is shown below. Please modify depending on your environment.

.. code-block:: shell
    
    #!/bin/sh
    #PBS -q default         ## name of queue that you can check with a command like "qstat -q".
    #PBS -l nodes=1:ppn=24  ## only nodes=1 is available
    #PBS -j oe
    #PBS -N test            ## job name
    
    export LANG=C
    export OMP_NUM_THREADS=1  ## Please set OMP_NUM_THREADS=1
    cd $PBS_O_WORKDIR
    
    ncores=24               ## ncores must be smaller than ppn, which is set above.
    
    mpid=mp-149             ## Si
    dir_db=${directory_of_downloaded_phoonondb}/${mpid}  ## This line must be modified.
    
    if [ ! -e $dir_db ]; then
        echo " Cannot find $dir_db"    
        exit
    fi
    
    akrun \
        --directory $dir_db \
        --material_name $mpid \
        --ncores $ncores

7. Tips

* You may get warning like below. While these messages will be removed, you can neglect them which do not affect the 
  calculation. These messages are shown because POTCAR files are generated by ASE, which addes a few information in the POTCAR 
  file, and these files are read by Pymatgen, which consideres that the additional information may be error.

.. code-block:: shell

    .../lib/python3.8/site-packages/pymatgen/io/vasp/inputs.py:1738: UserWarning: Ignoring unknown variable type SHA256 
    warnings.warn(f"Ignoring unknown variable type {key}")
    .../lib/python3.8/site-packages/pymatgen/io/vasp/inputs.py:1738: UserWarning: Ignoring unknown variable type COPYR
    warnings.warn(f"Ignoring unknown variable type {key}")



.. Installation of python libraries
.. ---------------------------------
.. 
.. .. code-block:: bash
.. 
..     $ conda create -n alm python=3.8
..     $ conda activate alm
..     $ pip install pymatgen 
..     $ conda install -c conda-forge phonopy
..     $ pip install ase
..     $ pip install seekpath
..     $ pip install custodian
..     $ conda install -c conda-forge eigen
..     $ conda install -c conda-forge gcc
..     $ pip install xmltodict
..     $ pip install f90nml
..     $
..     $ conda install -c conda-forge mkl
..     $
..     $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib
.. 
.. 
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


.. Setting for POTCAR with ASE
.. -----------------------------
.. 
.. Add the following line. In the directory, potpaw_PBE exists.
.. See the following pages for details:
.. `1 (ASE) <https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html>`_ and
.. `2 (pymatgen <https://pymatgen.org/installation.html#potcar-setup>`_.
.. 
.. .. code-block:: bash
..     
..     $ cat ~/.bash_profile
..     
..     ...
..     export VASP_PP_PATH=(directory in which potpaw_PBE is located.)
..     ...
.. 
.. .. code-block:: bash
..     
..     $ cat .pmgrc.yaml
..     
..     ...
..     PMG_VASP_PSP_DIR: (directory in which potpaw_PBE is located.)
..     PMG_MAPI_KEY: **********
..     ...

.. Installation of ALM
.. ----------------------
.. 
.. .. code-block:: bash
..     
..     $ source activate alm
..     $ git clone https://github.com/ttadano/ALM.git
..     $ cd ./ALM
..     $ git pull
..     $ cd ./python
..     $ python setup.py install
.. 
.. .. For Grand-Chariot, the following line may need to be added in setup.py.
.. .. 
.. .. .. code-block:: bash
.. .. 
.. ..     os.environ["CC"] = /usr/bin/gcc
.. 
