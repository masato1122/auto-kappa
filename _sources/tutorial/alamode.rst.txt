
Preparation of ALAMODE input files
========================================

The following Python scripts can be found in ``auto-kappa/examples/6_parser/alamode_inputs``.

.. contents:: Table of contents
   :local:
   :depth: 1

Force constants
---------------

Script to suggest displacement patterns:

.. literalinclude:: ./scripts/mksuggest.py

Script to calculate force constants after force calculations:

.. literalinclude:: ./scripts/mkoptimize.py

Phonon dispersion
------------------

.. literalinclude:: ./scripts/mkband.py

Thermal conductivity
--------------------

.. literalinclude:: ./scripts/mkkappa.py
