
Preparation of ALAMODE input files
========================================

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
