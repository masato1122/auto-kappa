
# Ver. 1.0.0 (July 7, 2025)

## New

- Update the version for release the package


# Ver. 0.4.0 (April 24, 2023)

## New

- Set the primitive cell tolerance for ALAMODE calculations to match that of Spglib.

## Fix

- Change the minimum number of displacement patterns for LASSO to 10.

# Ver. 0.3 (September 8, 2023)

## New

- New options ``--analyze_with_larger_supercell``, ``--delta_max_natoms``, and ``--max_loop_for_largesc`` were added.
These options enable the calculation of harmonic properties with larger supercells when negative frequencies appear with the default size of supercell.

- A log file ``ak.log`` is created in the working directory. Reading this log file, one can easily know the status of the calculation.

# Ver. 0.2 (April 3, 2023)

## New

- A new option ``--file_structure`` was added. This option enables the automation calculation with an arbitrary structure.

## Fix

- Fix a few bugs

- k-mesh densities are properly selected.

## Modifications

- Variable name: from ``dir_phdb`` to ``base_dir`` in ``ak_script.py``

