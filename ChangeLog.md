# Ver. 1.1.3 (Aug. 19, 2026)

- **The average over degenerate modes is now taken for the scattering rate (i.e. the self-energy
  Gamma), not for the lifetime, in agreement with ANPHON.** ANPHON averages the ph-ph damping
  `damping3` before it is combined with the isotope/boundary contributions and converted to a
  lifetime (`Conductivity::compute_kappa()`), whereas `Scattering` averaged the lifetime after all
  scattering processes had been summed. Since tau is proportional to 1/Gamma, the arithmetic mean
  of the lifetimes differs from the inverse of the mean rate whenever the rates of the degenerate
  modes are not exactly identical. The averaging now happens in `set_scattering_rate_phph()`.
  (`auto_kappa/io/scattering.py`, `tests/test_analyzer_degeneracy.py`)

- **Lifetime averaging over degenerate modes no longer overwrites the averaged values.**
  `get_average_at_degenerate_point()` had two defects. (1) The band loop ran over `range(nb-1)`,
  so the highest band was never used as the starting index of a degenerate group and kept its
  zero-initialized value whenever it was not degenerate with the band below it, zeroing that
  mode's lifetime and its contribution to kappa. (2) A group was re-detected starting from every
  band it contains, and each pass overwrote the value written before, so the last mode of every
  group was left unaveraged: for `omega = [0, 0, 0, 120, 120, 300]` and `tau = [3, 5, 10, 2, 8, 1]`
  the result was `[6, 7.5, 10, 5, 8, 1]` instead of `[6, 6, 6, 5, 5, 1]`, with the error growing
  with the degeneracy (e.g. triply degenerate optical modes at the Gamma point of cubic crystals).
  Each degenerate group is now detected once and the already-averaged modes are skipped, matching
  ANPHON's `Conductivity::average_self_energy_at_degenerate_point()`.
  (`auto_kappa/io/analyzer.py`, `tests/test_analyzer_degeneracy.py`)


# Ver. 1.1.2 (Jul. 3, 2026)

- Fixed a bug in BORNINFO generation where `StructureMatcher.get_transformation()`
  returned non-identity atom mappings for highly symmetric crystals (e.g. CdAs2,
  I4₁/amd), causing Born effective charge tensors to be assigned to the wrong atoms
  and reversing the sign of off-diagonal components for symmetry-related atom pairs.
  Replaced with direct nearest-neighbour matching in fractional coordinates.


# Ver. 1.1.1 (2026)

- Arrange the workflow related to --analyze_with_larger_sc, --scph, and --four options.


# Ver. 1.1.0 (2026)

- Excessive displacements in the random displacement method were limited to the smaller of 10% of the nearest-neighbor distance and 1.0 nm.

- The bug in determining KMESH_INTERPOLATE has been fixed.

- The default value of SCPH iterative was modified to 1,000, default value of ALAMODE.


# Ver. 1.0.1 (Dec. 20, 2025)

## New

- Add --config and --optimize_klength.

- Delete --vasp_parameters.

- VASP parameters can be modified with a user configuration file read by --config option.

- K-length can be automatically optimized using --optimize_klength option.

# Ver. 1.0.0 (Nov. 22, 2025)

## Fix

- Run dfc2 command properly when the file name contains brackets

- Modify the determination method of KMES_INTERPOLATE and KMESH_SCPH

# Ver. 1.0.0 (Sept. 28, 2025)

## New

- Release the package

- Add 'command_vasp_gam' for Gamma point calculation.

- Modify the structure relaxation with EOS. (Set the maximum number of iteration and strain range.)

- Read the previously-obtained relaxed structure for the volume relaxation (strict optimization).

- Add 'calculate_forces' option. When this option is set to be "0", forces are not calculate but structure with displacement patterns are prepared.

- Add 'four', 'command_anphon_ver2', and 'frac_kdensity_4ph' options.

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

