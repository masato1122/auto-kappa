"""
Tests for auto_kappa.io.born — Born effective charge mapping.

The key regression tested here: StructureMatcher.get_transformation() can
return non-identity atom mappings (e.g. [1,0,...]) for highly symmetric
crystals such as CdAs2 (I4_1/amd) even when the two structures are identical.
The fixed implementation uses direct nearest-position matching, which always
returns the position-preserving mapping.
"""
import numpy as np
import pytest
from unittest.mock import MagicMock, patch

from auto_kappa.io.born import BORNINFO, _transform_born_charges


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_mock_borninfo(vasp_frac, alm_frac, species, lattice_matrix):
    """Return a BORNINFO instance whose VASP and ALM structures are stubbed."""
    from pymatgen.core import Structure, Lattice

    lattice = Lattice(lattice_matrix)
    vasp_struct = Structure(lattice, species, vasp_frac)
    alm_struct  = Structure(lattice, species, alm_frac)

    b = object.__new__(BORNINFO)
    b.file_vasp = "dummy_vasprun.xml"
    b.file_fcs  = "dummy_fcs.xml"
    b.vasp = MagicMock()
    b.alm  = MagicMock()
    b.vasp.structure = vasp_struct
    b.alm.structure  = alm_struct
    return b


# CdAs2 primitive cell lattice (I4_1/amd, rhombohedral setting, in Angstrom)
CDAS2_LATTICE = np.array([
    [-3.9830465,  3.9830465,  2.3371877],
    [ 3.9830465, -3.9830465,  2.3371877],
    [ 3.9830465,  3.9830465, -2.3371877],
])

CDAS2_SPECIES = ["Cd", "Cd", "As", "As", "As", "As"]

CDAS2_FRAC = np.array([
    [0.25, 0.75, 0.5 ],   # Cd1
    [0.5,  0.5,  0.0 ],   # Cd2
    [0.311, 0.061, 0.75],  # As1
    [0.061, 0.311, 0.75],  # As2
    [0.689, 0.939, 0.25],  # As3
    [0.939, 0.689, 0.25],  # As4
])


# ---------------------------------------------------------------------------
# Tests for get_transformation
# ---------------------------------------------------------------------------

class TestGetTransformation:

    def test_identity_mapping_for_identical_structures(self):
        """Identical VASP and ALM structures must produce the identity mapping."""
        b = _make_mock_borninfo(CDAS2_FRAC, CDAS2_FRAC, CDAS2_SPECIES, CDAS2_LATTICE)
        map_v2a, R = b.get_transformation()
        assert map_v2a == list(range(len(CDAS2_SPECIES))), \
            f"Expected identity mapping, got {map_v2a}"

    def test_rotation_is_identity_for_same_lattice(self):
        """When both structures share the same lattice, R must be the identity."""
        b = _make_mock_borninfo(CDAS2_FRAC, CDAS2_FRAC, CDAS2_SPECIES, CDAS2_LATTICE)
        _, R = b.get_transformation()
        np.testing.assert_allclose(R, np.eye(3), atol=1e-10)

    def test_swapped_cd_atoms_detected(self):
        """If ALM has Cd1 and Cd2 swapped, the mapping should reflect that."""
        alm_frac = CDAS2_FRAC.copy()
        alm_frac[0], alm_frac[1] = CDAS2_FRAC[1].copy(), CDAS2_FRAC[0].copy()

        b = _make_mock_borninfo(CDAS2_FRAC, alm_frac, CDAS2_SPECIES, CDAS2_LATTICE)
        map_v2a, _ = b.get_transformation()

        # VASP Cd1 (idx 0) should map to ALM idx 1 (where Cd1 ended up)
        assert map_v2a[0] == 1
        # VASP Cd2 (idx 1) should map to ALM idx 0 (where Cd2 ended up)
        assert map_v2a[1] == 0
        # As atoms unchanged
        assert map_v2a[2:] == [2, 3, 4, 5]

    def test_no_structurematcher_import(self):
        """StructureMatcher must not be imported in born.py (it was the source of the bug)."""
        import auto_kappa.io.born as born_module
        assert not hasattr(born_module, "StructureMatcher"), \
            "StructureMatcher should not be imported in born.py"


# ---------------------------------------------------------------------------
# Tests for _transform_born_charges
# ---------------------------------------------------------------------------

class TestTransformBornCharges:

    def _make_born(self, n=2, seed=0):
        rng = np.random.default_rng(seed)
        return [rng.random((3, 3)) for _ in range(n)]

    def test_identity_mapping_and_rotation_leaves_charges_unchanged(self):
        born = self._make_born(n=3)
        R = np.eye(3)
        mapping = [0, 1, 2]
        result = _transform_born_charges(born, R, mapping)
        np.testing.assert_allclose(result, np.array(born), atol=1e-12)

    def test_swap_mapping_reorders_charges(self):
        """map_v2a = [1, 0] means VASP atom 0 → ALM atom 1 and vice versa.
        After inversion, BORNINFO atom 0 should get VASP atom 1's tensor."""
        born = self._make_born(n=2)
        R = np.eye(3)
        mapping = [1, 0]   # VASP Cd1 → ALM pos 1, VASP Cd2 → ALM pos 0
        result = _transform_born_charges(born, R, mapping)
        # inv_map = argsort([1,0]) = [1,0]
        # result[0] = born[1], result[1] = born[0]
        np.testing.assert_allclose(result[0], born[1], atol=1e-12)
        np.testing.assert_allclose(result[1], born[0], atol=1e-12)
