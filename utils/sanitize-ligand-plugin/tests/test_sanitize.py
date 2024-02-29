"""Test the sanitize_ligand plugin."""
import pytest
from rdkit import Chem
from src.polus.mm.utils.sanitize_ligand import attempt_fix_ligand


@pytest.mark.catch_error()
def test_kekulization_error_catch() -> None:
    """Test catching Kekulization error.

    Can't kekulize mol.  Unkekulized atoms: 6 7 8 9 10.
    """
    mol = Chem.MolFromSmiles("c1ccc(cc1)-c1nnc(n1)-c1ccccc1")
    valid_ligand, rdkit_mol = attempt_fix_ligand(mol)
    assert not valid_ligand


@pytest.mark.fix_ligand()
def test_fix_explicit_valence_error() -> None:
    """Test fixing explicit valence error.

    Explicit valence for atom # 1 C, 5, is greater than permitted
    """
    mol = Chem.MolFromSmiles("c1c(ccc2NC(CN=c(c21)(C)C)=O)O", sanitize=False)
    valid_ligand, rdkit_mol = attempt_fix_ligand(mol)
    assert valid_ligand
