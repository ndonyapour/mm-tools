from typing import List
from pathlib import Path
import openmm.app as omma

def extract_protein_openmm(input_pdb_path : str, output_pdb_path : str):
    '''extract_protein_openmm.

    Args:
        input_pdb_path_pattern: Input pdb file path, Type: string, File type: input, Accepted formats: pdb, Example file: https://github.com/bioexcel/biobb_structure_utils/raw/master/biobb_structure_utils/test/data/utils/cat_protein.pdb
        output_pdb_path: Output pdb file path, Type: string, File type: output, Accepted formats: pdb, Example file: https://github.com/bioexcel/biobb_structure_utils/raw/master/biobb_structure_utils/test/reference/utils/ref_cat_pdb.pdb
    Returns:
        None
    '''
    if input_pdb_path and not Path.exists(Path(input_pdb_path)):
        msg = f"Can not find file {input_pdb_path}"
        raise FileNotFoundError(msg)

    pdb = omma.PDBFile(input_pdb_path)
    bonded_atom_idxs: List[int] = []

    for bond in pdb.topology.bonds():
        bonded_atom_idxs.extend([bond.atom1.index, bond.atom2.index])

    atom_idxs = [atom.index for atom in pdb.topology.atoms()]
    stray_atom_idxs = list(set(atom_idxs) - set(bonded_atom_idxs))

    stray_atoms: List[omma.topology.Atom] = []
    for atom in pdb.topology.atoms():
        if atom.index in stray_atom_idxs:
            stray_atoms.append(atom)

    modeller = omma.Modeller(pdb.topology, pdb.positions)

     # Delete stary atoms
    modeller.delete(stray_atoms)

    # Delete water atoms
    modeller.deleteWater()

    # Delete non protein residues
    # The code has been adapted from https://github.com/openmm/pdbfixer/blob/master/pdbfixer/pdbfixer.py#L1016-L1026
    proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR', 'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
    rnaResidues = ['A', 'G', 'C', 'U', 'I']
    dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']
    keep = set(proteinResidues).union(dnaResidues).union(rnaResidues)

    # Adding N and UNK will preserve any unknown residues
    # keep.add('N')
    # keep.add('UNK')
    toDelete: List[omma.topology.Residue] = []
    for residue in modeller.topology.residues():
        if residue.name not in keep:
            toDelete.append(residue)
    modeller.delete(toDelete)


    pdb.topology = modeller.topology
    pdb.positions = modeller.positions

    with open(output_pdb_path, mode="w", encoding='utf-8') as wfile:
        omma.PDBFile.writeFile(pdb.topology, pdb.positions, wfile, keepIds=True)
