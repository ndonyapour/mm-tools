"""Function to remove terminal residue name prefixes."""


def remove_terminal_residue_name_prefixes(
    input_pdb_path: str,
    output_pdb_path: str,
) -> None:
    """remove_terminal_residue_name_prefixes.

    Args:
        script:
        input_pdb_path: Path to the input file
        output_pdb_path: Path to the output file
    Returns:
        None
    """
    # Some AmberTools operations will add N and C prefixes to the terminal residue
    # names. This causes problems when attempting to pass the results to other
    # software, e.g. gromacs.

    # See https://proteopedia.org/wiki/index.php/Amino_Acids
    # See https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/primary-sequences-and-the-pdb-format
    l_aminos_standard = [
        "ALA",
        "ARG",
        "ASN",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
    l_aminos_titrated = ["ASP", "HID", "HIE", "HIP", "CYX"]
    l_aminos = l_aminos_standard + l_aminos_titrated + ["PYL", "SEC", "UNL"]
    nter_l_aminos = ["N" + a for a in l_aminos]
    cter_l_aminos = ["C" + a for a in l_aminos]

    # TODO: Consider D-amino acids?

    with open(input_pdb_path, encoding="utf-8") as f:  # noqa: PTH123
        lines = f.readlines()

    lines_new = []
    for line in lines:
        current_l = line
        for nter, cter, resname in zip(nter_l_aminos, cter_l_aminos, l_aminos):
            # Remove the N and C prefixes
            current_l = current_l.replace(nter, resname + " ").replace(
                cter,
                resname + " ",
            )
        lines_new.append(current_l)

    with open(output_pdb_path, mode="w", encoding="utf-8") as f:  # noqa: PTH123
        f.writelines(lines_new)
