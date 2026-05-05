"""
Standalone La-Proteina motif locator for PyMOL visualization (uidx_aa case only).

Given a generated PDB and the original motif PDB, computes which residues in the
generated protein correspond to the motif using the same greedy aatype + RMSD
matching logic as La-Proteina's evaluation pipeline
(``proteinfoundation/utils/motif_utils.py::pad_motif_to_full_length_unindexed``).

Public API:
    get_motif_info(generated_pdb_path, motif_pdb_path, contig_string,
                   motif_only=True, atom_selection_mode="all_atom") -> dict

Dependencies: torch, numpy, biotite (>=0.36).
No openfold / loguru / einops required.
"""

import os
from functools import lru_cache
from typing import Dict, List, Literal, Optional, Tuple

import biotite.structure.io as strucio
import numpy as np
import torch
import yaml


# ---------------------------------------------------------------------------
# Motif benchmark configuration (configs/generation/motif_dict.yaml).
#
# Per-PDB parameters (contig_string, motif_only, motif_pdb_path,
# atom_selection_mode, ...) live in the YAML — never duplicate them here.
# Tasks are keyed by uppercase names like '5YUI_AA' / '5YUI_AA_TIP' /
# '1QJG_AA_NATIVE'; the default task for a PDB id is '{PDB}_AA'.
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
MOTIF_DICT_YAML_PATH = os.path.join(
    _PROJECT_ROOT, "configs", "generation", "motif_dict.yaml"
)
# `motif_pdb_path` entries in the YAML are relative to this directory.
MOTIF_DATA_ROOT = os.path.join(_PROJECT_ROOT, "data")


@lru_cache(maxsize=1)
def _load_motif_dict() -> Dict[str, Dict]:
    with open(MOTIF_DICT_YAML_PATH) as f:
        return yaml.safe_load(f)["dataset"]["motif_dict_cfg"]


def get_task_config(task_name: str) -> Dict:
    """Return the motif_dict.yaml entry for `task_name` (e.g. '5YUI_AA')."""
    cfg = _load_motif_dict()
    key = task_name.upper()
    if key not in cfg:
        raise KeyError(
            f"Unknown task '{task_name}'. {len(cfg)} tasks defined in "
            f"{MOTIF_DICT_YAML_PATH}."
        )
    return cfg[key]


def default_task_name(pdb_id: str) -> str:
    """Map a PDB id (e.g. '5YUI') to its default '{PDB}_AA' yaml task name."""
    return pdb_id.upper() + "_AA"


def _resolve_task_name(
    task_name: Optional[str], motif_pdb_path: Optional[str]
) -> Optional[str]:
    """Pick a yaml task key from an explicit name, or infer from a PDB filename.

    Filename inference uses the file stem (e.g. '5yui_aa.pdb' → '5YUI_AA') and
    will only resolve PDBs whose default '{PDB}_AA' variant exists in the yaml.
    Task variants like '_TIP' or '_NATIVE' must be passed via `task_name`.
    """
    if task_name is not None:
        return task_name.upper()
    if motif_pdb_path is not None:
        stem = os.path.splitext(os.path.basename(motif_pdb_path))[0].upper()
        if stem in _load_motif_dict():
            return stem
    return None


def is_benchmark_pdb(pdb_id: str) -> bool:
    """True if `pdb_id`'s default '{PDB}_AA' task exists in motif_dict.yaml."""
    return default_task_name(pdb_id) in _load_motif_dict()


def benchmark_motif_pdb_path(pdb_id: str) -> str:
    """Resolve a PDB id to its benchmark motif PDB file path (yaml-backed)."""
    cfg = get_task_config(default_task_name(pdb_id))
    return os.path.join(MOTIF_DATA_ROOT, cfg["motif_pdb_path"])


# ---------------------------------------------------------------------------
# Inlined constants (mirror openfold.np.residue_constants)
# ---------------------------------------------------------------------------

atom_types: List[str] = [
    "N", "CA", "C", "CB", "O",
    "CG", "CG1", "CG2", "OG", "OG1", "SG",
    "CD", "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD",
    "CE", "CE1", "CE2", "CE3", "NE", "NE1", "NE2", "OE1", "OE2",
    "CH2", "NH1", "NH2", "OH", "CZ", "CZ2", "CZ3", "NZ", "OXT",
]
assert len(atom_types) == 37
atom_order: Dict[str, int] = {a: i for i, a in enumerate(atom_types)}

restypes: List[str] = [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
]
restype_order: Dict[str, int] = {a: i for i, a in enumerate(restypes)}
restype_num: int = len(restypes)  # 20 → "unknown" sentinel

restype_3to1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

SIDECHAIN_TIP_ATOMS: Dict[str, List[str]] = {
    "ALA": ["CA", "CB"],
    "ARG": ["CD", "CZ", "NE", "NH1", "NH2"],
    "ASP": ["CB", "CG", "OD1", "OD2"],
    "ASN": ["CB", "CG", "ND2", "OD1"],
    "CYS": ["CA", "CB", "SG"],
    "GLU": ["CG", "CD", "OE1", "OE2"],
    "GLN": ["CG", "CD", "NE2", "OE1"],
    "GLY": [],
    "HIS": ["CB", "CG", "CD2", "CE1", "ND1", "NE2"],
    "ILE": ["CB", "CG1", "CG2", "CD1"],
    "LEU": ["CB", "CG", "CD1", "CD2"],
    "LYS": ["CE", "NZ"],
    "MET": ["CG", "CE", "SD"],
    "PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["CA", "CB", "CG", "CD", "N"],
    "SER": ["CA", "CB", "OG"],
    "THR": ["CA", "CB", "CG2", "OG1"],
    "TRP": ["CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "NE1"],
    "TYR": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL": ["CB", "CG1", "CG2"],
}


# ---------------------------------------------------------------------------
# Helpers (verbatim from La-Proteina, with einops/loguru stripped)
# ---------------------------------------------------------------------------

def mean_w_mask(a: torch.Tensor, mask: torch.Tensor, keepdim: bool = True) -> torch.Tensor:
    mask = mask[..., None]
    num_elements = torch.sum(mask, dim=-2, keepdim=True)
    num_elements = torch.where(num_elements == 0, torch.tensor(1.0), num_elements)
    a_masked = torch.masked_fill(a, ~mask, 0.0)
    mean = torch.sum(a_masked, dim=-2, keepdim=True) / num_elements
    mean = torch.masked_fill(mean, num_elements == 0, 0.0)
    if not keepdim:
        mean = mean.squeeze(-2)
    return mean


def _select_motif_atoms(
    available_atoms: List[int],
    atom_selection_mode: Literal["ca", "bb3o", "all_atom", "tip_atoms"] = "ca",
    residue_name: str = None,
) -> List[int]:
    backbone_atoms = [0, 1, 2, 4]  # N, CA, C, O in atom37
    ca_index = 1
    if atom_selection_mode == "ca":
        return [ca_index] if ca_index in available_atoms else []
    if atom_selection_mode == "bb3o":
        return [i for i in backbone_atoms if i in available_atoms]
    if atom_selection_mode == "all_atom":
        return available_atoms
    if atom_selection_mode == "tip_atoms":
        if residue_name is None:
            raise ValueError("residue_name required for tip_atoms mode")
        tip_names = SIDECHAIN_TIP_ATOMS.get(residue_name, [])
        return [
            atom_order[n] for n in tip_names
            if n in atom_order and atom_order[n] in available_atoms
        ]
    raise ValueError(f"Unknown atom_selection_mode: {atom_selection_mode}")


# ---------------------------------------------------------------------------
# Motif PDB parsing (residue/range branch — sufficient for uidx_aa tasks)
# ---------------------------------------------------------------------------

def extract_motif_from_pdb(
    position: str,
    pdb_path: str,
    motif_only: bool = False,
    atom_selection_mode: Literal["ca", "bb3o", "all_atom", "tip_atoms"] = "all_atom",
    coors_to_nm: bool = False,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Parse a benchmark motif PDB into atom37 (motif_mask, x_motif, residue_type).

    `position` is the contig template from motif_dict.yaml (e.g.
    "8-15/A92-99/16-30/..."). Only chain-letter-prefixed parts are used; bare
    integers (scaffold lengths) are ignored.
    """
    ALPHABET = "ABCDEFGHJKLMNOPQRSTUVWXYZ"
    parts = position.split("/")
    array = strucio.load_structure(pdb_path, model=1)

    motif_array = []
    seen = set()
    for p in parts:
        chain_id = p[0]
        if chain_id not in ALPHABET:
            continue
        atom_mask_arr = (array.chain_id == chain_id) & (array.hetero == False)
        if motif_only:
            if chain_id in seen:
                continue
            seen.add(chain_id)
        else:
            spec = p[1:]
            if "-" in spec:
                start, end = (int(x) for x in spec.split("-"))
            else:
                start = end = int(spec)
            atom_mask_arr = atom_mask_arr & (array.res_id >= start) & (array.res_id <= end)
        motif_array.append(array[atom_mask_arr])

    motif = motif_array[0]
    for extra in motif_array[1:]:
        motif += extra

    seen_pairs = set()
    unique_residues = []
    for chain, resid in zip(motif.chain_id, motif.res_id):
        key = (chain, resid)
        if key not in seen_pairs:
            seen_pairs.add(key)
            unique_residues.append(key)

    n_res = len(unique_residues)
    motif_mask = torch.zeros((n_res, 37), dtype=torch.bool)
    x_motif = torch.zeros((n_res, 37, 3), dtype=torch.float)
    residue_type = torch.full((n_res,), restype_num, dtype=torch.int64)

    for i, (chain_id, res_id) in enumerate(unique_residues):
        res_atoms = motif[(motif.chain_id == chain_id) & (motif.res_id == res_id)]
        res1 = restype_3to1.get(res_atoms[0].res_name, "UNK")
        residue_type[i] = restype_order.get(res1, restype_num)

        available = [
            atom_order[a.atom_name] for a in res_atoms if a.atom_name in atom_order
        ]
        if not available:
            continue
        selected = set(_select_motif_atoms(available, atom_selection_mode, res_atoms[0].res_name))
        for atom in res_atoms:
            if atom.atom_name not in atom_order:
                continue
            idx = atom_order[atom.atom_name]
            if idx not in selected:
                continue
            motif_mask[i, idx] = True
            coord = torch.as_tensor(atom.coord, dtype=torch.float)
            x_motif[i, idx] = coord / 10.0 if coors_to_nm else coord

    # Centering (matches La-Proteina behaviour)
    motif_center = mean_w_mask(
        x_motif.flatten(0, 1), motif_mask.flatten(0, 1)
    ).unsqueeze(0)
    x_motif = (x_motif - motif_center) * motif_mask[..., None]
    return motif_mask, x_motif, residue_type


# ---------------------------------------------------------------------------
# Greedy aatype + RMSD matching — verbatim from La-Proteina
# ---------------------------------------------------------------------------

def pad_motif_to_full_length_unindexed(
    motif_mask: torch.Tensor,
    x_motif: torch.Tensor,
    residue_type: torch.Tensor,
    gen_coors: torch.Tensor,
    gen_mask: torch.Tensor,
    gen_aa_type: torch.Tensor,
    match_aatype: bool = True,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    nres = gen_coors.shape[0]
    nres_motif = x_motif.shape[0]
    motif_index: List[int] = []

    for i in range(nres_motif):
        motif_mask_i = motif_mask[i]
        x_motif_i = x_motif[i]
        aatype_motif_i = residue_type[i]

        best_match_idx, best_rmsd = None, float("inf")
        for j in range(nres):
            gen_mask_j = gen_mask[j]
            gen_coors_j = gen_coors[j]
            aatype_gen_j = gen_aa_type[j]

            atom_overlap = motif_mask_i & gen_mask_j
            if atom_overlap.sum() == 0:
                continue

            diff = x_motif_i[atom_overlap] - gen_coors_j[atom_overlap]
            rmsd = torch.sqrt(torch.sum(diff ** 2, dim=1).mean())

            cond = rmsd < best_rmsd and j not in motif_index
            if match_aatype:
                cond = cond and aatype_motif_i == aatype_gen_j
            if cond:
                best_rmsd = rmsd
                best_match_idx = j
        motif_index.append(best_match_idx)

    if None in motif_index:
        # Codebase fallback: assign the first n_motif residues
        motif_index = list(range(nres_motif))

    motif_mask_full = torch.zeros((nres, 37), dtype=torch.bool)
    x_motif_full = torch.zeros((nres, 37, 3), dtype=torch.float)
    residue_type_full = torch.full((nres,), restype_num, dtype=torch.int64)
    motif_mask_full[motif_index] = motif_mask
    x_motif_full[motif_index] = x_motif
    residue_type_full[motif_index] = residue_type
    return motif_mask_full, x_motif_full, residue_type_full


# ---------------------------------------------------------------------------
# Standalone biotite-based atom37 PDB loader (replaces openfold from_pdb_string)
# ---------------------------------------------------------------------------

def load_pdb_atom37(pdb_path: str) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Load a PDB into (atom_positions [L,37,3] Å, atom_mask [L,37], aatype [L])."""
    array = strucio.load_structure(pdb_path, model=1)
    array = array[array.hetero == False]

    seen = set()
    residues: List[Tuple[str, int, str]] = []
    for chain, resid, resname in zip(array.chain_id, array.res_id, array.res_name):
        key = (chain, int(resid))
        if key not in seen:
            seen.add(key)
            residues.append((chain, int(resid), resname))

    L = len(residues)
    atom_positions = torch.zeros((L, 37, 3), dtype=torch.float)
    atom_mask = torch.zeros((L, 37), dtype=torch.bool)
    aatype = torch.full((L,), restype_num, dtype=torch.int64)

    for i, (chain, resid, resname) in enumerate(residues):
        res_atoms = array[(array.chain_id == chain) & (array.res_id == resid)]
        aa1 = restype_3to1.get(resname, "UNK")
        aatype[i] = restype_order.get(aa1, restype_num)
        for atom in res_atoms:
            if atom.atom_name in atom_order:
                idx = atom_order[atom.atom_name]
                atom_positions[i, idx] = torch.as_tensor(atom.coord, dtype=torch.float)
                atom_mask[i, idx] = True

    return atom_positions, atom_mask, aatype


def _atom_selection_from_mask(motif_mask_full: torch.Tensor) -> str:
    """Build a PyMOL selection string covering only the conditioning atoms.

    `motif_mask_full` is the (L, 37) per-atom bool mask returned by the matching
    step. The selection groups conditioning atom names per residue, so PyMOL can
    target the exact atoms (e.g. tip atoms) rather than whole residues.
    """
    parts = []
    for i in range(motif_mask_full.shape[0]):
        atom_idxs = torch.where(motif_mask_full[i])[0].tolist()
        if not atom_idxs:
            continue
        names = "+".join(atom_types[j] for j in atom_idxs)
        parts.append(f"(resi {i + 1} and name {names})")
    return " or ".join(parts)


# ---------------------------------------------------------------------------
# Top-level convenience function
# ---------------------------------------------------------------------------

def get_motif_info(
    generated_pdb_path: str,
    motif_pdb_path: Optional[str] = None,
    task_name: Optional[str] = None,
    contig_string: Optional[str] = None,
    motif_only: Optional[bool] = None,
    atom_selection_mode: Optional[Literal["ca", "bb3o", "all_atom", "tip_atoms"]] = None,
    match_aatype: bool = True,
) -> Dict:
    """Identify motif residues in a La-Proteina-generated PDB (uidx_aa case).

    All per-PDB params (contig_string, motif_only, motif_pdb_path,
    atom_selection_mode) default to the matching motif_dict.yaml entry. The
    entry is resolved from `task_name` if given, otherwise inferred from
    `motif_pdb_path`'s filename stem (e.g. '5yui_aa.pdb' → '5YUI_AA').
    Explicit kwargs always win over yaml values.

    Args:
        generated_pdb_path:  Path to the generated protein PDB.
        motif_pdb_path:      Path to the original motif PDB. Optional if
                             `task_name` is given (yaml supplies the path).
        task_name:           Yaml task key (e.g. '5YUI_AA', '5YUI_AA_TIP',
                             '1QJG_AA_NATIVE'). Required to disambiguate task
                             variants like _TIP / _NATIVE.
        contig_string:       Override yaml's contig_string.
        motif_only:          Override yaml's motif_only flag.
        atom_selection_mode: Override yaml's atom_selection_mode.
        match_aatype:        La-Proteina default is True.

    Returns:
        {
          "motif_resi":       [int, ...]     # 1-based residue numbers (PyMOL-friendly)
          "selection":        "12+13+14+..." # residue-level PyMOL selection
          "atom_selection":   "(resi 12 and name CG+ND1+...) or ..."  # atom-level
          "motif_mask_full":  Tensor (L, 37) # per-atom bool mask
          "n_motif":          int
          "n_total":          int
          "motif_aa":         str            # 1-letter motif sequence in matched order
        }
    """
    resolved = _resolve_task_name(task_name, motif_pdb_path)
    cfg = get_task_config(resolved) if resolved else None

    if cfg is not None:
        if contig_string is None:
            contig_string = cfg["contig_string"]
        if motif_only is None:
            motif_only = cfg["motif_only"]
        if atom_selection_mode is None:
            atom_selection_mode = cfg["atom_selection_mode"]
        if motif_pdb_path is None:
            motif_pdb_path = os.path.join(MOTIF_DATA_ROOT, cfg["motif_pdb_path"])

    if atom_selection_mode is None:
        atom_selection_mode = "all_atom"
    if motif_pdb_path is None or contig_string is None or motif_only is None:
        raise ValueError(
            "get_motif_info needs a yaml task or explicit "
            "motif_pdb_path/contig_string/motif_only; none could be resolved."
        )

    motif_mask, x_motif, residue_type = extract_motif_from_pdb(
        position=contig_string,
        pdb_path=motif_pdb_path,
        motif_only=motif_only,
        atom_selection_mode=atom_selection_mode,
        coors_to_nm=False,
    )
    gen_coors, gen_mask, gen_aa_type = load_pdb_atom37(generated_pdb_path)

    motif_mask_full, _, residue_type_full = pad_motif_to_full_length_unindexed(
        motif_mask=motif_mask,
        x_motif=x_motif,
        residue_type=residue_type,
        gen_coors=gen_coors,
        gen_mask=gen_mask,
        gen_aa_type=gen_aa_type,
        match_aatype=match_aatype,
    )

    motif_residue_mask = motif_mask_full.any(dim=-1)
    motif_resi = (motif_residue_mask.nonzero().flatten() + 1).tolist()

    inv_restype = {v: k for k, v in restype_order.items()}
    motif_aa = "".join(
        inv_restype.get(int(t), "X") for t in residue_type_full[motif_residue_mask]
    )

    return {
        "motif_resi": motif_resi,
        "selection": "+".join(str(r) for r in motif_resi),
        "atom_selection": _atom_selection_from_mask(motif_mask_full),
        "motif_mask_full": motif_mask_full,
        "n_motif": len(motif_resi),
        "n_total": int(motif_mask_full.shape[0]),
        "motif_aa": motif_aa,
    }
