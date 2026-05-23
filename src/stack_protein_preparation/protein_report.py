"""
Per-protein PDF report for the FRUTON protein preparation pipeline.

Generates one PDF per protein containing:
  - A PyMOL-rendered figure of the starting structure (centered on ligand
    when one is present)
  - A prose narrative for every pipeline step that ran, with inline citations
  - A bibliography section with full references
  - A shared BibTeX file written to the proteins directory

Entry point: generate_protein_report()
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

# ---------------------------------------------------------------------------
# Visual identity constants
# ---------------------------------------------------------------------------

# Palette sampled from the uploaded FRUTON and IDIS logos and cross-checked
# against the IDIS website SVG/icon accent colour (#8A730F).  The values are
# deliberately centralized so future report, workbook, and dossier outputs can
# share one visual language without changing the scientific reporting logic.
_IDIS_NAVY = "#28325A"
_IDIS_NAVY_2 = "#323264"
_IDIS_TEXT = "#252A3F"
_IDIS_MUTED = "#697086"
_IDIS_LINE = "#D8DCE7"
_IDIS_PANEL = "#F5F7FB"
_IDIS_PANEL_2 = "#EEF1F7"
_IDIS_GOLD = "#8A730F"
_FRUTON_RED = "#E60014"
_FRUTON_RED_SOFT = "#FDECEE"
_FRUTON_YELLOW = "#FAE646"
_FRUTON_YELLOW_SOFT = "#FFF8CF"
_FRUTON_GREEN = "#50961E"
_FRUTON_GREEN_SOFT = "#EAF4E3"
_AUDIT_GREY = "#F2F3F6"
_AUDIT_WARN = "#FFF3C4"
_AUDIT_FAIL = "#F9D7DB"
_AUDIT_OK = "#E7F3E1"

_FRUTON_EXPANSION = (
    "Framework for Reconstruction, UniProt alignment, and "
    "Topology-Oriented protein Normalization"
)

_IDIS_LOGO_NAMES = (
    "logo_idis.png",
    "logo_IDIS_2020-1.png",
    "logo_IDIS_2020-1(1).png",
    "idis_logo.png",
)
_FRUTON_LOGO_NAMES = (
    "logo.png",
    "logo(1).png",
    "fruton_logo.png",
    "FRUTON_logo.png",
)

# ---------------------------------------------------------------------------
# BibTeX bibliography (16 entries)
# ---------------------------------------------------------------------------

_BIBTEX = r"""@article{berman2000pdb,
  author  = {Berman, Helen M. and Westbrook, John and Feng, Zukang and
             Gilliland, Gary and Bhat, T. N. and Weissig, Helge and
             Shindyalov, Ilya N. and Bourne, Philip E.},
  title   = {The Protein Data Bank},
  journal = {Nucleic Acids Research},
  year    = {2000},
  volume  = {28},
  number  = {1},
  pages   = {235--242},
  doi     = {10.1093/nar/28.1.235},
}

@article{uniprot2023,
  author  = {{UniProt Consortium}},
  title   = {{UniProt}: the Universal Protein Database in 2023},
  journal = {Nucleic Acids Research},
  year    = {2023},
  volume  = {51},
  number  = {D1},
  pages   = {D523--D531},
  doi     = {10.1093/nar/gkac1052},
}

@article{cock2009biopython,
  author  = {Cock, Peter J. A. and Antao, Tiago and Chang, Jeffrey T. and
             Chapman, Brad A. and Cox, Cymon J. and Dalke, Andrew and
             Friedberg, Iddo and Hamelryck, Thomas and Kauff, Frank and
             Wilczynski, Bartek and de Hoon, Michiel J. L.},
  title   = {Biopython: freely available {Python} tools for computational
             molecular biology and bioinformatics},
  journal = {Bioinformatics},
  year    = {2009},
  volume  = {25},
  number  = {11},
  pages   = {1422--1423},
  doi     = {10.1093/bioinformatics/btp163},
}

@article{sali1993modeller,
  author  = {{\v{S}}ali, Andrej and Blundell, Tom L.},
  title   = {Comparative Protein Modelling by Satisfaction of Spatial
             Restraints},
  journal = {Journal of Molecular Biology},
  year    = {1993},
  volume  = {234},
  number  = {3},
  pages   = {779--815},
  doi     = {10.1006/jmbi.1993.1626},
}

@incollection{webb2016modeller,
  author    = {Webb, Benjamin and Sali, Andrej},
  title     = {Comparative Protein Structure Modeling Using {MODELLER}},
  booktitle = {Current Protocols in Bioinformatics},
  year      = {2016},
  volume    = {54},
  pages     = {5.6.1--5.6.37},
  doi       = {10.1002/cpbi.3},
}

@article{jumper2021alphafold,
  author  = {Jumper, John and Evans, Richard and Pritzel, Alexander and
             Green, Tim and Figurnov, Michael and Ronneberger, Olaf and
             Tunyasuvunakool, Kathryn and Bates, Russ and {\v{Z}}{\'\i}dek,
             Augustin and Potapenko, Anna and Bridgland, Alex and Meyer,
             Clemens and Kohl, Simon A. A. and Ballard, Andrew J. and
             Cowie, Andrew and Romera-Paredes, Bernardino and
             Nikolov, Stanimir and Jain, Rishub and Adler, Jonas and
             Back, Trevor and Petersen, Stig and Reiman, David and
             Clancy, Ellen and Zielinski, Michal and Steinegger, Martin and
             Pacholska, Michalina and Berghammer, Tamas and Bodenstein,
             Sebastian and Silver, David and Vinyals, Oriol and
             Senior, Andrew W. and Kavukcuoglu, Koray and Kohli, Pushmeet
             and Hassabis, Demis},
  title   = {Highly accurate protein structure prediction with {AlphaFold}},
  journal = {Nature},
  year    = {2021},
  volume  = {596},
  pages   = {583--589},
  doi     = {10.1038/s41586-021-03819-2},
}

@article{abraham2015gromacs,
  author  = {Abraham, Mark James and Murtola, Teemu and Schulz, Roland and
             P{\'a}ll, Szil{\'a}rd and Smith, Jeremy C. and Hess, Berk and
             Lindahl, Erik},
  title   = {{GROMACS}: High performance molecular simulations through
             multi-level parallelism from laptops to supercomputers},
  journal = {SoftwareX},
  year    = {2015},
  volume  = {1--2},
  pages   = {19--25},
  doi     = {10.1016/j.softx.2015.06.001},
}

@article{lindorff2010ff99sb,
  author  = {Lindorff-Larsen, Kresten and Piana, Stefano and Palmo, Kim and
             Maragakis, Paul and Klepeis, John L. and Dror, Ron O. and
             Shaw, David E.},
  title   = {Improved side-chain torsion potentials for the {Amber}
             ff99{SB} protein force field},
  journal = {Proteins: Structure, Function, and Bioinformatics},
  year    = {2010},
  volume  = {78},
  number  = {8},
  pages   = {1950--1958},
  doi     = {10.1002/prot.22711},
}

@article{jorgensen1983tip3p,
  author  = {Jorgensen, William L. and Chandrasekhar, Jayaraman and
             Madura, Jeffry D. and Impey, Roger W. and Klein, Michael L.},
  title   = {Comparison of simple potential functions for simulating liquid
             water},
  journal = {The Journal of Chemical Physics},
  year    = {1983},
  volume  = {79},
  number  = {2},
  pages   = {926--935},
  doi     = {10.1063/1.445869},
}

@article{li2016mcpb,
  author  = {Li, Pengfei and Merz, Kenneth M.},
  title   = {{MCPB.py}: A Python Based Metal Center Parameter Builder},
  journal = {Journal of Chemical Information and Modeling},
  year    = {2016},
  volume  = {56},
  number  = {4},
  pages   = {599--604},
  doi     = {10.1021/acs.jcim.5b00674},
}

@misc{frisch2016gaussian,
  author = {Frisch, M. J. and Trucks, G. W. and Schlegel, H. B. and
            Scuseria, G. E. and Robb, M. A. and Cheeseman, J. R. and
            Scalmani, G. and Barone, V. and Petersson, G. A. and
            Nakatsuji, H. and Li, X. and Caricato, M. and Marenich, A. V.
            and Bloino, J. and Janesko, B. G. and Gomperts, R. and
            Mennucci, B. and Hratchian, H. P. and Ortiz, J. V. and
            Izmaylov, A. F. and Sonnenberg, J. L. and Williams-Young, D.
            and Ding, F. and Lipparini, F. and Egidi, F. and Goings, J.
            and Peng, B. and Petrone, A. and Henderson, T. and Ranasinghe,
            D. and Zakrzewski, V. G. and Gao, J. and Rega, N. and Zheng,
            G. and Liang, W. and Hada, M. and Ehara, M. and Toyota, K.
            and Fukuda, R. and Hasegawa, J. and Ishida, M. and Nakajima,
            T. and Honda, Y. and Kitao, O. and Nakai, H. and Vreven, T.
            and Throssell, K. and Montgomery, Jr., J. A. and Peralta, J. E.
            and Ogliaro, F. and Bearpark, M. J. and Heyd, J. J. and
            Brothers, E. N. and Kudin, K. N. and Staroverov, V. N. and
            Keith, T. A. and Kobayashi, R. and Normand, J. and Raghavachari,
            K. and Rendell, A. P. and Burant, J. C. and Iyengar, S. S. and
            Tomasi, J. and Cossi, M. and Millam, J. M. and Klene, M. and
            Adamo, C. and Cammi, R. and Ochterski, J. W. and Martin, R. L.
            and Morokuma, K. and Farkas, O. and Foresman, J. B. and
            Fox, D. J.},
  title  = {{Gaussian~16, Revision C.01}},
  year   = {2016},
  note   = {Gaussian, Inc., Wallingford CT},
}

@article{bayly1993resp,
  author  = {Bayly, Christopher I. and Cieplak, Piotr and Cornell, Wendy and
             Kollman, Peter A.},
  title   = {A well-behaved electrostatic potential based method using charge
             restraints for deriving atomic charges: the {RESP} model},
  journal = {The Journal of Physical Chemistry},
  year    = {1993},
  volume  = {97},
  number  = {40},
  pages   = {10269--10280},
  doi     = {10.1021/j100142a004},
}

@article{lovell2003rama,
  author  = {Lovell, Simon C. and Davis, Ian W. and Arendall, III, W. Bryan
             and de Bakker, Paul I. W. and Word, J. Michael and Prisant,
             Michael G. and Richardson, Jane S. and Richardson, David C.},
  title   = {Structure validation by {C}alpha geometry: phi, psi and
             {C}beta deviation},
  journal = {Proteins: Structure, Function, and Genetics},
  year    = {2003},
  volume  = {50},
  number  = {3},
  pages   = {437--450},
  doi     = {10.1002/prot.10286},
}

@article{williams2018molprobity,
  author  = {Williams, Christopher J. and Headd, Jeffrey J. and
             Moriarty, Nigel W. and Prisant, Michael G. and Videau, Lizbeth
             L. and Deis, Lindsay N. and Verma, Vishal and Keedy, Daniel A.
             and Hintze, Bradley J. and Chen, Vincent B. and Jain, Swati and
             Lewis, Stephen M. and Arendall, III, W. Bryan and Snoeyink,
             Jack and Adams, Paul D. and Lovell, Simon C. and Richardson,
             Jane S. and Richardson, David C.},
  title   = {{MolProbity}: More and better reference data for improved
             all-atom structure validation},
  journal = {Protein Science},
  year    = {2018},
  volume  = {27},
  number  = {1},
  pages   = {293--315},
  doi     = {10.1002/pro.3330},
}

@article{laskowski1993procheck,
  author  = {Laskowski, Roman A. and MacArthur, Malcolm W. and Moss,
             David S. and Thornton, Janet M.},
  title   = {{PROCHECK}: a program to check the stereochemical quality of
             protein structures},
  journal = {Journal of Applied Crystallography},
  year    = {1993},
  volume  = {26},
  number  = {2},
  pages   = {283--291},
  doi     = {10.1107/S0021889892009944},
}

@misc{schrodinger2015pymol,
  author = {{Schr{\"o}dinger, LLC}},
  title  = {The {PyMOL} Molecular Graphics System, Version~2.0},
  year   = {2015},
  note   = {\url{https://www.pymol.org}},
}

@article{katoh2013mafft,
  author  = {Katoh, Kazutaka and Standley, Daron M.},
  title   = {{MAFFT} Multiple Sequence Alignment Software Version 7:
             Improvements in Performance and Usability},
  journal = {Molecular Biology and Evolution},
  year    = {2013},
  volume  = {30},
  number  = {4},
  pages   = {772--780},
  doi     = {10.1093/molbev/mst010},
}

@article{olsson2011propka,
  author  = {Olsson, Mats H. M. and S{\o}ndergaard, Chresten R. and
             Rostkowski, Michal and Jensen, Jan H.},
  title   = {{PROPKA3}: Consistent Treatment of Internal and Surface
             Residues in Empirical p{K}a Predictions},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2011},
  volume  = {7},
  number  = {2},
  pages   = {525--537},
  doi     = {10.1021/ct100578z},
}

@article{sondergaard2011propka,
  author  = {S{\o}ndergaard, Chresten R. and Olsson, Mats H. M. and
             Rostkowski, Michal and Jensen, Jan H.},
  title   = {Improved Treatment of Ligands and Coupling Effects in
             Empirical Calculation and Rationalization of p{K}$_a$ Values},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2011},
  volume  = {7},
  number  = {7},
  pages   = {2284--2295},
  doi     = {10.1021/ct200133y},
}

@article{hornak2006ff99sb,
  author  = {Hornak, Viktor and Abel, Robert and Okur, Asim and
             Strockbine, Bentley and Roitberg, Adrian and Simmerling, Carlos},
  title   = {Comparison of Multiple {Amber} Force Fields and Development of
             Improved Protein Backbone Parameters},
  journal = {Proteins: Structure, Function, and Bioinformatics},
  year    = {2006},
  volume  = {65},
  number  = {3},
  pages   = {712--725},
  doi     = {10.1002/prot.21123},
}

@article{wang2004gaff,
  author  = {Wang, Junmei and Wolf, Romain M. and Caldwell, James W. and
             Kollman, Peter A. and Case, David A.},
  title   = {Development and Testing of a General Amber Force Field},
  journal = {Journal of Computational Chemistry},
  year    = {2004},
  volume  = {25},
  number  = {9},
  pages   = {1157--1174},
  doi     = {10.1002/jcc.20035},
}

@article{becke1993b3lyp,
  author  = {Becke, Axel D.},
  title   = {Density-functional thermochemistry. {III}. The role of exact
             exchange},
  journal = {The Journal of Chemical Physics},
  year    = {1993},
  volume  = {98},
  number  = {7},
  pages   = {5648--5652},
  doi     = {10.1063/1.464913},
}

@article{lee1988lyp,
  author  = {Lee, Chengteh and Yang, Weitao and Parr, Robert G.},
  title   = {Development of the {Colle-Salvetti} correlation-energy formula
             into a functional of the electron density},
  journal = {Physical Review B},
  year    = {1988},
  volume  = {37},
  number  = {2},
  pages   = {785--789},
  doi     = {10.1103/PhysRevB.37.785},
}

@article{hehre1972basis,
  author  = {Hehre, Warren J. and Ditchfield, Robert and Pople, John A.},
  title   = {Self-Consistent Molecular Orbital Methods. {XII}. Further
             Extensions of {Gaussian}-Type Basis Sets for Use in Molecular
             Orbital Studies of Organic Molecules},
  journal = {The Journal of Chemical Physics},
  year    = {1972},
  volume  = {56},
  number  = {5},
  pages   = {2257--2261},
  doi     = {10.1063/1.1677527},
}

@article{francl1982basis,
  author  = {Francl, Michelle M. and Pietro, William J. and Hehre, Warren J. and
             Binkley, J. Stephen and Gordon, Mark S. and DeFrees, Douglas J. and
             Pople, John A.},
  title   = {Self-consistent molecular orbital methods. {XXIII}. A
             polarization-type basis set for second-row elements},
  journal = {The Journal of Chemical Physics},
  year    = {1982},
  volume  = {77},
  number  = {7},
  pages   = {3654--3665},
  doi     = {10.1063/1.444267},
}
"""

# ---------------------------------------------------------------------------
# Inline reference texts (one entry per key)
# ---------------------------------------------------------------------------

_INLINE_REFS: dict[str, str] = {
    "berman2000pdb":
        "Berman HM et al. (2000) The Protein Data Bank. "
        "Nucleic Acids Res 28:235-242. doi:10.1093/nar/28.1.235",
    "uniprot2023":
        "UniProt Consortium (2023) UniProt: the Universal Protein Database in 2023. "
        "Nucleic Acids Res 51:D523-D531. doi:10.1093/nar/gkac1052",
    "cock2009biopython":
        "Cock PJA et al. (2009) Biopython: freely available Python tools for "
        "computational molecular biology and bioinformatics. "
        "Bioinformatics 25:1422-1423. doi:10.1093/bioinformatics/btp163",
    "sali1993modeller":
        "Sali A, Blundell TL (1993) Comparative Protein Modelling by Satisfaction "
        "of Spatial Restraints. J Mol Biol 234:779-815. "
        "doi:10.1006/jmbi.1993.1626",
    "webb2016modeller":
        "Webb B, Sali A (2016) Comparative Protein Structure Modeling Using MODELLER. "
        "Curr Protoc Bioinformatics 54:5.6.1-5.6.37. doi:10.1002/cpbi.3",
    "jumper2021alphafold":
        "Jumper J et al. (2021) Highly accurate protein structure prediction with "
        "AlphaFold. Nature 596:583-589. doi:10.1038/s41586-021-03819-2",
    "abraham2015gromacs":
        "Abraham MJ et al. (2015) GROMACS: High performance molecular simulations "
        "through multi-level parallelism from laptops to supercomputers. "
        "SoftwareX 1-2:19-25. doi:10.1016/j.softx.2015.06.001",
    "lindorff2010ff99sb":
        "Lindorff-Larsen K et al. (2010) Improved side-chain torsion potentials for "
        "the Amber ff99SB protein force field. "
        "Proteins 78:1950-1958. doi:10.1002/prot.22711",
    "jorgensen1983tip3p":
        "Jorgensen WL et al. (1983) Comparison of simple potential functions for "
        "simulating liquid water. J Chem Phys 79:926-935. doi:10.1063/1.445869",
    "li2016mcpb":
        "Li P, Merz KM (2016) MCPB.py: A Python Based Metal Center Parameter Builder. "
        "J Chem Inf Model 56:599-604. doi:10.1021/acs.jcim.5b00674",
    "frisch2016gaussian":
        "Frisch MJ et al. (2016) Gaussian 16, Revision C.01. "
        "Gaussian, Inc., Wallingford CT.",
    "bayly1993resp":
        "Bayly CI et al. (1993) A well-behaved electrostatic potential based method "
        "using charge restraints for deriving atomic charges: the RESP model. "
        "J Phys Chem 97:10269-10280. doi:10.1021/j100142a004",
    "lovell2003rama":
        "Lovell SC et al. (2003) Structure validation by Calpha geometry: phi, psi "
        "and Cbeta deviation. Proteins 50:437-450. doi:10.1002/prot.10286",
    "williams2018molprobity":
        "Williams CJ et al. (2018) MolProbity: More and better reference data for "
        "improved all-atom structure validation. "
        "Protein Sci 27:293-315. doi:10.1002/pro.3330",
    "laskowski1993procheck":
        "Laskowski RA et al. (1993) PROCHECK: a program to check the stereochemical "
        "quality of protein structures. "
        "J Appl Crystallogr 26:283-291. doi:10.1107/S0021889892009944",
    "schrodinger2015pymol":
        "Schrodinger LLC (2015) The PyMOL Molecular Graphics System, Version 2.0. "
        "https://www.pymol.org",
    "katoh2013mafft":
        "Katoh K, Standley DM (2013) MAFFT Multiple Sequence Alignment Software "
        "Version 7: Improvements in Performance and Usability. "
        "Mol Biol Evol 30:772-780. doi:10.1093/molbev/mst010",
    "olsson2011propka":
        "Olsson MHM et al. (2011) PROPKA3: Consistent Treatment of Internal and "
        "Surface Residues in Empirical pKa Predictions. "
        "J Chem Theory Comput 7:525-537. doi:10.1021/ct100578z",
    "sondergaard2011propka":
        "Sondergaard CR et al. (2011) Improved Treatment of Ligands and Coupling "
        "Effects in Empirical Calculation and Rationalization of pKa Values. "
        "J Chem Theory Comput 7:2284-2295. doi:10.1021/ct200133y",
    "hornak2006ff99sb":
        "Hornak V et al. (2006) Comparison of Multiple Amber Force Fields and "
        "Development of Improved Protein Backbone Parameters. "
        "Proteins 65:712-725. doi:10.1002/prot.21123",
    "wang2004gaff":
        "Wang J et al. (2004) Development and Testing of a General Amber Force "
        "Field. J Comput Chem 25:1157-1174. doi:10.1002/jcc.20035",
    "becke1993b3lyp":
        "Becke AD (1993) Density-functional thermochemistry. III. The role of "
        "exact exchange. J Chem Phys 98:5648-5652. doi:10.1063/1.464913",
    "lee1988lyp":
        "Lee C, Yang W, Parr RG (1988) Development of the Colle-Salvetti "
        "correlation-energy formula into a functional of the electron density. "
        "Phys Rev B 37:785-789. doi:10.1103/PhysRevB.37.785",
    "hehre1972basis":
        "Hehre WJ, Ditchfield R, Pople JA (1972) Self-Consistent Molecular "
        "Orbital Methods. XII. Further Extensions of Gaussian-Type Basis Sets. "
        "J Chem Phys 56:2257-2261. doi:10.1063/1.1677527",
    "francl1982basis":
        "Francl MM et al. (1982) Self-consistent molecular orbital methods. XXIII. "
        "A polarization-type basis set for second-row elements. "
        "J Chem Phys 77:3654-3665. doi:10.1063/1.444267",
}


# ---------------------------------------------------------------------------
# PDB crystal paper lookup via RCSB REST API
# ---------------------------------------------------------------------------

def _fetch_pdb_citation(pdb_id: str) -> dict[str, str] | None:
    """
    Query the RCSB REST API for the primary citation of *pdb_id*.

    Returns a dict with keys: authors, title, journal, year, volume, pages, doi.
    Returns None on network error or missing data.
    """
    import urllib.request
    import json as _json

    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = _json.loads(resp.read().decode())
    except Exception:
        return None

    cit = data.get("rcsb_primary_citation")
    if not cit:
        return None

    authors_raw = cit.get("rcsb_authors") or []
    if isinstance(authors_raw, list):
        author_str = ", ".join(str(a) for a in authors_raw[:6])
        if len(authors_raw) > 6:
            author_str += " et al."
    else:
        author_str = str(authors_raw)

    return {
        "authors": author_str,
        "title":   str(cit.get("title", "") or ""),
        "journal": str(cit.get("journal_abbrev", "") or cit.get("journal_full", "") or ""),
        "year":    str(cit.get("year", "") or ""),
        "volume":  str(cit.get("journal_volume", "") or ""),
        "pages":   str(cit.get("page_first", "") or ""),
        "doi":     str(cit.get("pdbx_database_id_doi", "") or ""),
    }


def _format_pdb_citation_inline(cit: dict[str, str], pdb_id: str) -> str:
    """Format a PDB citation dict into a compact inline string for prose."""
    parts = []
    if cit.get("authors"):
        parts.append(cit["authors"])
    if cit.get("year"):
        parts.append(f"({cit['year']})")
    title = cit.get("title", "").rstrip(".")
    if title:
        parts.append(title + ".")
    journal_parts = []
    if cit.get("journal"):
        journal_parts.append(cit["journal"])
    vol = cit.get("volume", "")
    pg = cit.get("pages", "")
    if vol and pg:
        journal_parts.append(f" {vol}:{pg}")
    elif vol:
        journal_parts.append(f" {vol}")
    elif pg:
        journal_parts.append(f" :{pg}")
    if journal_parts:
        parts.append("".join(journal_parts) + ".")
    if cit.get("doi"):
        parts.append(f"doi:{cit['doi']}")
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Citation key collector — determines which refs to include in the bibliography
# ---------------------------------------------------------------------------

def _collect_ref_keys(pipeline_record: dict[str, str]) -> list[str]:
    """Return an ordered, deduplicated list of citation keys for steps that ran."""
    keys: list[str] = []
    seen: set[str] = set()

    def add(*ks: str) -> None:
        for k in ks:
            if k not in seen:
                seen.add(k)
                keys.append(k)

    def s(col: str) -> str:
        return str(pipeline_record.get(col, "")).strip()

    add("schrodinger2015pymol")  # always included: figure section

    if s("pdb_sync_done") in ("success", "warning"):
        add("berman2000pdb")

    if s("fasta_files_done") in ("success", "warning"):
        add("uniprot2023", "cock2009biopython")

    if s("sequence_alignment_done") in ("success", "warning"):
        add("katoh2013mafft", "cock2009biopython")

    if s("filler.status") in ("success", "warning"):
        src = s("filler.model_source").lower()
        if "modeller" in src:
            add("sali1993modeller", "webb2016modeller")
        elif "alphafold" in src:
            add("jumper2021alphafold")

    if s("protonation.status") in ("success", "warning"):
        add("olsson2011propka", "sondergaard2011propka",
            "hornak2006ff99sb", "lindorff2010ff99sb",
            "abraham2015gromacs", "jorgensen1983tip3p")

    if s("model_eval.status") in ("success", "warning"):
        add("lovell2003rama", "williams2018molprobity", "laskowski1993procheck")

    if s("metall_params.status") in ("success", "warning"):
        add("li2016mcpb", "becke1993b3lyp", "lee1988lyp",
            "hehre1972basis", "francl1982basis", "frisch2016gaussian")

    if s("nonstd_residue_params.status") in ("success", "warning"):
        add("wang2004gaff", "bayly1993resp",
            "hehre1972basis", "francl1982basis", "frisch2016gaussian")

    return keys


# ---------------------------------------------------------------------------
# Step narrative builder
# ---------------------------------------------------------------------------

def _cite(*keys: str) -> str:
    return f" [{'; '.join(keys)}]" if keys else ""


def _build_narratives(
    pdb_id: str,
    record: dict[str, str],
    pdb_citation_inline: str = "",
) -> list[tuple[int, str]]:
    """
    Return (fruton_step_number, prose_text) tuples for every completed step.

    The step number matches the fruton pipeline step (e.g. step 5 = fasta_files)
    and is used by the PDF builder to insert the 'input modified by user' marker
    before the first paragraph whose step number >= rerun.from_step.
    """
    paras: list[tuple[int, str]] = []

    def s(col: str) -> str:
        return str(record.get(col, "")).strip()

    def add(step: int, text: str) -> None:
        paras.append((step, text))

    done = {"success", "warning"}

    if s("pdb_sync_done") in done:
        crystal_ref = (
            f" The deposited structure was originally reported in: "
            f"{pdb_citation_inline}."
            if pdb_citation_inline else ""
        )
        add(1,
            f"The crystal structure of {pdb_id} was retrieved from the "
            f"Protein Data Bank{_cite('berman2000pdb')} and stored locally "
            f"for subsequent processing.{crystal_ref}"
        )

    if s("fasta_files_done") in done:
        uniprot = s("uniprot_id")
        uid_str = f" (UniProt accession: {uniprot})" if uniprot else ""
        add(5,
            f"The canonical amino acid sequence{uid_str} was retrieved from "
            f"UniProt{_cite('uniprot2023')}. Sequence input/output used "
            f"Biopython{_cite('cock2009biopython')}."
        )

    if s("sequence_alignment_done") in done:
        add(6,
            f"Pairwise sequence alignment between the crystallographic chain "
            f"and the UniProt reference sequence was performed with "
            f"MAFFT{_cite('katoh2013mafft')} via Biopython{_cite('cock2009biopython')} "
            f"to identify insertions, deletions, and missing terminal residues."
        )

    if s("insertion_codes_done") in done:
        add(7,
            f"PDB insertion codes were resolved and residue numbering was "
            f"re-indexed to produce a monotonically numbered structure "
            f"suitable for GROMACS topology builders."
        )

    n_gaps = s("n_gaps")
    gap_sizes = s("gap_sizes")
    has_gaps = s("has_gaps").lower()
    if n_gaps and n_gaps != "0":
        gap_str = f" with sizes {gap_sizes}" if gap_sizes else ""
        add(9,
            f"Gap analysis revealed {n_gaps} missing-residue region(s){gap_str} "
            f"in the deposited model relative to the reference sequence."
        )
    elif has_gaps in ("false", "no"):
        add(9,
            "Gap analysis confirmed that the deposited model contains no "
            "missing internal residues; gap-filling was not required."
        )

    filler_status = s("filler.status")
    filler_source = s("filler.model_source").lower()
    if filler_status in done:
        if "modeller" in filler_source:
            add(10,
                f"Missing residues were reconstructed by comparative modelling "
                f"with MODELLER{_cite('sali1993modeller', 'webb2016modeller')}, "
                f"using the best-scoring model from a multi-model run. "
                f"The filled structure was retained for downstream protonation."
            )
        elif "alphafold" in filler_source:
            add(10,
                f"Missing regions were grafted from an AlphaFold2 "
                f"predicted structure{_cite('jumper2021alphafold')}, "
                f"superimposed onto the crystallographic frame by least-squares "
                f"fitting of the shared backbone atoms."
            )
        elif filler_source:
            add(10, f"Gap-filling was performed using source: {filler_source}.")

    if s("protonation.status") in done:
        src = s("protonation.input_source")
        src_str = f" using {src} as input" if src else ""
        ph_raw = s("protonation.ph") or "7.4"
        try:
            ph_str = f"{float(ph_raw):.1f}"
        except (ValueError, TypeError):
            ph_str = str(ph_raw)
        add(11,
            f"Histidine pKa values were predicted at pH {ph_str} using "
            f"PROPKA 3{_cite('olsson2011propka', 'sondergaard2011propka')}. "
            f"Histidines with predicted pKa above pH {ph_str} were assigned "
            f"as HIP (doubly protonated, net charge +1); all others were "
            f"assigned as HIE (epsilon-tautomer, neutral). "
            f"Hydrogen atoms were then placed by GROMACS pdb2gmx "
            f"(-ignh){_cite('abraham2015gromacs')}{src_str}. "
            f"Partial charges and bonded parameters were assigned with the "
            f"AMBER ff99SB force field{_cite('hornak2006ff99sb')} as extended "
            f"by the ff99SB-ILDN side-chain torsion corrections"
            f"{_cite('lindorff2010ff99sb')} and TIP3P water "
            f"geometry{_cite('jorgensen1983tip3p')}."
        )

    if s("internal_capping.status") in done:
        add(11,
            f"Artificial chain termini introduced by missing internal regions "
            f"were neutralised by attaching ACE (acetyl) and NME "
            f"(N-methylamide) blocking groups."
        )

    prep_status = s("prepared_structure.status")
    prep_variant = s("prepared_structure.variant")
    available = s("available_models")
    if prep_status in done:
        variant_str = f" (variant: {prep_variant})" if prep_variant else ""
        avail_str = f" Available model families: {available}." if available else ""
        add(12,
            f"The final prepared structure was assembled{variant_str} "
            f"and written to the prepared/ directory.{avail_str}"
        )

    eval_status = s("model_eval.status")
    if eval_status in done:
        pct_fav = s("model_eval.rama_pct_favored")
        pct_out = s("model_eval.rama_pct_outlier")
        clash = s("model_eval.clashscore")
        quality = s("model_eval.overall_quality")
        filler_src = s("filler.model_source").lower()
        struct_note = (
            "the Modeller-filled structure" if "modeller" in filler_src
            else "the prepared structure"
        )
        add(13,
            f"Stereochemical quality of {struct_note} was assessed using "
            f"local PROCHECK{_cite('laskowski1993procheck')} and "
            f"MolProbity{_cite('lovell2003rama', 'williams2018molprobity')} "
            f"equivalents (see Ramachandran scatter plot). "
            f"Ramachandran: {pct_fav}% favored, {pct_out}% outlier. "
            f"Clashscore: {clash}. Overall quality: {quality}."
        )

    metals_status = s("metall_params.status")
    metals_sites = s("metall_params.site_count")
    if metals_status in done:
        site_str = f"{metals_sites} " if metals_sites else ""
        b3lyp_cite = _cite(
            "li2016mcpb", "becke1993b3lyp", "lee1988lyp",
            "hehre1972basis", "francl1982basis", "frisch2016gaussian",
        )
        if metals_status == "success":
            add(14,
                f"Metal-site AMBER-compatible parameters were derived for "
                f"{site_str}coordination site(s) using the MCPB.py "
                f"workflow{b3lyp_cite}. "
                f"Gaussian 16 single-point calculations were carried out at the "
                f"B3LYP{_cite('becke1993b3lyp', 'lee1988lyp')}/6-31G*"
                f"{_cite('hehre1972basis', 'francl1982basis')} level of theory "
                f"to derive bonded and non-bonded metal-site parameters. "
                f"Acceptance should be verified against the MCPB, Gaussian, "
                f"frcmod, and mol2 evidence paths rather than inferred from the "
                f"presence of a metal figure alone."
            )
        else:
            add(14,
                f"Metal-site parameterization evidence was recorded with warning "
                f"status for {site_str}coordination site(s). The B3LYP"
                f"{_cite('becke1993b3lyp', 'lee1988lyp')}/6-31G*"
                f"{_cite('hehre1972basis', 'francl1982basis')}/Gaussian 16"
                f"{_cite('frisch2016gaussian')} calculation via MCPB.py"
                f"{_cite('li2016mcpb')} should be reviewed before accepting "
                f"the generated force-field files."
            )

    nonstd_status = s("nonstd_residue_params.status")
    nonstd_n = s("nonstd_residue_params.n_residues")
    if nonstd_status in done:
        n_str = f"{nonstd_n} " if nonstd_n else ""
        add(15,
            f"AMBER-compatible force-field parameters were derived for "
            f"{n_str}non-standard residue type(s). "
            f"Atom types and bonded terms were assigned from the General AMBER "
            f"Force Field (GAFF){_cite('wang2004gaff')}. "
            f"Partial charges were computed by the RESP method"
            f"{_cite('bayly1993resp')} using electrostatic potential "
            f"grids calculated with Gaussian 16{_cite('frisch2016gaussian')} "
            f"at the HF/6-31G*{_cite('hehre1972basis', 'francl1982basis')} "
            f"level of theory, following the two-stage RESP fitting protocol."
        )

    return paras


# ---------------------------------------------------------------------------
# PyMOL launcher discovery
# ---------------------------------------------------------------------------

# Env vars required by the Debian pymol package (wrapper script values).
_DEBIAN_PYMOL_ENV: dict[str, str] = {
    "PYMOL_PATH":    "/usr/share/pymol",
    "PYMOL_DATA":    "/usr/share/pymol/data",
    "PYMOL_SCRIPTS": "/usr/share/pymol/scripts",
    "CHEMPY_DATA":   "/usr/share/pymol/data/chempy",
}

_pymol_cmd_cache: list[str] | None = None   # [python3_path, "-m", "pymol.__init__"]
_pymol_env_cache: dict[str, str] | None = None


def _find_pymol_launcher() -> tuple[list[str], dict[str, str]] | None:
    """
    Return ([python3, '-m', 'pymol.__init__'], extra_env) for a working PyMOL,
    or None if no working interpreter is found.

    Tries the current interpreter, system /usr/bin/python3, and 'python3' on
    PATH; each with and without the Debian env-var block.  Result is cached.
    """
    global _pymol_cmd_cache, _pymol_env_cache
    if _pymol_cmd_cache is not None:
        return _pymol_cmd_cache, _pymol_env_cache  # type: ignore[return-value]

    candidates = [sys.executable, "/usr/bin/python3", "python3"]
    for py in candidates:
        for extra in ({}, _DEBIAN_PYMOL_ENV):
            env = {**os.environ, **extra}
            try:
                r = subprocess.run(
                    [py, "-c", "import pymol"],
                    env=env, capture_output=True, timeout=15,
                )
                if r.returncode == 0:
                    cmd = [py, "-m", "pymol.__init__"]
                    _pymol_cmd_cache = cmd
                    _pymol_env_cache = extra
                    return cmd, extra
            except (FileNotFoundError, subprocess.TimeoutExpired):
                pass

    _pymol_cmd_cache = []   # sentinel: searched, not found
    _pymol_env_cache = {}
    return None


def _run_pml(pml_lines: list[str], output_png: Path, timeout: int = 180) -> str | None:
    """
    Write *pml_lines* to a temp .pml file, run PyMOL headless, and check that
    *output_png* was created.  Returns an error string or None on success.
    """
    launcher = _find_pymol_launcher()
    if launcher is None:
        return "PyMOL not found in any Python environment."
    cmd, extra_env = launcher

    env = {**os.environ, **extra_env}

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".pml", delete=False, encoding="utf-8"
    ) as tmp:
        tmp.write("\n".join(pml_lines) + "\n")
        pml_path = Path(tmp.name)

    try:
        subprocess.run(
            cmd + ["-cq", str(pml_path)],
            env=env, capture_output=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return f"PyMOL rendering timed out after {timeout} s."
    except Exception as exc:
        return f"PyMOL subprocess error: {exc}"
    finally:
        pml_path.unlink(missing_ok=True)

    if not output_png.exists():
        return "PyMOL ran but did not produce the expected PNG."
    return None


# ---------------------------------------------------------------------------
# PyMOL figure renderers
# ---------------------------------------------------------------------------

def _render_pymol_figure(
    protein_dir: Path,
    pdb_id: str,
    pipeline_record: dict[str, str],
    output_png: Path,
) -> str | None:
    """
    Render a PNG overview of the starting structure via headless PyMOL.

    Searches for {pdb_id}_delins.pdb first (post-insertion-code step), then
    the raw downloaded PDB.  Returns an error string on failure, None on success.
    """
    # Prefer the delins PDB (cleanest starting point post-renumbering)
    candidates = sorted(protein_dir.rglob(f"*{pdb_id.upper()}_delins.pdb"))
    if not candidates:
        candidates = sorted(protein_dir.rglob(f"*{pdb_id.lower()}_delins.pdb"))
    if not candidates:
        candidates = sorted(protein_dir.rglob(f"{pdb_id.upper()}.pdb"))
    if not candidates:
        candidates = sorted(protein_dir.rglob(f"{pdb_id.lower()}.pdb"))
    if not candidates:
        return f"No input PDB found for {pdb_id} in {protein_dir}."

    input_pdb = candidates[0]
    has_ligands = str(pipeline_record.get("has_ligands", "")).strip().lower()
    center_on_ligand = has_ligands in ("yes", "true")

    pml = [
        f"load {input_pdb}, s",
        "bg_color white",
        "hide everything, s",
        "set_color fruton_navy, [0.157, 0.196, 0.353]",
        "set_color fruton_gold, [0.541, 0.451, 0.059]",
        "set_color fruton_red, [0.902, 0.000, 0.078]",
        "set_color fruton_softgrey, [0.720, 0.735, 0.770]",
        "show cartoon, s and polymer",
        "color fruton_navy, s and ss h",
        "color fruton_gold, s and ss s",
        "color fruton_softgrey, s and ss l+''",
    ]
    if center_on_ligand:
        pml += [
            "select lig, s and hetatm and not (resn HOH or resn WAT)",
            "show sticks, lig",
            "color atomic, lig and not elem C",
            "color fruton_red, lig and elem C",
            "zoom lig, 12",
        ]
    else:
        pml += ["orient s", "zoom s"]
    pml += [
        "set ray_opaque_background, 1",
        "set antialias, 2",
        "set ray_shadows, 0",
        f"png {output_png}, width=1200, height=900, dpi=150, ray=1",
        "quit",
    ]
    return _run_pml(pml, output_png)


def _ion_element(ion_type: str) -> str:
    """
    Extract a 1–2 character PDB element symbol from an ion type string.

    Examples: 'Ca2+' → 'CA', 'Zn2+' → 'ZN', 'Fe' → 'FE', 'MG' → 'MG'.
    Used to build a reliable PyMOL ``elem`` selection that avoids the
    CA (alpha-carbon vs calcium) ambiguity of atom-name–based selections.
    """
    import re
    clean = re.sub(r"[0-9+\-\s]", "", ion_type).upper()
    return clean[:2] if clean else ion_type.upper()[:2]


def _render_metal_pocket_figures(
    pdb_path: Path,
    ion_type: str,
    output_dir: Path,
    pdb_id: str,
) -> dict[str, Path | None]:
    """
    Render two metal-pocket PNGs for one ion type in *pdb_path*:

    ``with_h``    — coordination shell with H visible (protonation state)
    ``distances`` — coordination shell without H, with labelled M-donor distances

    Returns dict with keys 'with_h' and 'distances'; value is Path or None.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    ion_elem = _ion_element(ion_type)   # reliable element symbol, e.g. 'CA', 'ZN'
    ion_tag  = ion_elem.replace(" ", "")

    png_h = output_dir / f"{pdb_id}_metal_{ion_tag}_with_h.png"
    png_d = output_dir / f"{pdb_id}_metal_{ion_tag}_distances.png"
    results: dict[str, Path | None] = {"with_h": None, "distances": None}

    # Coordination cutoffs: use 3.0 Å for residue shell, 2.8 Å for direct bonds.
    # 2.8 Å excludes beta-carbons (typically 3.1–3.3 Å from transition metals)
    # and water second-shell atoms while capturing all N/O/S/Se direct donors.
    COORD_SHELL_CUTOFF = 3.5   # residue selection: grab full sidechains
    DONOR_CUTOFF       = 2.8   # strict donor atoms (direct bonds only)

    # Common setup block — elem selector avoids the CA (alpha-C) / Ca2+ collision
    setup = [
        f"load {pdb_path}, s",
        "bg_color white",
        "hide everything",
        "set_color fruton_navy, [0.157, 0.196, 0.353]",
        "set_color fruton_gold, [0.541, 0.451, 0.059]",
        "set_color fruton_red, [0.902, 0.000, 0.078]",
        "set_color fruton_yellow, [0.980, 0.902, 0.275]",
        "set_color fruton_softgrey, [0.760, 0.770, 0.790]",
        # Thin transparent cartoon for context
        "show cartoon, s and polymer",
        "color fruton_softgrey, s and polymer",
        "set cartoon_transparency, 0.65",
        # Metal: HETATM record with matching element
        f"select metal_at, (s and hetatm and elem {ion_elem})",
        "show sphere, metal_at",
        "color fruton_red, metal_at",
        "set sphere_scale, 0.45, metal_at",
        # Coordinating residues: full sidechain within COORD_SHELL_CUTOFF of metal
        f"select coord_res, byres ((s and polymer) within {COORD_SHELL_CUTOFF} of metal_at)",
        "show sticks, coord_res",
        "color fruton_gold, coord_res and elem C",
        "color atomic, coord_res and not elem C",
        "center metal_at",
        "zoom (metal_at or coord_res), 3",
    ]

    # ---- Scene 1: with hydrogens ----
    pml_h = setup + [
        "show sticks, coord_res and elem H",
        "color white, coord_res and elem H",
        "set ray_opaque_background, 1",
        "set antialias, 2",
        "set ray_shadows, 0",
        f"png {png_h}, width=900, height=900, dpi=150, ray=1",
        "quit",
    ]
    err = _run_pml(pml_h, png_h)
    if err is None:
        results["with_h"] = png_h

    # ---- Scene 2: coordination distances ----
    # Select only proper donor atoms (N, O, S, Se) at strict cutoff.
    # Pass the same cutoff to `distance` so PyMOL only draws bonds, not all pairs.
    pml_d = setup + [
        "hide sticks, coord_res and elem H",
        # Direct donor heavy atoms (N/O/S/Se) within the strict bond cutoff.
        # NOTE: 'donors' is a reserved PyMOL keyword — use 'metal_donors'.
        (f"select metal_donors, (s within {DONOR_CUTOFF} of metal_at) and not metal_at "
         "and not elem H and (elem N or elem O or elem S or elem SE)"),
        # distance name, sel1, sel2, cutoff — cutoff prevents extra cross-pairs
        f"distance coord_dist, metal_at, metal_donors, {DONOR_CUTOFF}",
        "show dashes, coord_dist",
        "enable coord_dist",
        "set dash_gap, 0.15",
        "set dash_radius, 0.07",
        "set label_size, 16",
        "set label_color, black",
        "color fruton_yellow, coord_dist",
        "set ray_opaque_background, 1",
        "set antialias, 2",
        "set ray_shadows, 0",
        f"png {png_d}, width=900, height=900, dpi=150, ray=1",
        "quit",
    ]
    err = _run_pml(pml_d, png_d)
    if err is None:
        results["distances"] = png_d

    return results


# ---------------------------------------------------------------------------
# ReportLab visual helpers
# ---------------------------------------------------------------------------

def _safe_pdf_text(value: Any) -> str:
    """Return text escaped for ReportLab Paragraph input."""
    return escape(str(value or "").strip())


def _record_text(record: dict[str, str], key: str, default: str = "") -> str:
    """Read one pipeline-record value as stripped text."""
    value = str(record.get(key, default) or default).strip()
    return value


def _boolish_yes(value: str) -> bool:
    """Return True for the yes/true/1 forms used in flattened pipeline state."""
    return value.strip().lower() in {"yes", "true", "1", "y"}


def _read_fasta_records(fasta_path: Path) -> list[dict[str, str]]:
    """Read FASTA entries from *fasta_path* without adding a BioPython dependency here.

    The report module is intentionally a presentation layer.  It should consume the
    files produced by ``fasta_files.py`` and ``sequence_alignment.py`` but should not
    re-run sequence logic or assume one exact file naming convention.  A tiny FASTA
    parser is sufficient for report evidence because it only needs headers,
    sequences, lengths, and source paths.
    """
    records: list[dict[str, str]] = []
    header = ""
    seq_parts: list[str] = []

    try:
        lines = fasta_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return records

    def flush() -> None:
        nonlocal header, seq_parts
        if not header and not seq_parts:
            return
        sequence = "".join(part.strip() for part in seq_parts if part.strip()).replace(" ", "")
        if sequence:
            records.append({
                "header": header.strip() or fasta_path.stem,
                "sequence": sequence,
                "path": str(fasta_path),
            })
        header = ""
        seq_parts = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            flush()
            header = stripped[1:].strip()
        else:
            seq_parts.append(stripped)
    flush()
    return records


def _classify_fasta_source(path: Path, header: str) -> str:
    """Return a human-readable source label for a FASTA record."""
    text = f"{path.name} {header}".lower()
    if "uniprot" in text:
        return "UniProt canonical sequence"
    if "seqres" in text:
        return "PDB SEQRES deposited sequence"
    if "atom" in text:
        return "PDB ATOM resolved sequence"
    if "alignment" in [part.lower() for part in path.parts] or "aln" in text:
        return "Alignment FASTA"
    return "FASTA sequence"


def _wrap_sequence_for_pdf(sequence: str, line_width: int = 60, max_chars: int = 1400) -> tuple[str, bool]:
    """Wrap a sequence for a ReportLab paragraph and report whether it was truncated."""
    clean = "".join(sequence.split())
    truncated = len(clean) > max_chars
    if truncated:
        keep_tail = min(180, max_chars // 5)
        keep_head = max_chars - keep_tail
        clean = clean[:keep_head] + "..." + clean[-keep_tail:]
    lines = [clean[i:i + line_width] for i in range(0, len(clean), line_width)] or [""]
    return "<br/>".join(_safe_pdf_text(line) for line in lines), truncated


def _discover_fasta_records(protein_dir: Path, pdb_id: str, max_records: int = 8) -> list[dict[str, str]]:
    """Find the most useful FASTA records produced for one FRUTON protein folder."""
    fasta_roots = [protein_dir / "fasta", protein_dir]
    candidates: list[Path] = []
    seen: set[Path] = set()

    for root in fasta_roots:
        if not root.exists():
            continue
        for pattern in ("*.fasta", "*.fa", "*.faa"):
            for hit in root.rglob(pattern):
                lower_name = hit.name.lower()
                lower_parts = {part.lower() for part in hit.parts}
                if "alignments" in lower_parts:
                    continue
                if lower_name.endswith(".aln.fasta") or ".aln." in lower_name:
                    continue
                if "alignment" in lower_name or "_vs_" in lower_name:
                    continue
                try:
                    resolved = hit.resolve()
                except OSError:
                    resolved = hit
                if resolved not in seen:
                    seen.add(resolved)
                    candidates.append(hit)

    def sort_key(path: Path) -> tuple[int, str]:
        name = path.name.lower()
        if "uniprot" in name:
            priority = 0
        elif "atom" in name:
            priority = 1
        elif "seqres" in name:
            priority = 2
        else:
            priority = 3
        return (priority, str(path))

    records: list[dict[str, str]] = []
    for fasta_path in sorted(candidates, key=sort_key):
        for rec in _read_fasta_records(fasta_path):
            rec_path = Path(rec["path"])
            rec["source"] = _classify_fasta_source(rec_path, rec["header"])
            rec["length"] = str(len(rec["sequence"]))
            rec["unknown_x"] = str(rec["sequence"].upper().count("X"))
            records.append(rec)
            if len(records) >= max_records:
                return records
    return records


def _discover_alignment_assets(protein_dir: Path, max_images: int = 3) -> dict[str, list[Path]]:
    """Return alignment PNG and mapping TSV evidence generated by FRUTON."""
    alignment_roots = [protein_dir / "fasta" / "alignments", protein_dir / "alignments"]
    images: list[Path] = []
    mappings: list[Path] = []
    fastas: list[Path] = []

    for root in alignment_roots:
        if not root.exists():
            continue
        images.extend(sorted(root.glob("*.png")))
        mappings.extend(sorted(root.glob("*.mapping.tsv")))
        mappings.extend(sorted(root.glob("*mapping*.tsv")))
        fastas.extend(sorted(root.glob("*.aln.fasta")))

    # De-duplicate while preserving order.
    def unique(paths: list[Path]) -> list[Path]:
        out: list[Path] = []
        seen: set[Path] = set()
        for path in paths:
            try:
                resolved = path.resolve()
            except OSError:
                resolved = path
            if resolved not in seen:
                seen.add(resolved)
                out.append(path)
        return out

    return {
        "images": unique(images)[:max_images],
        "mappings": unique(mappings),
        "alignments": unique(fastas),
    }


def _short_path_for_report(path: str | Path, protein_dir: Path) -> str:
    """Return a compact path string for PDF tables."""
    p = Path(path)
    try:
        return str(p.relative_to(protein_dir))
    except ValueError:
        try:
            return str(p.relative_to(protein_dir.parent))
        except ValueError:
            return str(p)


def _collect_path_evidence(record: dict[str, str], protein_dir: Path, max_rows: int = 18) -> list[tuple[str, str]]:
    """Collect non-empty path-like fields from the flattened FRUTON record."""
    preferred = [
        "pdb_path",
        "processed_pdb_path",
        "fasta.seqres_path",
        "fasta.atom_path",
        "fasta.uniprot_path",
        "sequence_alignment.mapping_path",
        "components_dir",
        "components.protein_path",
        "components.water_path",
        "components.ligand_path",
        "components.metals_path",
        "filler.final_model_path",
        "protonation.output_path",
        "prepared_structure.output_path",
        "model_eval.rama_plot_path",
        "metall_params.output_dir",
        "metall_params.frcmod_path",
        "metall_params.mol2_path",
        "metall_params.gaussian_input_path",
        "metall_params.gaussian_log_path",
    ]
    rows: list[tuple[str, str]] = []
    seen_keys: set[str] = set()

    for key in preferred:
        value = _record_text(record, key)
        if value:
            rows.append((key, _short_path_for_report(value, protein_dir)))
            seen_keys.add(key)

    for key in sorted(record):
        if key in seen_keys:
            continue
        key_l = key.lower()
        if not any(token in key_l for token in ("path", "dir")):
            continue
        value = _record_text(record, key)
        if not value:
            continue
        rows.append((key, _short_path_for_report(value, protein_dir)))
        if len(rows) >= max_rows:
            break

    return rows[:max_rows]


def _collect_prefixed_evidence(record: dict[str, str], prefixes: tuple[str, ...]) -> list[tuple[str, str]]:
    """Collect non-empty flattened record values matching any prefix."""
    rows: list[tuple[str, str]] = []
    for key in sorted(record):
        if not key.startswith(prefixes):
            continue
        value = _record_text(record, key)
        if value:
            rows.append((key, value))
    return rows


def _preparation_decision(record: dict[str, str]) -> tuple[str, str]:
    """Return a conservative audit decision and its reason."""
    prep_status = _record_text(record, "prepared_structure.status").lower()
    quality = _record_text(record, "model_eval.overall_quality").lower()
    metals = _record_text(record, "has_metals").lower()
    metal_status = _record_text(record, "metall_params.status").lower()
    nonstd = _record_text(record, "has_nonstandard_residues").lower()
    nonstd_status = _record_text(record, "nonstd_residue_params.status").lower()

    if prep_status in {"failed", "failure", "error"}:
        return "blocked", "prepared structure generation failed"
    if metals in {"yes", "true"} and metal_status not in {"success", "ready", "passed"}:
        return "requires metal-parameter review", "metal detected without accepted metal-parameter evidence"
    if nonstd in {"yes", "true"} and nonstd_status not in {"success", "ready", "passed"}:
        return "requires non-standard-residue review", "non-standard residue detected without accepted parameter evidence"
    if quality in {"poor", "bad", "failed", "fail"}:
        return "requires structure-quality review", "Ramachandran/clash metrics are outside the preferred range"
    if prep_status in {"success", "warning"}:
        return "prepared with caveats", "review the evidence tables before downstream MD/GBSA"
    return "not enough evidence", "prepared-structure status or validation metrics are missing"


def _asset_search_roots(protein_dir: Path) -> list[Path]:
    """Return likely asset directories for report logos.

    The report generator should not require a package-data installation just to
    render a PDF.  This helper therefore checks local report/protein folders,
    the canonical project data folder (`data/logo_idis.png`), common project
    asset folders, the current working directory, the module directory, and
    `/mnt/data` for notebook/sandbox runs.  Missing logos are not fatal: the page
    header falls back to text labels, preserving the public API.
    """
    module_dir = Path(__file__).resolve().parent
    roots = [
        protein_dir / "report",
        protein_dir,
        protein_dir.parent,
        protein_dir.parent / "assets",
        protein_dir.parent.parent,
        protein_dir.parent.parent / "assets",
        protein_dir.parent.parent / "data",
        Path.cwd(),
        Path.cwd() / "data",
        Path.cwd() / "assets",
        module_dir,
        module_dir / "assets",
        module_dir.parent / "assets",
        module_dir.parent / "data",
        Path("/mnt/data"),
    ]
    seen: set[Path] = set()
    unique: list[Path] = []
    for root in roots:
        try:
            resolved = root.resolve()
        except OSError:
            resolved = root
        if resolved not in seen:
            seen.add(resolved)
            unique.append(root)
    return unique


def _find_report_asset(
    *,
    env_var: str,
    names: tuple[str, ...],
    protein_dir: Path,
) -> Path | None:
    """Find the first existing image asset from an env var or standard names."""
    env_value = os.environ.get(env_var, "").strip()
    if env_value:
        env_path = Path(env_value).expanduser()
        if env_path.exists():
            return env_path

    for root in _asset_search_roots(protein_dir):
        for name in names:
            candidate = root / name
            if candidate.exists():
                return candidate
    return None


def _status_palette_for_value(field: str, value: str) -> tuple[str, str]:
    """Return background/text colours for compact audit cells.

    Status fields are process outcomes; semantic fields are domain flags.  This
    distinction follows the FRUTON workbook/audit policy: ligands are mild,
    metals and non-standard residues require review, and gaps are structural
    caution rather than automatic failure.
    """
    raw = value.strip().lower()
    field_l = field.strip().lower()

    if field_l in {"has_ligands", "ligands"}:
        return (_FRUTON_GREEN_SOFT, _IDIS_TEXT) if raw in {"yes", "true"} else ("#FFFFFF", _IDIS_MUTED)
    if field_l in {"has_metals", "metals", "has_nonstandard_residues", "nonstandard"}:
        return (_AUDIT_FAIL, _IDIS_TEXT) if raw in {"yes", "true"} else (_AUDIT_OK, _IDIS_TEXT)
    if field_l in {"has_gaps", "gaps"}:
        return (_AUDIT_WARN, _IDIS_TEXT) if raw in {"yes", "true"} else (_AUDIT_OK, _IDIS_TEXT)
    if field_l in {"overall_quality", "quality"}:
        if raw in {"good", "excellent", "acceptable", "pass", "passed"}:
            return (_AUDIT_OK, _IDIS_TEXT)
        if raw in {"poor", "bad", "failed", "fail"}:
            return (_AUDIT_FAIL, _IDIS_TEXT)
        if raw:
            return (_AUDIT_WARN, _IDIS_TEXT)
        return (_AUDIT_GREY, _IDIS_MUTED)

    if raw in {"success", "passed", "pass", "ready", "yes", "true"}:
        return (_AUDIT_OK, _IDIS_TEXT)
    if raw in {"warning", "required", "review", "caution"}:
        return (_AUDIT_WARN, _IDIS_TEXT)
    if raw in {"failed", "failure", "error", "no", "false"}:
        return (_AUDIT_FAIL, _IDIS_TEXT)
    if raw in {"skipped", "skip", "none", ""}:
        return (_AUDIT_GREY, _IDIS_MUTED)
    return (_IDIS_PANEL, _IDIS_TEXT)


def _quality_summary(record: dict[str, str]) -> str:
    """Build one compact model-quality sentence for the title card."""
    favored = _record_text(record, "model_eval.rama_pct_favored")
    outlier = _record_text(record, "model_eval.rama_pct_outlier")
    clash = _record_text(record, "model_eval.clashscore")
    quality = _record_text(record, "model_eval.overall_quality")
    parts: list[str] = []
    if quality:
        parts.append(f"quality: {quality}")
    if favored:
        parts.append(f"Rama favored: {favored}%")
    if outlier:
        parts.append(f"outlier: {outlier}%")
    if clash:
        parts.append(f"clashscore: {clash}")
    return " · ".join(parts) if parts else "quality metrics not recorded"


def _make_page_chrome(
    *,
    pdb_id: str,
    idis_logo: Path | None,
    fruton_logo: Path | None,
):
    """Return a ReportLab page callback with IDIS/FRUTON header and footer."""
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import cm

    page_w, page_h = A4

    def _draw(canvas, doc) -> None:  # type: ignore[no-untyped-def]
        canvas.saveState()
        left = doc.leftMargin
        right = page_w - doc.rightMargin
        top = page_h - 0.75 * cm

        if idis_logo and idis_logo.exists():
            canvas.drawImage(
                str(idis_logo),
                left,
                page_h - 1.45 * cm,
                width=2.1 * cm,
                height=0.85 * cm,
                preserveAspectRatio=True,
                mask="auto",
            )
        else:
            canvas.setFillColor(colors.HexColor(_IDIS_NAVY))
            canvas.setFont("Helvetica-Bold", 18)
            canvas.drawString(left, page_h - 1.20 * cm, "idis")

        canvas.setFillColor(colors.HexColor(_IDIS_NAVY))
        canvas.setFont("Helvetica-Bold", 8.5)
        canvas.drawCentredString(page_w / 2, top - 0.05 * cm, "FRUTON protein preparation audit")
        canvas.setFillColor(colors.HexColor(_IDIS_MUTED))
        canvas.setFont("Helvetica", 6.6)
        canvas.drawCentredString(page_w / 2, top - 0.35 * cm, f"per-protein report · {pdb_id}")

        if fruton_logo and fruton_logo.exists():
            canvas.drawImage(
                str(fruton_logo),
                right - 1.45 * cm,
                page_h - 1.45 * cm,
                width=1.45 * cm,
                height=0.95 * cm,
                preserveAspectRatio=True,
                mask="auto",
            )
        else:
            canvas.setFillColor(colors.HexColor(_FRUTON_RED))
            canvas.setFont("Helvetica-Bold", 9)
            canvas.drawRightString(right, page_h - 1.05 * cm, "FRUTON")

        y = page_h - 1.63 * cm
        span = right - left
        canvas.setLineWidth(1.0)
        canvas.setStrokeColor(colors.HexColor(_FRUTON_RED))
        canvas.line(left, y, left + 0.25 * span, y)
        canvas.setStrokeColor(colors.HexColor(_IDIS_GOLD))
        canvas.line(left + 0.26 * span, y, left + 0.49 * span, y)
        canvas.setStrokeColor(colors.HexColor(_IDIS_NAVY))
        canvas.line(left + 0.50 * span, y, right, y)

        footer_y = 0.95 * cm
        canvas.setStrokeColor(colors.HexColor(_IDIS_LINE))
        canvas.setLineWidth(0.45)
        canvas.line(left, footer_y + 0.22 * cm, right, footer_y + 0.22 * cm)
        canvas.setFillColor(colors.HexColor(_IDIS_MUTED))
        canvas.setFont("Helvetica", 6.5)
        canvas.drawString(left, footer_y, f"FRUTON - {_FRUTON_EXPANSION}")
        canvas.drawRightString(right, footer_y, f"page {doc.page}")
        canvas.restoreState()

    return _draw


def _section_block(title: str, subtitle: str | None, styles: dict[str, Any], width: float):
    """Build a compact section title with a FRUTON red side rule."""
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph, Table, TableStyle

    content = [Paragraph(_safe_pdf_text(title), styles["SectionHeading"])]
    if subtitle:
        content.append(Paragraph(_safe_pdf_text(subtitle), styles["SectionSubtitle"]))
    block = Table([[[*content]]], colWidths=[width])
    block.setStyle(TableStyle([
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LINEBEFORE", (0, 0), (0, 0), 3.0, colors.HexColor(_FRUTON_RED)),
    ]))
    return block


def _build_summary_table(record: dict[str, str], styles: dict[str, Any], width: float):
    """Build the per-protein triage card shown on the first page."""
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph, Table, TableStyle

    variant = _record_text(record, "prepared_structure.variant") or _record_text(record, "prepared_structure.final_type") or "not recorded"
    prep_status = _record_text(record, "prepared_structure.status") or "not recorded"
    models = _record_text(record, "available_models") or _record_text(record, "final_model_type") or "not recorded"
    n_gaps = _record_text(record, "n_gaps")
    gap_sizes = _record_text(record, "gap_sizes")
    gaps_value = "yes" if (n_gaps and n_gaps != "0") or _boolish_yes(_record_text(record, "has_gaps")) else "no"
    gaps_display = f"{n_gaps or '0'}" + (f" ({gap_sizes})" if gap_sizes else "")
    metal_value = _record_text(record, "has_metals") or "no"
    ligand_value = _record_text(record, "has_ligands") or "no"
    nonstd_value = _record_text(record, "has_nonstandard_residues") or "no"
    quality = _record_text(record, "model_eval.overall_quality") or "not recorded"

    rows = [
        ["prepared", prep_status, "variant", variant],
        ["models", models, "gaps", gaps_display],
        ["ligands", ligand_value, "metals", metal_value],
        ["nonstandard", nonstd_value, "quality", quality],
    ]
    data = []
    for label_a, value_a, label_b, value_b in rows:
        data.append([
            Paragraph(_safe_pdf_text(label_a.upper()), styles["MetricLabel"]),
            Paragraph(_safe_pdf_text(value_a), styles["MetricValue"]),
            Paragraph(_safe_pdf_text(label_b.upper()), styles["MetricLabel"]),
            Paragraph(_safe_pdf_text(value_b), styles["MetricValue"]),
        ])

    table = Table(data, colWidths=[2.35 * width / 14, 4.65 * width / 14, 2.35 * width / 14, 4.65 * width / 14])
    style_cmds: list[tuple[Any, ...]] = [
        ("BOX", (0, 0), (-1, -1), 0.75, colors.HexColor(_IDIS_LINE)),
        ("INNERGRID", (0, 0), (-1, -1), 0.35, colors.HexColor(_IDIS_LINE)),
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#FFFFFF")),
        ("BACKGROUND", (0, 0), (0, -1), colors.HexColor(_IDIS_PANEL_2)),
        ("BACKGROUND", (2, 0), (2, -1), colors.HexColor(_IDIS_PANEL_2)),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]
    semantic_cells = {
        (1, 0): ("prepared", prep_status),
        (3, 1): ("gaps", gaps_value),
        (1, 2): ("has_ligands", ligand_value),
        (3, 2): ("has_metals", metal_value),
        (1, 3): ("has_nonstandard_residues", nonstd_value),
        (3, 3): ("quality", quality),
    }
    for cell, (field, value) in semantic_cells.items():
        bg, fg = _status_palette_for_value(field, value)
        style_cmds.append(("BACKGROUND", cell, cell, colors.HexColor(bg)))
        style_cmds.append(("TEXTCOLOR", cell, cell, colors.HexColor(fg)))
    table.setStyle(TableStyle(style_cmds))
    return table


def _build_narrative_table(paragraphs: list[Any], width: float):
    """Wrap pipeline narrative paragraphs in a light audit panel."""
    from reportlab.lib import colors
    from reportlab.platypus import Table, TableStyle

    table = Table([[paragraphs]], colWidths=[width])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#FFFFFF")),
        ("BOX", (0, 0), (-1, -1), 0.75, colors.HexColor(_IDIS_LINE)),
        ("LINEBEFORE", (0, 0), (-1, -1), 2.0, colors.HexColor(_IDIS_GOLD)),
        ("LEFTPADDING", (0, 0), (-1, -1), 10),
        ("RIGHTPADDING", (0, 0), (-1, -1), 10),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
    ]))
    return table

def _build_simple_evidence_table(
    rows: list[tuple[str, str]],
    styles: dict[str, Any],
    width: float,
    *,
    left_header: str = "Field",
    right_header: str = "Value",
):
    """Build a two-column evidence table for paths, parameters, and decisions."""
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph, Table, TableStyle

    if not rows:
        rows = [("status", "no evidence recorded")]

    data = [[
        Paragraph(_safe_pdf_text(left_header), styles["EvidenceHeader"]),
        Paragraph(_safe_pdf_text(right_header), styles["EvidenceHeader"]),
    ]]
    for key, value in rows:
        data.append([
            Paragraph(_safe_pdf_text(key), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(value), styles["EvidenceCell"]),
        ])

    table = Table(data, colWidths=[0.32 * width, 0.68 * width], repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(_IDIS_PANEL_2)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor(_IDIS_NAVY)),
        ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(_IDIS_LINE)),
        ("INNERGRID", (0, 0), (-1, -1), 0.3, colors.HexColor(_IDIS_LINE)),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    return table


def _collect_residue_renames(
    input_pdb: Path,
    output_pdb: Path,
) -> list[dict[str, str]]:
    """Return per-residue rename records by diffing protonation input vs output.

    Both files are read and residues are matched by (chain, resseq, icode).
    A record is emitted only when the residue name changes. Returns an empty
    list when either file is missing or unreadable.

    Each record has keys: chain, res_num, icode, from_name, to_name, source.
    source is 'PROPKA' for HIP assignments, 'GROMACS' for others (e.g.
    disulfide CYS2, MSE->MET).
    """
    def _parse_residues(pdb_path: Path) -> dict[tuple[str, int, str], str]:
        res: dict[tuple[str, int, str], str] = {}
        try:
            with pdb_path.open("r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
                        continue
                    try:
                        res_name = line[17:20].strip()
                        chain = line[21]
                        res_num = int(line[22:26])
                        icode = line[26].strip() if len(line) > 26 else ""
                        key = (chain, res_num, icode)
                        if key not in res:
                            res[key] = res_name
                    except (ValueError, IndexError):
                        pass
        except OSError:
            pass
        return res

    if not input_pdb.is_file() or not output_pdb.is_file():
        return []

    before = _parse_residues(input_pdb)
    after = _parse_residues(output_pdb)

    _HIP_KEYS = {"HIP"}
    renames: list[dict[str, str]] = []
    for key, from_name in sorted(before.items(), key=lambda kv: kv[0]):
        to_name = after.get(key)
        if to_name is None or to_name == from_name:
            continue
        chain, res_num, icode = key
        source = "PROPKA" if to_name in _HIP_KEYS else "GROMACS"
        renames.append({
            "chain": chain,
            "res_num": str(res_num),
            "icode": icode or "",
            "from_name": from_name,
            "to_name": to_name,
            "source": source,
        })
    return renames


def _build_residue_rename_table(
    renames: list[dict[str, str]],
    styles: dict[str, Any],
    width: float,
) -> Any:
    """Build a table of protonation-step residue renamings."""
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph, Table, TableStyle

    data: list[list[Any]] = [[
        Paragraph("Chain", styles["EvidenceHeader"]),
        Paragraph("Res.", styles["EvidenceHeader"]),
        Paragraph("Ins.", styles["EvidenceHeader"]),
        Paragraph("From", styles["EvidenceHeader"]),
        Paragraph("To", styles["EvidenceHeader"]),
        Paragraph("Driver", styles["EvidenceHeader"]),
        Paragraph("Note", styles["EvidenceHeader"]),
    ]]

    _notes = {
        "HIP": "doubly protonated His; net charge +1",
        "CYS2": "disulfide-bridged Cys",
        "MET": "selenomethionine converted to methionine",
    }

    propka_rows: set[int] = set()
    for r in renames:
        row_idx = len(data)
        note = _notes.get(r["to_name"], "")
        data.append([
            Paragraph(_safe_pdf_text(r["chain"]), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(r["res_num"]), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(r["icode"] or "-"), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(r["from_name"]), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(r["to_name"]), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(r["source"]), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(note), styles["EvidenceCell"]),
        ])
        if r["source"] == "PROPKA":
            propka_rows.add(row_idx)

    col_w = width
    col_widths = [
        0.07 * col_w,   # Chain
        0.07 * col_w,   # Res.
        0.06 * col_w,   # Ins.
        0.10 * col_w,   # From
        0.10 * col_w,   # To
        0.12 * col_w,   # Driver
        0.48 * col_w,   # Note
    ]

    table = Table(data, colWidths=col_widths, repeatRows=1)
    style_cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(_IDIS_PANEL_2)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor(_IDIS_NAVY)),
        ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(_IDIS_LINE)),
        ("INNERGRID", (0, 0), (-1, -1), 0.3, colors.HexColor(_IDIS_LINE)),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]
    for row_idx in propka_rows:
        style_cmds.append(
            ("BACKGROUND", (0, row_idx), (-1, row_idx), colors.HexColor(_FRUTON_YELLOW_SOFT))
        )
    table.setStyle(TableStyle(style_cmds))
    return table


def _build_sequence_evidence_flowables(
    fasta_records: list[dict[str, str]],
    styles: dict[str, Any],
    width: float,
    protein_dir: Path,
) -> list[Any]:
    """Build flowables showing FASTA sources and wrapped sequence content."""
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph, Spacer, Table, TableStyle
    from reportlab.lib.units import cm

    flowables: list[Any] = []
    summary_rows: list[list[Any]] = [[
        Paragraph("Source", styles["EvidenceHeader"]),
        Paragraph("Record", styles["EvidenceHeader"]),
        Paragraph("Length", styles["EvidenceHeader"]),
        Paragraph("X", styles["EvidenceHeader"]),
        Paragraph("File", styles["EvidenceHeader"]),
    ]]

    for idx, rec in enumerate(fasta_records, start=1):
        header = rec.get("header", "")
        short_header = header if len(header) <= 70 else header[:67] + "..."
        summary_rows.append([
            Paragraph(_safe_pdf_text(rec.get("source", "FASTA")), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(short_header), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(rec.get("length", "")), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(rec.get("unknown_x", "0")), styles["EvidenceCell"]),
            Paragraph(_safe_pdf_text(_short_path_for_report(rec.get("path", ""), protein_dir)), styles["EvidenceCell"]),
        ])

    summary = Table(summary_rows, colWidths=[0.23 * width, 0.30 * width, 0.08 * width, 0.05 * width, 0.34 * width], repeatRows=1)
    summary.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(_IDIS_PANEL_2)),
        ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(_IDIS_LINE)),
        ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.HexColor(_IDIS_LINE)),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
    ]))
    flowables.extend([summary, Spacer(1, 0.13 * cm)])

    for idx, rec in enumerate(fasta_records, start=1):
        wrapped, truncated = _wrap_sequence_for_pdf(rec.get("sequence", ""))
        trunc_note = " — sequence truncated for report length" if truncated else ""
        title = (
            f"FASTA {idx}: {rec.get('source', 'FASTA')} · "
            f"length {rec.get('length', '')} · X={rec.get('unknown_x', '0')}{trunc_note}"
        )
        seq_block = Table([[
            [
                Paragraph(_safe_pdf_text(title), styles["EvidenceHeader"]),
                Paragraph(wrapped, styles["MonoSequence"]),
            ]
        ]], colWidths=[width])
        seq_block.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_IDIS_PANEL)),
            ("BOX", (0, 0), (-1, -1), 0.45, colors.HexColor(_IDIS_LINE)),
            ("LEFTPADDING", (0, 0), (-1, -1), 7),
            ("RIGHTPADDING", (0, 0), (-1, -1), 7),
            ("TOPPADDING", (0, 0), (-1, -1), 5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ]))
        flowables.extend([seq_block, Spacer(1, 0.10 * cm)])

    return flowables


# ---------------------------------------------------------------------------
# PDF builder (lazy reportlab import)
# ---------------------------------------------------------------------------

def _build_pdf(
    pdb_id: str,
    pipeline_record: dict[str, str],
    output_pdf: Path,
    figure_png: Path | None,
    ref_keys: list[str],
    metal_png_h: Path | None = None,
    metal_png_dist: Path | None = None,
    ion_type: str = "",
    pdb_citation_inline: str = "",
) -> None:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import cm
    from reportlab.platypus import (
        HRFlowable,
        Image,
        KeepTogether,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    report_dir = output_pdf.parent
    protein_dir = report_dir.parent
    idis_logo = _find_report_asset(
        env_var="IDIS_LOGO_PATH",
        names=_IDIS_LOGO_NAMES,
        protein_dir=protein_dir,
    )
    fruton_logo = _find_report_asset(
        env_var="FRUTON_LOGO_PATH",
        names=_FRUTON_LOGO_NAMES,
        protein_dir=protein_dir,
    )

    doc = SimpleDocTemplate(
        str(output_pdf),
        pagesize=A4,
        leftMargin=1.75 * cm,
        rightMargin=1.75 * cm,
        topMargin=2.25 * cm,
        bottomMargin=1.70 * cm,
        title=f"FRUTON protein preparation audit: {pdb_id}",
        author="FRUTON pipeline",
        subject="Per-protein preparation report",
    )
    usable_width = A4[0] - doc.leftMargin - doc.rightMargin
    styles = getSampleStyleSheet()

    styles.add(ParagraphStyle(
        "Eyebrow",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=7.5,
        leading=9,
        textColor=colors.HexColor(_IDIS_GOLD),
        spaceAfter=2,
        uppercase=True,
    ))
    styles.add(ParagraphStyle(
        "ReportTitle",
        parent=styles["Title"],
        fontName="Helvetica-Bold",
        fontSize=21,
        leading=24,
        alignment=0,
        spaceBefore=0,
        spaceAfter=4,
        textColor=colors.HexColor(_IDIS_NAVY),
    ))
    styles.add(ParagraphStyle(
        "ReportSubtitle",
        parent=styles["Normal"],
        fontSize=8.8,
        leading=12,
        textColor=colors.HexColor(_IDIS_MUTED),
        spaceAfter=8,
    ))
    styles.add(ParagraphStyle(
        "SectionHeading",
        parent=styles["Heading2"],
        fontName="Helvetica-Bold",
        fontSize=12.4,
        leading=15,
        spaceBefore=0,
        spaceAfter=1,
        textColor=colors.HexColor(_IDIS_NAVY),
    ))
    styles.add(ParagraphStyle(
        "SectionSubtitle",
        parent=styles["Normal"],
        fontSize=7.8,
        leading=10,
        textColor=colors.HexColor(_IDIS_MUTED),
        spaceAfter=0,
    ))
    styles.add(ParagraphStyle(
        "BodyAudit",
        parent=styles["Normal"],
        fontSize=9.2,
        leading=13.4,
        spaceAfter=6,
        textColor=colors.HexColor(_IDIS_TEXT),
    ))
    styles.add(ParagraphStyle(
        "CaptionAudit",
        parent=styles["Normal"],
        fontSize=7.5,
        leading=9.5,
        textColor=colors.HexColor(_IDIS_MUTED),
        spaceAfter=2,
    ))
    styles.add(ParagraphStyle(
        "ReferenceAudit",
        parent=styles["Normal"],
        fontSize=7.1,
        leading=9.2,
        leftIndent=8,
        firstLineIndent=-8,
        spaceAfter=2.4,
        textColor=colors.HexColor(_IDIS_TEXT),
    ))
    styles.add(ParagraphStyle(
        "MetricLabel",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=6.5,
        leading=8,
        textColor=colors.HexColor(_IDIS_MUTED),
    ))
    styles.add(ParagraphStyle(
        "MetricValue",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=8.3,
        leading=10.5,
        textColor=colors.HexColor(_IDIS_TEXT),
    ))
    styles.add(ParagraphStyle(
        "SmallNote",
        parent=styles["Normal"],
        fontSize=7.4,
        leading=9.5,
        textColor=colors.HexColor(_IDIS_MUTED),
    ))

    styles.add(ParagraphStyle(
        "MonoSequence",
        parent=styles["Normal"],
        fontName="Courier",
        fontSize=6.7,
        leading=8.2,
        textColor=colors.HexColor(_IDIS_TEXT),
        spaceAfter=3,
    ))
    styles.add(ParagraphStyle(
        "EvidenceHeader",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=7.0,
        leading=8.5,
        textColor=colors.HexColor(_IDIS_NAVY),
    ))
    styles.add(ParagraphStyle(
        "EvidenceCell",
        parent=styles["Normal"],
        fontSize=7.0,
        leading=9.0,
        textColor=colors.HexColor(_IDIS_TEXT),
    ))

    thin_rule = HRFlowable(
        width="100%",
        thickness=0.6,
        color=colors.HexColor(_IDIS_LINE),
        spaceBefore=8,
        spaceAfter=6,
    )

    story: list[Any] = []

    # ------------------------------------------------------------------
    # Title and summary block
    # ------------------------------------------------------------------
    story.append(Paragraph("FRUTON · IDIS visual report", styles["Eyebrow"]))
    story.append(Paragraph(f"Protein preparation audit: {_safe_pdf_text(pdb_id)}", styles["ReportTitle"]))

    uniprot = _record_text(pipeline_record, "uniprot_id")
    residue_range = _record_text(pipeline_record, "range")
    subtitle_parts = []
    if uniprot:
        subtitle_parts.append(f"UniProt: {_safe_pdf_text(uniprot)}")
    if residue_range:
        subtitle_parts.append(f"Range: {_safe_pdf_text(residue_range)}")
    subtitle_parts.append(_safe_pdf_text(_quality_summary(pipeline_record)))
    story.append(Paragraph(" &nbsp;·&nbsp; ".join(subtitle_parts), styles["ReportSubtitle"]))

    story.append(_build_summary_table(pipeline_record, styles, usable_width))
    story.append(Spacer(1, 0.32 * cm))

    # ------------------------------------------------------------------
    # Structure figure
    # ------------------------------------------------------------------
    story.append(_section_block(
        "Starting structure",
        "Initial coordinate evidence rendered for visual inspection before chemistry-specific interpretation.",
        styles,
        usable_width,
    ))
    story.append(Spacer(1, 0.08 * cm))
    if figure_png and figure_png.exists():
        ligand_note = (
            "centered on ligand"
            if _boolish_yes(_record_text(pipeline_record, "has_ligands"))
            else "default orientation"
        )
        img = Image(str(figure_png), width=13.7 * cm, height=10.25 * cm)
        fig_table = Table([[img]], colWidths=[usable_width])
        fig_table.setStyle(TableStyle([
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_IDIS_PANEL)),
            ("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(_IDIS_LINE)),
            ("TOPPADDING", (0, 0), (-1, -1), 8),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
        ]))
        story.append(fig_table)
        story.append(Paragraph(
            f"<i>Figure 1. PyMOL{_cite('schrodinger2015pymol')} rendering of the starting "
            f"structure ({_safe_pdf_text(pdb_id)}, {_safe_pdf_text(ligand_note)}). "
            "Helices use the IDIS navy; strands use the IDIS gold accent; "
            "loops are muted grey; non-water heteroatoms are shown as sticks.</i>",
            styles["CaptionAudit"],
        ))
    else:
        story.append(Paragraph(
            "Structure figure not available. PyMOL rendering was skipped or failed, but the textual audit remains usable.",
            styles["BodyAudit"],
        ))
    story.append(Spacer(1, 0.28 * cm))

    # ------------------------------------------------------------------
    # Sequence and alignment evidence
    # ------------------------------------------------------------------
    figure_num = 2
    fasta_records = _discover_fasta_records(protein_dir, pdb_id)
    story.append(_section_block(
        "Sequence evidence: FASTA records",
        "Canonical, deposited, and coordinate-resolved sequences used by FRUTON before alignment and gap interpretation.",
        styles,
        usable_width,
    ))
    story.append(Spacer(1, 0.08 * cm))
    if fasta_records:
        story.extend(_build_sequence_evidence_flowables(fasta_records, styles, usable_width, protein_dir))
    else:
        story.append(Paragraph(
            "No FASTA files were found under the protein folder. This usually means the report was generated without the fasta_files.py outputs or from a reduced artifact bundle.",
            styles["BodyAudit"],
        ))
    story.append(Spacer(1, 0.22 * cm))

    alignment_assets = _discover_alignment_assets(protein_dir)
    if alignment_assets["images"] or alignment_assets["mappings"] or alignment_assets["alignments"]:
        story.append(_section_block(
            "Alignment evidence",
            "Files produced by the PDB-vs-UniProt alignment layer; figures expose matches, gaps, and unresolved regions.",
            styles,
            usable_width,
        ))
        story.append(Spacer(1, 0.08 * cm))
        alignment_rows = [
            ("alignment PNG files", ", ".join(_short_path_for_report(path, protein_dir) for path in alignment_assets["images"]) or "none"),
            ("mapping TSV files", ", ".join(_short_path_for_report(path, protein_dir) for path in alignment_assets["mappings"][:6]) or "none"),
            ("alignment FASTA files", ", ".join(_short_path_for_report(path, protein_dir) for path in alignment_assets["alignments"][:6]) or "none"),
        ]
        story.append(_build_simple_evidence_table(alignment_rows, styles, usable_width, left_header="Evidence", right_header="Files"))
        story.append(Spacer(1, 0.12 * cm))
        for alignment_png in alignment_assets["images"]:
            story.append(KeepTogether([
                Image(str(alignment_png), width=13.7 * cm, height=5.4 * cm),
                Paragraph(
                    f"<i>Figure {figure_num}. Sequence-alignment visualization generated by FRUTON. "
                    "The figure should be used to verify sequence coverage, terminal padding, internal gaps, and chain-specific PDB-to-UniProt consistency before interpreting gap-filling decisions.</i>",
                    styles["CaptionAudit"],
                ),
            ]))
            figure_num += 1
            story.append(Spacer(1, 0.14 * cm))
        story.append(Spacer(1, 0.18 * cm))

    # ------------------------------------------------------------------
    # Pipeline narrative
    # ------------------------------------------------------------------
    story.append(_section_block(
        "Preparation pipeline",
        "Step narrative follows the FRUTON runner order and keeps citations local to the action being reported.",
        styles,
        usable_width,
    ))
    story.append(Spacer(1, 0.08 * cm))

    try:
        rerun_from = int(str(pipeline_record.get("rerun.from_step", "") or "0"))
    except ValueError:
        rerun_from = 0
    rerun_ts = _record_text(pipeline_record, "rerun.timestamp")

    narrative_flowables: list[Any] = []
    tagged = _build_narratives(pdb_id, pipeline_record, pdb_citation_inline=pdb_citation_inline)
    marker_inserted = False
    if tagged:
        for step_num, para_text in tagged:
            if rerun_from > 0 and not marker_inserted and step_num >= rerun_from:
                ts_note = f" ({_safe_pdf_text(rerun_ts)})" if rerun_ts else ""
                marker_table = Table(
                    [[Paragraph(
                        f"Input modified by user — pipeline re-run from step {rerun_from}{ts_note}",
                        styles["MetricValue"],
                    )]],
                    colWidths=[usable_width - 0.9 * cm],
                )
                marker_table.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(_FRUTON_RED_SOFT)),
                    ("TEXTCOLOR", (0, 0), (-1, -1), colors.HexColor(_FRUTON_RED)),
                    ("BOX", (0, 0), (-1, -1), 0.8, colors.HexColor(_FRUTON_RED)),
                    ("LEFTPADDING", (0, 0), (-1, -1), 8),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                    ("TOPPADDING", (0, 0), (-1, -1), 5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
                ]))
                narrative_flowables.extend([marker_table, Spacer(1, 0.08 * cm)])
                marker_inserted = True
            narrative_flowables.append(Paragraph(para_text, styles["BodyAudit"]))
    else:
        narrative_flowables.append(Paragraph("No completed pipeline steps were recorded.", styles["BodyAudit"]))
    story.append(_build_narrative_table(narrative_flowables, usable_width))
    story.append(Spacer(1, 0.30 * cm))

    # ------------------------------------------------------------------
    # Protonation residue rename table
    # ------------------------------------------------------------------
    prot_in_path_str = _record_text(pipeline_record, "protonation.input_path")
    prot_out_path_str = _record_text(pipeline_record, "protonation.output_path")
    if prot_in_path_str and prot_out_path_str:
        prot_in = Path(prot_in_path_str)
        prot_out = Path(prot_out_path_str)
        residue_renames = _collect_residue_renames(prot_in, prot_out)
        if residue_renames:
            story.append(_section_block(
                "Protonation: residue assignments",
                "All residue name changes introduced by the PROPKA + gmx pdb2gmx protonation step. "
                "PROPKA-driven rows (highlighted) reflect pKa-based histidine state selection. "
                "GROMACS-driven rows reflect force-field template substitutions (e.g. disulfide "
                "bridges, selenomethionine conversion).",
                styles,
                usable_width,
            ))
            story.append(Spacer(1, 0.08 * cm))
            story.append(_build_residue_rename_table(residue_renames, styles, usable_width))
            story.append(Spacer(1, 0.22 * cm))

    # ------------------------------------------------------------------
    # Explicit audit decision and file evidence
    # ------------------------------------------------------------------
    decision, reason = _preparation_decision(pipeline_record)
    decision_rows = [
        ("preparation decision", decision),
        ("reason", reason),
        ("prepared status", _record_text(pipeline_record, "prepared_structure.status") or "not recorded"),
        ("prepared variant", _record_text(pipeline_record, "prepared_structure.variant") or "not recorded"),
        ("Ramachandran favored", (_record_text(pipeline_record, "model_eval.rama_pct_favored") + "%") if _record_text(pipeline_record, "model_eval.rama_pct_favored") else "not recorded"),
        ("Ramachandran outlier", (_record_text(pipeline_record, "model_eval.rama_pct_outlier") + "%") if _record_text(pipeline_record, "model_eval.rama_pct_outlier") else "not recorded"),
        ("clashscore", _record_text(pipeline_record, "model_eval.clashscore") or "not recorded"),
    ]
    story.append(_section_block(
        "Preparation decision",
        "Conservative audit status derived from prepared-structure, quality, metal, and non-standard-residue evidence.",
        styles,
        usable_width,
    ))
    story.append(Spacer(1, 0.08 * cm))
    story.append(_build_simple_evidence_table(decision_rows, styles, usable_width, left_header="Decision field", right_header="Recorded value"))
    story.append(Spacer(1, 0.22 * cm))

    path_rows = _collect_path_evidence(pipeline_record, protein_dir)
    if path_rows:
        story.append(_section_block(
            "Output and file evidence",
            "Important paths recorded by the flattened FRUTON state; these are the files a reviewer should open before accepting a structure.",
            styles,
            usable_width,
        ))
        story.append(Spacer(1, 0.08 * cm))
        story.append(_build_simple_evidence_table(path_rows, styles, usable_width, left_header="Record key", right_header="Path"))
        story.append(Spacer(1, 0.22 * cm))

    metal_evidence_rows = _collect_prefixed_evidence(pipeline_record, ("metals.", "metall_params."))
    if _boolish_yes(_record_text(pipeline_record, "has_metals")) or metal_evidence_rows:
        story.append(_section_block(
            "Metal and parameterization evidence",
            "Metal detection, coordination-site interpretation, and parameter-generation status are separated from generic structure preparation.",
            styles,
            usable_width,
        ))
        story.append(Spacer(1, 0.08 * cm))
        if not metal_evidence_rows:
            metal_evidence_rows = [("has_metals", _record_text(pipeline_record, "has_metals") or "yes")]
        story.append(_build_simple_evidence_table(metal_evidence_rows[:18], styles, usable_width, left_header="Record key", right_header="Value"))
        story.append(Paragraph(
            "Interpretation rule: a rendered metal pocket is evidence of detection and inspection, not proof that bonded metal parameters were accepted. MCPB.py/Gaussian acceptance should be based on explicit status fields and generated force-field files.",
            styles["SmallNote"],
        ))
        story.append(Spacer(1, 0.22 * cm))

    # ------------------------------------------------------------------
    # Metal coordination pocket figures
    # ------------------------------------------------------------------
    if metal_png_h or metal_png_dist:
        ion_label = f" ({_safe_pdf_text(ion_type)})" if ion_type else ""
        story.append(_section_block(
            f"Metal coordination pocket{ion_label}",
            "Local metal-site view for protonation-state and donor-geometry inspection.",
            styles,
            usable_width,
        ))
        story.append(Spacer(1, 0.08 * cm))

        if metal_png_h and metal_png_h.exists():
            story.append(KeepTogether([
                Image(str(metal_png_h), width=9.6 * cm, height=9.6 * cm),
                Paragraph(
                    f"<i>Figure {figure_num}. Metal coordination pocket{ion_label} with hydrogen atoms visible. "
                    "The protein backbone is a transparent grey cartoon, the metal is a FRUTON-red sphere, "
                    "coordinating residues use the IDIS-gold carbon accent, and hydrogens are white. "
                    "Hydrogen positions expose donor protonation states such as histidine δ/ε tautomer choices.</i>",
                    styles["CaptionAudit"],
                ),
            ]))
            figure_num += 1
            story.append(Spacer(1, 0.20 * cm))

        if metal_png_dist and metal_png_dist.exists():
            story.append(KeepTogether([
                Image(str(metal_png_dist), width=9.6 * cm, height=9.6 * cm),
                Paragraph(
                    f"<i>Figure {figure_num}. Metal coordination geometry{ion_label}. "
                    "Dashed FRUTON-yellow lines show distances (Å) between the metal ion and donor heavy atoms "
                    "within the direct-donor cutoff. Hydrogen atoms are hidden for clarity.</i>",
                    styles["CaptionAudit"],
                ),
            ]))
            figure_num += 1
            story.append(Spacer(1, 0.20 * cm))

        story.append(Spacer(1, 0.20 * cm))

    # ------------------------------------------------------------------
    # Ramachandran scatter plot
    # ------------------------------------------------------------------
    rama_plot_path_str = _record_text(pipeline_record, "model_eval.rama_plot_path")
    rama_plot_file = Path(rama_plot_path_str) if rama_plot_path_str else None
    if rama_plot_file and rama_plot_file.exists():
        story.append(_section_block(
            "Structure quality: Ramachandran plot",
            "Stereochemical quality evidence for the prepared or model-completed structure.",
            styles,
            usable_width,
        ))
        story.append(Spacer(1, 0.08 * cm))
        story.append(Image(str(rama_plot_file), width=9.8 * cm, height=9.8 * cm))
        pct_fav_val = _record_text(pipeline_record, "model_eval.rama_pct_favored")
        pct_out_val = _record_text(pipeline_record, "model_eval.rama_pct_outlier")
        cs_val = _record_text(pipeline_record, "model_eval.clashscore")
        quality_val = _record_text(pipeline_record, "model_eval.overall_quality")
        stats_line = ""
        if pct_fav_val:
            stats_line = (
                f" Ramachandran: {_safe_pdf_text(pct_fav_val)}% favored, "
                f"{_safe_pdf_text(pct_out_val)}% outlier. Clashscore: {_safe_pdf_text(cs_val)}. "
                f"Overall quality: {_safe_pdf_text(quality_val)}."
            )
        story.append(Paragraph(
            f"<i>Figure {figure_num}. Ramachandran φ/ψ scatter plot"
            f"{_cite('lovell2003rama', 'williams2018molprobity', 'laskowski1993procheck')}. "
            "Background: green = favored region, yellow = allowed region. "
            "Residue markers: circle = general, triangle = Gly, diamond = Pro. "
            "Dot colour: dark green = favored, orange = allowed, red = outlier."
            + stats_line + "</i>",
            styles["CaptionAudit"],
        ))
        figure_num += 1
        story.append(Spacer(1, 0.30 * cm))

    # ------------------------------------------------------------------
    # References
    # ------------------------------------------------------------------
    story.append(thin_rule)
    story.append(_section_block(
        "References",
        "Only methods triggered by the recorded pipeline state are included here.",
        styles,
        usable_width,
    ))
    seen: set[str] = set()
    for key in ref_keys:
        if key in _INLINE_REFS and key not in seen:
            seen.add(key)
            story.append(Paragraph(
                f"[{_safe_pdf_text(key)}]&nbsp;&nbsp;{_safe_pdf_text(_INLINE_REFS[key])}",
                styles["ReferenceAudit"],
            ))
    if not seen:
        story.append(Paragraph("No references collected.", styles["ReferenceAudit"]))

    chrome = _make_page_chrome(pdb_id=pdb_id, idis_logo=idis_logo, fruton_logo=fruton_logo)
    doc.build(story, onFirstPage=chrome, onLaterPages=chrome)


# ---------------------------------------------------------------------------
# BibTeX writer
# ---------------------------------------------------------------------------

def _write_bib(bib_path: Path) -> None:
    bib_path.parent.mkdir(parents=True, exist_ok=True)
    bib_path.write_text(_BIBTEX, encoding="utf-8")


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def generate_protein_report(
    protein_dir: Path,
    pdb_id: str,
    pipeline_record: dict[str, str],
    *,
    global_bib_path: Path | None = None,
) -> dict[str, Any]:
    """
    Generate a per-protein PDF report with a PyMOL figure, step narratives,
    and bibliography.

    Returns
    -------
    dict with keys: status, report_pdf_path, figure_png_path, bib_path, message.
    """
    report_dir = protein_dir / "report"
    report_dir.mkdir(parents=True, exist_ok=True)

    output_pdf = report_dir / f"{pdb_id}_report.pdf"
    output_png = report_dir / f"{pdb_id}_structure.png"

    messages: list[str] = []
    result: dict[str, Any] = {
        "pdb_id":           pdb_id,
        "status":           "failed",
        "message":          "",
        "report_pdf_path":  str(output_pdf),
        "figure_png_path":  str(output_png),
        "bib_path":         str(global_bib_path) if global_bib_path else "",
    }

    # 0. Fetch PDB crystal paper citation from RCSB (best-effort, non-blocking)
    pdb_citation_inline = ""
    try:
        pdb_cit = _fetch_pdb_citation(pdb_id)
        if pdb_cit:
            pdb_citation_inline = _format_pdb_citation_inline(pdb_cit, pdb_id)
    except Exception as exc:
        messages.append(f"PDB citation fetch: {exc}")

    # 1. Render PyMOL overview figure (best-effort)
    fig_error = _render_pymol_figure(protein_dir, pdb_id, pipeline_record, output_png)
    if fig_error:
        messages.append(f"Figure: {fig_error}")

    # 2. Render metal pocket figures when metals were parametrised (best-effort)
    metal_png_h: Path | None = None
    metal_png_dist: Path | None = None
    ion_type_str = ""
    has_metals = str(pipeline_record.get("has_metals", "")).strip().lower()
    prepared_out = str(pipeline_record.get("prepared_structure.output_path", "")).strip()
    if has_metals in ("yes", "true") and prepared_out:
        prepared_pdb = Path(prepared_out)
        if not prepared_pdb.exists():
            # Fall back: look for any PDB under prepared/
            prep_dir = protein_dir / "prepared"
            hits = sorted(prep_dir.rglob("*.pdb")) if prep_dir.exists() else []
            prepared_pdb = hits[0] if hits else prepared_pdb

        if prepared_pdb.exists():
            ion_type_str = str(pipeline_record.get("metals.ion_type", "")).strip().upper()
            if not ion_type_str:
                ion_type_str = "METAL"
            metal_out_dir = report_dir / "metal_figures"
            try:
                pocket = _render_metal_pocket_figures(
                    pdb_path=prepared_pdb,
                    ion_type=ion_type_str,
                    output_dir=metal_out_dir,
                    pdb_id=pdb_id,
                )
                metal_png_h    = pocket.get("with_h")
                metal_png_dist = pocket.get("distances")
            except Exception as exc:
                messages.append(f"Metal figures: {exc}")

    # 3. Collect citation keys for the bibliography section
    ref_keys = _collect_ref_keys(pipeline_record)

    # 4. Build the PDF
    try:
        _build_pdf(
            pdb_id=pdb_id,
            pipeline_record=pipeline_record,
            output_pdf=output_pdf,
            figure_png=output_png if output_png.exists() else None,
            ref_keys=ref_keys,
            metal_png_h=metal_png_h,
            metal_png_dist=metal_png_dist,
            ion_type=ion_type_str,
            pdb_citation_inline=pdb_citation_inline,
        )
    except Exception as exc:
        messages.append(f"PDF generation failed: {exc}")
        result["message"] = " | ".join(messages)
        return result

    # 4. Write/update shared BibTeX file
    if global_bib_path:
        try:
            _write_bib(global_bib_path)
        except Exception as exc:
            messages.append(f"BibTeX write warning: {exc}")

    result["status"] = "success"
    result["message"] = " | ".join(messages) if messages else f"Report: {output_pdf.name}"
    return result