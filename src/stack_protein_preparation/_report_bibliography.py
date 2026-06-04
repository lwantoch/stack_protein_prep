"""BibTeX bibliography, inline reference texts, and citation utilities for FRUTON reports."""
from __future__ import annotations

import urllib.request
import json as _json

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
# BibTeX writer
# ---------------------------------------------------------------------------

def _write_bib(bib_path) -> None:
    from pathlib import Path
    bib_path = Path(bib_path)
    bib_path.parent.mkdir(parents=True, exist_ok=True)
    bib_path.write_text(_BIBTEX, encoding="utf-8")
