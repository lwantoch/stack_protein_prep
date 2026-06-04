"""Step narrative builder for FRUTON per-protein PDF reports."""
from __future__ import annotations


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
    available = s("available_models")
    if prep_status in done:
        variants = [v.strip() for v in available.split("|") if v.strip()] if available else []
        if not variants:
            fallback = s("prepared_structure.variant")
            variants = [fallback] if fallback else []
        n = len(variants)
        if n == 0:
            add(12,
                "The prepared structure was assembled and written to the "
                "prepared/ directory."
            )
        elif n == 1:
            add(12,
                f"One prepared model variant was produced ({variants[0]}) "
                f"and written to the prepared/ directory."
            )
        else:
            label_list = ", ".join(variants)
            add(12,
                f"{n} prepared model variants were produced "
                f"({label_list}) and written to the prepared/ directory. "
                f"Each variant is described below."
            )
        _variant_prose = {
            "initial": (
                "The initial variant retains the deposited crystal structure "
                "without any gap-filling. Missing internal residues remain as "
                "structural chain breaks. This model is selected when the protein "
                "sequence is complete or when gap-filling is not applicable."
            ),
            "gaps": (
                "The gaps variant preserves all chain breaks present in the "
                "deposited coordinates. Each artificial terminus introduced by a "
                "missing-residue region was capped with ACE (acetyl) and NME "
                "(N-methylamide) blocking groups to prevent spurious terminal "
                "charges during simulation."
            ),
            "modeller": (
                f"The modeller variant provides a gap-free, fully connected chain: "
                f"missing residues were reconstructed by MODELLER"
                f"{_cite('sali1993modeller', 'webb2016modeller')} "
                f"(see gap-filling step above). The filled model was subsequently "
                f"protonated and assembled into the final prepared topology."
            ),
            "alphafold": (
                f"The alphafold variant provides a gap-free, fully connected chain: "
                f"missing regions were grafted from an AlphaFold2"
                f"{_cite('jumper2021alphafold')} predicted structure "
                f"(see gap-filling step above). The grafted model was subsequently "
                f"protonated and assembled into the final prepared topology."
            ),
        }
        if n > 1:
            for v in variants:
                prose = _variant_prose.get(
                    v, f"The {v} variant was assembled and written to the prepared/ directory."
                )
                add(12, prose)

        _req_range = s("range")
        _act_range = s("prepared_structure.actual_range")
        if _act_range:
            if _req_range and _act_range != _req_range:
                add(12,
                    f"The requested residue range was {_req_range}; the final "
                    f"prepared structure spans residues {_act_range} (ATOM records "
                    f"of the representative variant). The discrepancy typically "
                    f"reflects residues removed due to incomplete heavy atoms or "
                    f"re-numbering introduced by insertion-code removal."
                )
            else:
                add(12,
                    f"The prepared structure spans residues {_act_range} "
                    f"(ATOM records of the representative variant)."
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
