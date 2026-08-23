#!/usr/bin/env python3
"""Aggregate per-PDB JSONs from fruton_generic_bench into audit CSV + summary.

Usage:
    python fruton_generic_aggregate.py <BENCH_OUT_DIR>

Reads every ``<PDB>.json`` in the directory, extracts the
ProteinDeliveryReport, emits the reviewer-facing audit CSV +
markdown summary.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, '/mnt/netapp2/Store_uni/home/otras/hcx/lwa/repos/FRUTON/stack_protein_prep/src')

from stack_protein_preparation._component_confidence import (
    Confidence,
    ComponentConfidence,
    ProteinDeliveryReport,
    summarise_delivery,
)


def _reload_delivery(payload: dict) -> ProteinDeliveryReport:
    """Rehydrate a ProteinDeliveryReport from its .to_dict() JSON."""
    d = payload.get("delivery") or {}
    components = []
    for c in d.get("components", []):
        try:
            components.append(ComponentConfidence(
                component_type=c["component_type"],
                name=c["name"],
                confidence=Confidence(c["confidence"]),
                reason=c["reason"],
                suggested_action=c.get("suggested_action", ""),
                method=c.get("method", ""),
                details=dict(c.get("details", {})),
            ))
        except Exception:
            continue
    return ProteinDeliveryReport(
        pdb=d.get("pdb", payload.get("pdb", "?")),
        components=components,
        model_written=bool(d.get("model_written", True)),
        tleap_loads=d.get("tleap_loads"),
        md_deck_written=bool(d.get("md_deck_written", False)),
        notes=d.get("notes", ""),
    )


def _row_to_csv(payload: dict, delivery: ProteinDeliveryReport) -> dict:
    ct = delivery.component_type_counts()
    def _sum(bucket, level):
        return sum(v.get(level.value, 0) for v in ct.values())
    return {
        "pdb": delivery.pdb,
        "overall_status": delivery.overall_status,
        "overall_confidence": delivery.overall_confidence.value,
        "model_written": delivery.model_written,
        "n_gaps": payload.get("n_gaps", ""),
        "delta_n": payload.get("delta_n", ""),
        "rolled_gaps": payload.get("rolled_gaps", ""),
        "clash_pairs": payload.get("n_clash_pairs", ""),
        "omega_np": payload.get("n_omega_non_planar", ""),
        "omega_cis_nonpro": payload.get("n_omega_cis_nonpro", ""),
        "rama_outliers": payload.get("n_rama_outlier", ""),
        "peptide_bonds_broken": payload.get("n_peptide_bonds_broken", ""),
        "gate_pass": payload.get("gate_pass", ""),
        "refine_seconds": payload.get("refine_seconds", ""),
        "n_high": _sum(ct, Confidence.HIGH),
        "n_medium": _sum(ct, Confidence.MEDIUM),
        "n_low": _sum(ct, Confidence.LOW),
        "n_failed": _sum(ct, Confidence.FAILED),
        "component_types": ";".join(sorted(ct.keys())),
        "action_items": " | ".join(delivery.action_items()),
        "notes": delivery.notes,
        "error": payload.get("error", ""),
    }


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: fruton_generic_aggregate.py <BENCH_OUT_DIR>", file=sys.stderr)
        return 2
    out_dir = Path(argv[1])
    if not out_dir.is_dir():
        print(f"[aggregate] not a directory: {out_dir}", file=sys.stderr)
        return 2

    per_pdb = sorted(
        p for p in out_dir.glob("*.json")
        if p.name not in {"combined_results.json"}
        and not p.name.startswith("audit_")
    )
    if not per_pdb:
        print(f"[aggregate] no per-PDB JSONs in {out_dir}", file=sys.stderr)
        return 2

    payloads = []
    deliveries: list[ProteinDeliveryReport] = []
    for p in per_pdb:
        try:
            payload = json.loads(p.read_text())
        except Exception as e:
            print(f"[aggregate] SKIP {p.name}: {e}", file=sys.stderr)
            continue
        payloads.append(payload)
        deliveries.append(_reload_delivery(payload))

    # audit_report.csv
    csv_path = out_dir / "audit_report.csv"
    fieldnames = list(_row_to_csv(payloads[0], deliveries[0]).keys())
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for pl, d in zip(payloads, deliveries):
            w.writerow(_row_to_csv(pl, d))

    # audit_summary.md
    s = summarise_delivery(deliveries)
    lines = [
        "# FRUTON audit summary",
        "",
        f"**{s['n_proteins']} proteins** processed",
        "",
        "| Overall status | n | % |",
        "|---|---:|---:|",
    ]
    for status in ("delivered_full_confidence", "delivered_with_notes",
                   "delivered_needs_review", "not_delivered"):
        n = s["by_overall_status"].get(status, 0)
        pct = s["by_overall_status_pct"].get(status, 0.0)
        lines.append(f"| {status} | {n} | {pct:.1f} |")
    lines.append("")
    lines.append(f"**Total components:** {s['component_totals']['total']} — "
                 f"HIGH={s['component_totals']['high']}, "
                 f"MEDIUM={s['component_totals']['medium']}, "
                 f"LOW={s['component_totals']['low']}, "
                 f"FAILED={s['component_totals']['failed']}")
    lines.append("")

    # List proteins needing review with reasons
    needs = [d for d in deliveries if d.overall_status == "delivered_needs_review"]
    if needs:
        lines.append(f"## delivered_needs_review ({len(needs)})")
        lines.append("")
        for d in needs:
            reason_line = "; ".join(d.action_items()[:2]) or d.notes
            lines.append(f"- **{d.pdb}** — {reason_line}")
        lines.append("")

    not_delivered = [d for d in deliveries if d.overall_status == "not_delivered"]
    if not_delivered:
        lines.append(f"## not_delivered ({len(not_delivered)})")
        lines.append("")
        for d in not_delivered:
            lines.append(f"- **{d.pdb}** — {d.notes}")
        lines.append("")

    (out_dir / "audit_summary.md").write_text("\n".join(lines) + "\n")

    print(f"[aggregate] {s['n_proteins']} proteins")
    for k, v in s["by_overall_status"].items():
        print(f"[aggregate]   {k:<30} {v}")
    print(f"[aggregate] -> {csv_path}")
    print(f"[aggregate] -> {out_dir / 'audit_summary.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
