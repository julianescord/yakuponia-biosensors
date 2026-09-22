#!/usr/bin/env python3
"""Valida que cada constructo contenga sus partes modulares declaradas.

Comprueba, de forma reproducible y sin intervencion manual:
  1. cada parte esta presente en el constructo que la declara
  2. el orden de las partes es el esperado y no hay solapamientos
  3. cada CDS empieza en ATG, termina en stop y esta en marco
  4. el espaciado RBS->ATG cae en el rango funcional
  5. sitios de restriccion de ensamblaje (BsaI/EcoRI/XbaI/SpeI/PstI)

Salida: tabla por constructo y exit code 1 si algo falla.

Uso:  .venv/bin/python scripts/validate_assembly.py
"""
from __future__ import annotations

import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

SEQ_DIR = Path(__file__).resolve().parent.parent / "sequences" / "annotated"

# Orden esperado de partes por constructo (labels tal como los escribe reannotate.py)
EXPECTED_ORDER = {
    "insert_ftagen": ["P_pht + regulator", "RBS medium (BBa_B0032)", "GFP", "BBa_B1002"],
    "insert_metalgen": ["P_mer + regulator", "RBS strong (BBa_B0034)",
                        "arsR-like regulator", "mCherry", "BBa_B1002"],
    "insert_killswitch": [
        "J23117", "RBS medium (BBa_B0032)", "mazF", "BBa_B0015",
        "J23100", "RBS strong (BBa_B0034)", "mazE", "BBa_B0015",
    ],
}

# Enzimas relevantes para los estandares de ensamblaje en juego
SITES = {
    "BsaI": "GGTCTC",
    "EcoRI": "GAATTC",
    "XbaI": "TCTAGA",
    "SpeI": "ACTAGT",
    "PstI": "CTGCAG",
    "NotI": "GCGGCCGC",
}

# Rango de espaciado RBS -> ATG que se considera funcional
SPACER_MIN, SPACER_MAX = 5, 12


def count_sites(seq: Seq) -> dict[str, int]:
    s = str(seq).upper()
    rc = str(seq.reverse_complement()).upper()
    out = {}
    for name, site in SITES.items():
        n = s.count(site) + (rc.count(site) if site != str(Seq(site).reverse_complement()) else 0)
        if n:
            out[name] = n
    return out


def check_cds(feat, record) -> list[str]:
    sub = feat.extract(record.seq)
    errs = []
    label = feat.qualifiers.get("label", ["?"])[0]
    if len(sub) % 3:
        errs.append(f"{label}: longitud {len(sub)} no es multiplo de 3")
    if not str(sub).upper().startswith("ATG"):
        errs.append(f"{label}: no empieza en ATG (empieza en {str(sub)[:3].upper()})")
    prot = sub.translate()
    if not prot.endswith("*"):
        errs.append(f"{label}: no termina en codon stop")
    # Un stop doble/triple al final (TAATAA) es diseno deliberado, no un
    # stop interno: solo cuentan los stops antes de la cola terminal.
    body = str(prot).rstrip("*")
    internal = body.count("*")
    if internal:
        errs.append(f"{label}: {internal} codon(es) stop interno(s)")
    return errs


def main() -> int:
    failures = 0
    files = sorted(SEQ_DIR.glob("*.gb"))
    if not files:
        print(f"No hay archivos en {SEQ_DIR}. Corre antes scripts/reannotate.py",
              file=sys.stderr)
        return 1

    for path in files:
        rec = SeqIO.read(path, "genbank")
        print(f"=== {path.stem}  ({len(rec.seq)} bp, {rec.annotations['topology']}) ===")
        errs: list[str] = []

        feats = sorted((f for f in rec.features if f.type != "source"),
                       key=lambda f: int(f.location.start))
        labels = [f.qualifiers.get("label", ["?"])[0] for f in feats]

        # 1-2. orden y solapamiento
        expected = EXPECTED_ORDER.get(path.stem)
        if expected:
            if labels != expected:
                errs.append(f"orden inesperado:\n      obtenido={labels}\n      esperado={expected}")
            else:
                print(f"  orden de partes OK ({len(labels)} partes)")

        for a, b in zip(feats, feats[1:]):
            if int(a.location.end) > int(b.location.start):
                errs.append(f"solapamiento: {a.qualifiers['label'][0]} y {b.qualifiers['label'][0]}")

        # 3. integridad de CDS
        for f in feats:
            if f.type == "CDS":
                sub = f.extract(rec.seq)
                cds_errs = check_cds(f, rec)
                errs.extend(cds_errs)
                if not cds_errs:
                    aa = len(sub.translate()) - 1
                    print(f"  CDS {f.qualifiers['label'][0]:<12} {len(sub):>5} bp -> {aa:>3} aa  OK")

        # 4. espaciado RBS -> ATG
        for rbs, nxt in zip(feats, feats[1:]):
            if rbs.type == "RBS" and nxt.type == "CDS":
                # Un RBS solo alimenta un CDS de su misma hebra; si el vecino
                # inmediato es divergente, el CDS relevante es el siguiente.
                if nxt.location.strand == -1:
                    later = [f for f in feats
                             if f.type == "CDS" and f.location.strand != -1
                             and int(f.location.start) >= int(rbs.location.end)]
                    if not later:
                        continue
                    nxt = later[0]
                gap = int(nxt.location.start) - int(rbs.location.end)
                ok = SPACER_MIN <= gap <= SPACER_MAX
                tag = "OK" if ok else "FUERA DE RANGO"
                line = (f"  RBS->ATG {rbs.qualifiers['label'][0]:<24} "
                        f"{gap:>3} nt  {tag}")
                if ok:
                    print(line)
                else:
                    errs.append(line.strip()
                                + f" (esperado {SPACER_MIN}-{SPACER_MAX} nt) "
                                  "-- ver README, limitacion conocida")

        # 5. sitios de restriccion
        sites = count_sites(rec.seq)
        if sites:
            print(f"  sitios: {', '.join(f'{k}x{v}' for k, v in sorted(sites.items()))}")

        for e in errs:
            print(f"  [FALLA] {e}")
        failures += len(errs)
        print()

    print("=" * 60)
    if failures:
        print(f"{failures} problema(s) encontrados")
    else:
        print("Todos los constructos pasan la validacion")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
