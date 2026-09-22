#!/usr/bin/env python3
"""Reanota los GenBank exportados de Benchling con features tipadas.

Benchling exporta casi todo como `misc_feature` sin metadatos. Este script
reconstruye la anotacion buscando cada parte modular dentro de cada
constructo y escribiendo features con tipo correcto (promoter/RBS/CDS/
terminator), etiqueta, y referencia al registro iGEM cuando aplica.

Uso:  python3 scripts/reannotate.py [--check]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

SEQ_DIR = Path(__file__).resolve().parent.parent / "sequences"
PART_DIR = SEQ_DIR / "parts"
RAW_DIR = SEQ_DIR / "raw"
OUT_DIR = SEQ_DIR / "annotated"

# Catalogo de partes: archivo fuente -> (tipo GenBank, label, nota)
PARTS = {
    "J23117": ("promoter", "J23117", "Anderson constitutive promoter, weak (iGEM BBa_J23117)"),
    "J23100": ("promoter", "J23100", "Anderson constitutive promoter, strong (iGEM BBa_J23100)"),
    "RBS_medium_B0032": ("RBS", "RBS medium (BBa_B0032)", "ribosome binding site, medium strength"),
    "RBS_strong_B0034": ("RBS", "RBS strong (BBa_B0034)", "ribosome binding site, strong"),
    "mazF_toxin": ("CDS", "mazF", "MazF endoribonuclease (toxin), E. coli MazEF type II TA system"),
    "mazE_antitoxin": ("CDS", "mazE", "MazE antitoxin, E. coli MazEF type II TA system"),
    "terminator_B0015": ("terminator", "BBa_B0015", "double terminator rrnB T1 + T7TE"),
    "terminator_B1002": ("terminator", "BBa_B1002", "rrnB-like transcriptional terminator"),
    "P_pht_regulator": ("regulatory", "P_pht + regulator", "phthalate-responsive regulatory region with cognate regulator CDS"),
    "P_mer_regulator": ("regulatory", "P_mer + regulator", "heavy-metal-responsive regulatory region with cognate regulator CDS"),
    "reporter_GFP": ("CDS", "GFP", "green fluorescent protein (Aequorea-type) -- SEE README: chromoprotein swap pending"),
    "reporter_mCherry": ("CDS", "mCherry", "red fluorescent protein -- SEE README: chromoprotein swap pending"),
    "biobrick_prefix_RFC10": ("misc_feature", "BioBrick prefix RFC[10]", "EcoRI-NotI-XbaI standard assembly prefix"),
}

# Constructos a reanotar
# Features que no son partes independientes en el export de Benchling pero
# que el analisis de marcos de lectura identifico dentro de un constructo.
# Se declaran por secuencia para que el mapeo siga siendo automatico.
EXTRA_CDS = {
    # regulador tipo ArsR/SmtB, 107 aa, codificado en la hebra reversa
    # inmediatamente 5' de la reportera de MetalGen
    "arsR_like": {
        "source": "arsR_like_regulator",
        "kind": "CDS",
        "label": "arsR-like regulator",
        "note": ("ArsR/SmtB-family metalloregulatory repressor, 107 aa, "
                 "reverse strand; HTH motif TVSHQLRFL. Identificado por "
                 "analisis de marcos de lectura, no anotado en Benchling"),
        "strand": -1,
    },
}

CONSTRUCTS = [
    "insert_ftagen", "insert_metalgen", "insert_killswitch",
    "ftagen_phthalate_biosensor", "metalgen_heavy_metal_biosensor",
    "killswitch_mazEF",
    "ftagen_sensor_module", "arsR_like_regulator",
]

# LOCUS name abreviado (<=16 chars) por constructo.
LOCUS_NAME = {
    "ftagen_phthalate_biosensor": "pYAK_FtaGen",
    "metalgen_heavy_metal_biosensor": "pYAK_MetalGen",
    "killswitch_mazEF": "pYAK_KillSw",
    "insert_ftagen": "ins_FtaGen",
    "insert_metalgen": "ins_MetalGen",
    "insert_killswitch": "ins_KillSw",
    "ftagen_sensor_module_biobrick": "FtaGen_sensor",
    "metalgen_regulator_spacer": "MetalGen_spcr",
}

NAME_MAP = {
    "ftagen_phthalate_biosensor": "ftagen_phthalate_biosensor",
    "metalgen_heavy_metal_biosensor": "metalgen_heavy_metal_biosensor",
    "killswitch_mazEF": "killswitch_mazEF",
    "insert_ftagen": "insert_ftagen",
    "insert_metalgen": "insert_metalgen",
    "insert_killswitch": "insert_killswitch",
    "ftagen_sensor_module": "ftagen_sensor_module_biobrick",
    "arsR_like_regulator": "metalgen_regulator_spacer",
}


def read_seq(path: Path) -> str:
    text = path.read_text(encoding="utf8", errors="replace")
    origin = text.split("ORIGIN")[1]
    return "".join(re.findall(r"[acgtnACGTN]", origin.replace("//", ""))).upper()


def is_circular(path: Path) -> bool:
    return "circular" in path.read_text(encoding="utf8", errors="replace").split("\n")[0]


def find_all(haystack: str, needle: str) -> list[tuple[int, int]]:
    hits, start = [], 0
    while True:
        i = haystack.find(needle, start)
        if i < 0:
            return hits
        hits.append((i + 1, i + len(needle)))
        start = i + 1


def wrap_note(text: str, width: int = 58) -> list[str]:
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width and cur:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return lines


def render_feature(kind: str, start: int, end: int, label: str, note: str,
                   strand: int = 1) -> str:
    loc = f"{start}..{end}" if strand == 1 else f"complement({start}..{end})"
    out = [f"     {kind:<16}{loc}"]
    out.append(f'                     /label="{label}"')
    note_lines = wrap_note(note)
    out.append(f'                     /note="{note_lines[0]}' + ('"' if len(note_lines) == 1 else ""))
    for extra in note_lines[1:-1]:
        out.append(f"                     {extra}")
    if len(note_lines) > 1:
        out.append(f'                     {note_lines[-1]}"')
    return "\n".join(out)


def render_origin(seq: str) -> str:
    lines = []
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60].lower()
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{i + 1:>9} {groups}")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="solo valida el ensamblaje, no escribe archivos")
    args = ap.parse_args()

    parts = {}
    for name in PARTS:
        p = PART_DIR / f"{name}.gb"
        if not p.exists():
            print(f"FALTA parte: {name}.gb", file=sys.stderr)
            return 1
        parts[name] = read_seq(p)

    if not args.check:
        OUT_DIR.mkdir(exist_ok=True)

    problems = 0
    for cname in CONSTRUCTS:
        cpath = RAW_DIR / f"{cname}.gb"
        if not cpath.exists():
            cpath = PART_DIR / f"{cname}.gb"
        if not cpath.exists():
            print(f"FALTA constructo: {cname}.gb", file=sys.stderr)
            problems += 1
            continue
        seq = read_seq(cpath)
        circular = is_circular(cpath)

        feats = []
        for pname, (kind, label, note) in PARTS.items():
            for (s, e) in find_all(seq, parts[pname]):
                feats.append((s, e, kind, label, note, 1))
        for spec in EXTRA_CDS.values():
            src = read_seq(PART_DIR / f"{spec['source']}.gb")
            for (s, e) in find_all(seq, src):
                feats.append((s, e, spec["kind"], spec["label"],
                              spec["note"], spec["strand"]))
        # BBa_B1002 es el extremo 3' de BBa_B0015: si esta contenido en otro
        # terminador ya anotado, es la misma feature vista dos veces.
        feats = [f for f in feats if not (
            f[3] == "BBa_B1002" and any(
                g[3] == "BBa_B0015" and g[0] <= f[0] and f[1] <= g[1] for g in feats)
        )]
        feats.sort(key=lambda f: (f[0], -f[1]))

        print(f"{NAME_MAP[cname]:38s} {len(seq):>6} bp  {len(feats):>2} features"
              f"  {'circular' if circular else 'linear'}")
        for s, e, kind, label, _, st in feats:
            arrow = "->" if st == 1 else "<-"
            print(f"    {s:>6}..{e:<6} {arrow} {kind:<12} {label}")

        if args.check:
            continue

        newname = NAME_MAP[cname]
        # El LOCUS name no puede exceder 16 chars sin romper las columnas
        # fijas del formato; se abrevia ahi y el nombre completo vive en
        # ACCESSION/DEFINITION.
        locus = LOCUS_NAME.get(newname, newname)[:16]
        topology = "circular" if circular else "linear  "
        header = (f"LOCUS       {locus:<16} {len(seq):>11} bp    DNA     "
                  f"{topology} SYN 22-SEP-2026\n"
                  f"DEFINITION  Yakuponia / {newname}. Reanotado desde export Benchling.\n"
                  f"ACCESSION   {newname}\n"
                  f"VERSION     {newname}\n"
                  f"KEYWORDS    .\n"
                  f"SOURCE      synthetic DNA construct\n"
                  f"  ORGANISM  synthetic DNA construct\n"
                  "FEATURES             Location/Qualifiers\n")
        body = "\n".join(render_feature(k, s, e, lb, nt, st)
                         for s, e, k, lb, nt, st in feats)
        text = header + body + "\nORIGIN\n" + render_origin(seq) + "\n//\n"
        (OUT_DIR / f"{newname}.gb").write_text(text, encoding="utf8")

    print()
    print("OK" if not problems else f"{problems} problema(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
