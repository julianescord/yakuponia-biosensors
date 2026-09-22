#!/usr/bin/env python3
"""Construye las versiones corregidas de los insertos (v2 y v3).

Aplica cuatro correcciones documentadas en el README y cuantificadas por el
modelo cinetico:

  1. GFP -> eforRed y mCherry -> amilCP (cromoproteinas: color a simple vista
     sin excitacion)
  2. espaciador de 7 nt en cada union RBS->ATG que estaba a 0 nt
  3. RBS propio para el reportero de MetalGen, que no tenia ninguno: el unico
     RBS del casete estaba a 324 nt, separado por el arsR reverso
  4. tag ssrA en el C-terminal de MazF, para desestabilizarla y subir la cota
     de viabilidad del kill switch

La v3 anade una quinta correccion sobre la v2:

  5. RBS con Shine-Dalgarno optimizado para B. subtilis. Los BBa_B0032 y
     B0034 estan caracterizados en E. coli y solo aparean 4 y 5 nt con el 3'
     del 16S de B. subtilis, que exige un SD mas largo y contiguo.

Cada inserto v2 se ensambla desde las partes, no se parchea sobre el v1: asi
la construccion es reproducible y el diff es legible.

Uso:  .venv/bin/python scripts/build_v2.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PARTS = ROOT / "sequences" / "parts"
OUT = {"v2": ROOT / "sequences" / "v2", "v3": ROOT / "sequences" / "v3"}


def read_seq(path: Path) -> str:
    text = path.read_text(encoding="utf8", errors="replace")
    origin = text.split("ORIGIN")[1]
    return "".join(re.findall(r"[acgtnACGTN]", origin.replace("//", ""))).upper()


def part(name: str) -> str:
    return read_seq(PARTS / f"{name}.gb")


# Cada inserto v2 como lista de (parte, tipo, etiqueta, nota, hebra).
# La secuencia sale de concatenar las partes en este orden.
LAYOUTS: dict[str, list[tuple]] = {
    "insert_ftagen_v2": [
        ("P_pht_regulator", "regulatory", "P_pht + regulator",
         "phthalate-responsive regulatory region with cognate regulator CDS", 1),
        ("RBS_medium_B0032", "RBS", "RBS medium (BBa_B0032)",
         "ribosome binding site, medium strength", 1),
        ("spacer_RBS_7nt", "misc_feature", "spacer 7 nt",
         "restores functional RBS-ATG distance; was 0 nt in v1", 1),
        ("chromoprotein_eforRed", "CDS", "eforRed",
         "red chromoprotein (BBa_K592012); replaces GFP, visible without "
         "excitation", 1),
        ("terminator_B1002", "terminator", "BBa_B1002",
         "rrnB-like transcriptional terminator", 1),
    ],
    "insert_metalgen_v2": [
        ("P_mer_regulator", "regulatory", "P_mer + regulator",
         "heavy-metal-responsive regulatory region with cognate regulator CDS", 1),
        ("RBS_strong_B0034", "RBS", "RBS strong (BBa_B0034)",
         "ribosome binding site for the arsR-like regulator", 1),
        ("arsR_like_regulator", "CDS", "arsR-like regulator",
         "ArsR/SmtB-family metalloregulatory repressor, 107 aa, reverse "
         "strand; HTH motif TVSHQLRFL", -1),
        ("RBS_strong_B0034", "RBS", "RBS strong (BBa_B0034)",
         "ribosome binding site for the reporter; absent in v1", 1),
        ("spacer_RBS_7nt", "misc_feature", "spacer 7 nt",
         "functional RBS-ATG distance for the reporter", 1),
        ("chromoprotein_amilCP", "CDS", "amilCP",
         "blue-violet chromoprotein (BBa_K592009); replaces mCherry, visible "
         "without excitation", 1),
        ("terminator_B1002", "terminator", "BBa_B1002",
         "rrnB-like transcriptional terminator", 1),
    ],
    "insert_killswitch_v2": [
        ("J23117", "promoter", "J23117",
         "Anderson constitutive promoter, weak (iGEM BBa_J23117)", 1),
        ("RBS_medium_B0032", "RBS", "RBS medium (BBa_B0032)",
         "ribosome binding site, medium strength", 1),
        ("spacer_RBS_7nt", "misc_feature", "spacer 7 nt",
         "restores functional RBS-ATG distance; was 0 nt in v1", 1),
        ("mazF_toxin_ssrA", "CDS", "mazF-ssrA",
         "MazF endoribonuclease with C-terminal AANDENYALAA degradation tag; "
         "destabilised to raise the viability bound", 1),
        ("terminator_B0015", "terminator", "BBa_B0015",
         "double terminator rrnB T1 + T7TE", 1),
        ("J23100", "promoter", "J23100",
         "Anderson constitutive promoter, strong (iGEM BBa_J23100)", 1),
        ("RBS_strong_B0034", "RBS", "RBS strong (BBa_B0034)",
         "ribosome binding site, strong", 1),
        ("spacer_RBS_7nt", "misc_feature", "spacer 7 nt",
         "restores functional RBS-ATG distance; was 0 nt in v1", 1),
        ("mazE_antitoxin", "CDS", "mazE",
         "MazE antitoxin, E. coli MazEF type II TA system", 1),
        ("terminator_B0015", "terminator", "BBa_B0015",
         "double terminator rrnB T1 + T7TE", 1),
    ],
}

LOCUS_NAME = {
    "insert_ftagen_v2": "ins_FtaGen_v2",
    "insert_metalgen_v2": "ins_MetalGn_v2",
    "insert_killswitch_v2": "ins_KillSw_v2",
    "insert_ftagen_v3": "ins_FtaGen_v3",
    "insert_metalgen_v3": "ins_MetalGn_v3",
    "insert_killswitch_v3": "ins_KillSw_v3",
}


def build_mazF_ssrA() -> None:
    """Fusiona el tag ssrA al C-terminal de MazF, quitando su stop."""
    mazf = part("mazF_toxin")
    tag = part("ssrA_tag_MazF")
    # MazF termina en TAATAA (stop doble): se retira para fusionar en marco.
    while mazf.endswith("TAA") or mazf.endswith("TGA") or mazf.endswith("TAG"):
        mazf = mazf[:-3]
    fused = mazf + tag
    assert len(fused) % 3 == 0, "la fusion rompe el marco de lectura"

    body = "\n".join(
        f"{i + 1:>9} " + " ".join(fused[i:i + 60].lower()[j:j + 10]
                                  for j in range(0, len(fused[i:i + 60]), 10))
        for i in range(0, len(fused), 60)
    )
    (PARTS / "mazF_toxin_ssrA.gb").write_text(
        f"LOCUS       {'mazF_ssrA':<16} {len(fused):>11} bp    DNA     "
        "linear   SYN 22-SEP-2026\n"
        "DEFINITION  Yakuponia / mazF_toxin_ssrA.\n"
        "ACCESSION   mazF_toxin_ssrA\nVERSION     mazF_toxin_ssrA\n"
        "KEYWORDS    .\nSOURCE      synthetic DNA construct\n"
        "  ORGANISM  synthetic DNA construct\n"
        "FEATURES             Location/Qualifiers\n"
        f"     CDS             1..{len(fused)}\n"
        '                     /label="mazF-ssrA"\n'
        '                     /note="MazF fused in frame to the AANDENYALAA\n'
        '                     ssrA degradation tag"\n'
        "ORIGIN\n" + body + "\n//\n", encoding="utf8")
    print(f"  mazF_toxin_ssrA.gb            {len(fused):>5} bp "
          f"({len(fused) // 3 - 1} aa)")


def render(name: str, layout: list[tuple], ver: str = "v2") -> str:
    seq, feats, pos = "", [], 0
    for pname, kind, label, note, strand in layout:
        # Las partes se almacenan ya en la orientacion del inserto, tal como
        # las extrajo reannotate.py. `strand` solo indica como anotarlas: no
        # hay que complementar la secuencia o se invertiria dos veces.
        s = part(pname)
        feats.append((pos + 1, pos + len(s), kind, label, note, strand))
        seq += s
        pos += len(s)

    lines = []
    for start, end, kind, label, note, strand in feats:
        loc = f"{start}..{end}" if strand == 1 else f"complement({start}..{end})"
        lines.append(f"     {kind:<16}{loc}")
        lines.append(f'                     /label="{label}"')
        words, wrapped, cur = note.split(), [], ""
        for w in words:
            if len(cur) + len(w) + 1 > 58 and cur:
                wrapped.append(cur)
                cur = w
            else:
                cur = f"{cur} {w}".strip()
        wrapped.append(cur)
        if len(wrapped) == 1:
            lines.append(f'                     /note="{wrapped[0]}"')
        else:
            lines.append(f'                     /note="{wrapped[0]}')
            for x in wrapped[1:-1]:
                lines.append(f"                     {x}")
            lines.append(f'                     {wrapped[-1]}"')

    origin = "\n".join(
        f"{i + 1:>9} " + " ".join(seq[i:i + 60].lower()[j:j + 10]
                                  for j in range(0, len(seq[i:i + 60]), 10))
        for i in range(0, len(seq), 60)
    )
    extra = (", mas RBS optimizado para B. subtilis."
             if ver == "v3" else ".")
    return (
        f"LOCUS       {LOCUS_NAME[name]:<16} {len(seq):>11} bp    DNA     "
        "linear   SYN 22-SEP-2026\n"
        f"DEFINITION  Yakuponia / {name}. Correcciones: cromoproteinas, "
        f"espaciadores RBS y tag ssrA{extra}\n"
        f"ACCESSION   {name}\nVERSION     {name}\nKEYWORDS    .\n"
        "SOURCE      synthetic DNA construct\n"
        "  ORGANISM  synthetic DNA construct\n"
        "FEATURES             Location/Qualifiers\n"
        + "\n".join(lines) + "\nORIGIN\n" + origin + "\n//\n"
    )


# La v3 es la v2 con los RBS sustituidos por versiones aptas para B. subtilis.
RBS_V3 = {
    "RBS_medium_B0032": ("RBS_Bsub_medium", "RBS medium (B. subtilis)",
                         "SD AAAGGAGG, 6 nt de apareamiento con el 16S de "
                         "B. subtilis; BBa_B0032 solo apareaba 4 nt"),
    "RBS_strong_B0034": ("RBS_Bsub_strong", "RBS strong (B. subtilis)",
                         "SD AAAGGAGGTG, 7 nt de apareamiento; BBa_B0034 solo "
                         "apareaba 5 nt"),
}


def to_v3(layout: list[tuple]) -> list[tuple]:
    """Devuelve el layout de la v2 con los RBS sustituidos."""
    out = []
    for pname, kind, label, note, strand in layout:
        if pname in RBS_V3:
            new_part, new_label, new_note = RBS_V3[pname]
            out.append((new_part, kind, new_label, new_note, strand))
        else:
            out.append((pname, kind, label, note, strand))
    return out


def main() -> int:
    print("Construyendo partes derivadas")
    build_mazF_ssrA()

    for ver in ("v2", "v3"):
        OUT[ver].mkdir(exist_ok=True)
        print()
        print(f"Construyendo insertos {ver}")
        for name, layout in LAYOUTS.items():
            vname = name.replace("_v2", f"_{ver}")
            lay = layout if ver == "v2" else to_v3(layout)
            text = render(vname, lay, ver)
            (OUT[ver] / f"{vname}.gb").write_text(text, encoding="utf8")
            n = len("".join(re.findall(r"[acgt]", text.split("ORIGIN")[1])))
            print(f"  {vname:<24} {n:>6} bp  {len(lay)} partes")

    print()
    print("Escritos en sequences/v2/ y sequences/v3/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
