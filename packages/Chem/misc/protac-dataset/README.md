# PROTAC dataset

Builds `data/demo/chem/protac_degraders.csv` — 2,792 PROTAC degraders split into warhead,
linker and E3 ligand, used to demo the SAR Matrix "Use existing R-groups" mode (`Linker` is
the core, `Warhead` and `E3 Ligand` are the R-group columns). The dataset's own
documentation is `protac_degraders.README.md` next to the CSV.

```
pip install rdkit pandas openpyxl
python build.py
```

## Where the data comes from

Structures, targets and patent metadata come from PROTAC-PatentDB
([Scientific Data 12, 1840 (2025)](https://doi.org/10.1038/s41597-025-06136-9),
[figshare](https://doi.org/10.6084/m9.figshare.29351321)), released under CC BY 4.0 and so
redistributable with attribution. PROTAC-DB carries a ready-made split but forbids
redistribution, so it is deliberately not used.

## How the split is derived

The source has whole molecules only, so `split.py` derives the three parts. A PROTAC is
linear — warhead–linker–E3 ligand — and both joins are bridge bonds:

1. Match an E3 ligand scaffold (IMiD/CRBN or the VH032 family) and cut the bridge bond that
   leaves the *smallest* fragment still containing it. That converges on the same exit
   vector every time, so one E3 ligand always yields one SMILES.
2. From the exposed atom, walk outward through linker-like atoms — acyclic, or small
   saturated rings and click triazoles — and stop at the first warhead ring system.
3. Reject anything ambiguous: a branched tether, more than one exit from the linker, or a
   warhead with fewer than two rings.

Every surviving row is then reassembled with `molzip` and kept only if it reproduces the
original molecule exactly. About 48% of the source splits cleanly; the current build drops
0 rows at the round-trip check.

The VHL exit vector sits on the *tert*-leucine amine and the CRBN one on the aryl ring, so
a linker SMILES is usually specific to one E3 family. Linkers are therefore selected per
family (`CRITERIA` in `build.py`), otherwise the selection goes almost entirely CRBN.
