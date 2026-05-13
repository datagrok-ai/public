---
name: chem-intro
description: Orient cheminformatics work in Datagrok — pick the right notation, validate it, and choose between `grok.chem.*` and the Bio package
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# chem-intro

## When to use

You're starting cheminformatics work in a Datagrok package or Console
session and need to confirm which molecular notation your data uses,
whether it's tractable as a small molecule, and which API surface to
call. Triggers: "is this SMILES or molblock?", "convert SDF to MOL",
"validate SMARTS before substructure search", "do I use Chem or Bio
for this peptide?".

## Prerequisites

- A logged-in Datagrok session (`Tools → Console`) or a package with
  `import * as grok from 'datagrok-api/grok';` and
  `import * as DG from 'datagrok-api/dg';` (knowledge `DG-FACT-224`).
- The `Chem` package installed — `chem.checkSmiles`/`isSmarts`/`convert`
  delegate to `Chem:isSmiles` / `Chem:isSmarts` / `Chem:convertMolNotation`
  via `Func.find` and throw "function not found" without it
  (knowledge `DG-FACT-253`).
- A column or sample string of molecules. Datagrok semantically types
  small molecules as `SEMTYPE.MOLECULE`; biopolymers (proteins,
  oligonucleotides, peptide chains) as `SEMTYPE.MACROMOLECULE` — the
  latter are handled by the `Bio` package, NOT `grok.chem.*`
  (knowledge `DG-FACT-252`).

## Steps

1. **Identify the notation.**
   ```typescript
   const s = '...';
   if (grok.chem.isMolBlock(s))      console.log('molblock');
   else if (grok.chem.isSmarts(s))   console.log('smarts');
   else if (grok.chem.checkSmiles(s)) console.log('smiles');
   else                              console.log('unknown');
   ```
   Expected: exactly one branch fires for a well-formed string.
   `isMolBlock` is pure JS (marker check `s.includes('M  END')`,
   `js-api/src/chem.ts:68`); `checkSmiles` and `isSmarts` are sync but
   route through the Chem package (knowledge `DG-FACT-253`). The full
   notation enum is `DG.chem.Notation.{Smiles, CxSmiles, Smarts,
   CxSmarts, MolBlock, V3KMolBlock, Unknown}` — seven members
   (knowledge `DG-FACT-251`).

2. **Pick the storage format from the use case.**
   The article documents five purposes; the table maps each to a
   `DG.chem.Notation` (knowledge `DG-FACT-251`). Note that SDF is a
   *file*, not a notation — it's concatenated `M  END\n$$$$\n` records
   loaded by `grok.data.loadTable`, one row per molecule.

   | Need | Format | `chem.Notation` |
   |---|---|---|
   | Compact column storage (single line) | SMILES | `Smiles` |
   | 2D coordinates / atom-level annotations | MOL V2000 | `MolBlock` |
   | >999 atoms, enhanced stereo, R-groups | MOL V3000 | `V3KMolBlock` |
   | Multi-record file (compounds + activity) | SDF (file) | (load as table) |
   | Substructure pattern with logical operators | SMARTS | `Smarts` |
   | Chemical-reaction transformation | SMIRKS | (use as SMARTS) |

3. **Convert between formats when needed.**
   ```typescript
   const mol = grok.chem.convert(
     'CCO', DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
   const smi = grok.chem.convert(
     mol, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles);
   ```
   Expected: `mol` ends with `M  END`; `smi === 'CCO'` modulo
   canonicalisation. `convert` is sync (calls `Chem:convertMolNotation`
   via `Func.find`); use `Notation.Unknown` as source when the input
   may be SMILES or molblock and let RDKit detect (knowledge
   `DG-FACT-253`). SMARTS → SMILES is meaningless — patterns are not
   real molecules; the wrapper warns via
   `chem.smilesFromSmartsWarning()` and returns `''`.

4. **Gate long strings before RDKit calls.**
   ```typescript
   const MAX_SMILES_LENGTH = 5000;
   if (!grok.chem.isMolBlock(molecule) &&
       molecule?.length > MAX_SMILES_LENGTH)
     return;
   ```
   Expected: structures whose SMILES exceeds 5000 characters are
   rejected. This is the convention used by `Chem/src/widgets/{scaffold-tree,
   pharmacophore-features,structural-alerts}.ts` and SureChembl's
   `patentDataSearch` to keep RDKit from choking on macromolecules
   smuggled in as SMILES (knowledge `DG-FACT-254`). Molblocks bypass
   the length check because their atom count is encoded in the
   counts-line header, not in line length.

5. **Choose the API namespace (Chem vs Bio).**
   ```typescript
   import * as DG from 'datagrok-api/dg';
   const semType = column.semType;
   if (semType === DG.SEMTYPE.MOLECULE)            // grok.chem.*
     console.log('chem');
   else if (semType === DG.SEMTYPE.MACROMOLECULE)  // Bio package
     console.log('bio');
   else if (semType === DG.SEMTYPE.MOLECULE3D)     // 3D poses → Bio/Docking
     console.log('molecule3D');
   ```
   Expected: cheminformatics calls (`searchSubstructure`, `findSimilar`,
   `descriptors`, `mcs`, `rGroup`) only operate on `Molecule`-typed
   columns. Proteins, oligonucleotides and peptide chains carry
   `Macromolecule` (knowledge `DG-FACT-252`); use the `Bio` package
   (e.g. `to-atomic-level-widget`) to expand into a `Molecule` column
   before any `grok.chem` call. If `column.semType` is empty,
   run `await grok.data.detectSemanticTypes(t)` first (knowledge
   `DG-FACT-249`).

6. **Verify the toolkit is reachable.**
   ```typescript
   grok.shell.info(
     `Notations: ${Object.values(DG.chem.Notation).length}, ` +
     `RDKit ok: ${grok.chem.checkSmiles('CCO')}`);
   ```
   Expected: a balloon reading `Notations: 7, RDKit ok: true`. Anything
   else (`function not found`, `false` for canonical ethanol) means
   the Chem package isn't installed or the bundle didn't load.

## Common failure modes

- **`grok.chem.checkSmiles is not a function` / `Function 'isSmiles'
  not found`.** The Chem package isn't installed in this environment.
  `chem.checkSmiles`/`isSmarts`/`convert` are JS-API wrappers that
  `Func.find` against the Chem package (knowledge `DG-FACT-253`). Fix:
  install Chem (`grok install Chem`) or restrict yourself to
  `chem.isMolBlock` for marker-based detection.
- **SDF loaded as a single SMILES string.** SDF is a *file format*
  (concatenated `M  END\n$$$$\n` records), not a single notation —
  there is no `DG.chem.Notation.Sdf` (knowledge `DG-FACT-251`). Fix:
  `grok.data.loadTable(path)` splits records into one row per
  molecule; consume the per-row `MolBlock` afterward.
- **`grok.chem.convert(smarts, Smarts, Smiles)` returns `''` and a
  yellow balloon.** SMARTS encodes atom/bond *queries* (e.g. "C or O",
  "any aromatic bond"), not real molecules — there's no canonical
  SMILES for a pattern. Fix: keep SMARTS as SMARTS for
  `searchSubstructure`'s `molBlockFailover` and `mcs(returnSmarts=true)`
  outputs; render via `grok.chem.svgMol` if a picture is needed.
- **Substructure search of a peptide column finds nothing.** Peptides
  ≥ ~50 residues carry `SEMTYPE.MACROMOLECULE`; `grok.chem.*` is
  silently no-op on non-`Molecule` columns (knowledge `DG-FACT-252`).
  Fix: route through the `Bio` package or expand to an atomic-level
  `Molecule` column first.
- **RDKit "Cannot kekulize molecule" / hangs on a long string.**
  Macromolecules (≥5000 char SMILES) smuggle in past UI guards. Fix:
  apply the 5000-char gate from step 4 before any
  `searchSubstructure` / `findSimilar` / `mcs` call (knowledge
  `DG-FACT-254`).

## See also

- Source articles:
  - `help/develop/domains/chem/intro.md` (this overview).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-235` (chem namespace), `DG-FACT-236` (5 RDKit notations),
  `DG-FACT-249` (detectSemanticTypes), `DG-FACT-251` (format taxonomy
  → `Notation` enum), `DG-FACT-252` (Molecule vs Macromolecule scope),
  `DG-FACT-253` (format-detection trio + `convert`), `DG-FACT-254`
  (5000-char SMILES guard).
- Reference packages:
  - `packages/Chem/src/widgets/{scaffold-tree,pharmacophore-features,
    structural-alerts}.ts` — 5000-char gate.
  - `packages/SureChembl/src/package.ts:193` — same gate, named
    constant `MAX_SMILES_LENGTH`.
  - `packages/Chembl/src/package.ts:137` — `Notation.Unknown → Smiles`
    canonicalisation.
  - `packages/Bio/src/utils/monomer-lib/monomer-manager/monomer-manager.ts`
    — validate-then-convert flow against `DG.chem.Notation`.
- Related skills:
  - `cheminformatics` — actual API calls (`searchSubstructure`,
    `findSimilar`, `mcs`, `descriptors`, `svgMol`/`canvasMol`).
  - `docking` — AutoDock pipeline (3D pose generation, `Molecule3D`).
  - `hit-triage-system` — annotate + curate compound subsets.
  - `register-identifiers` — register custom semantic-type renderers.
