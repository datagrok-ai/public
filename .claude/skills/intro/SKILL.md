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

1. **Identify the notation.** `isMolBlock` is pure JS; `checkSmiles`
   and `isSmarts` route through the Chem package (`DG-FACT-253`).
   The notation enum has 7 members (`DG-FACT-251`).
   ```typescript
   const s = '...';
   if (grok.chem.isMolBlock(s))      console.log('molblock');
   else if (grok.chem.isSmarts(s))   console.log('smarts');
   else if (grok.chem.checkSmiles(s)) console.log('smiles');
   else                              console.log('unknown');
   ```

2. **Pick the storage format from the use case** (`DG-FACT-251`). SDF
   is a *file*, not a notation — load it as a table.

   | Need | Format | `chem.Notation` |
   |---|---|---|
   | Compact column storage (single line) | SMILES | `Smiles` |
   | 2D coordinates / atom-level annotations | MOL V2000 | `MolBlock` |
   | >999 atoms, enhanced stereo, R-groups | MOL V3000 | `V3KMolBlock` |
   | Multi-record file (compounds + activity) | SDF (file) | (load as table) |
   | Substructure pattern with logical operators | SMARTS | `Smarts` |
   | Chemical-reaction transformation | SMIRKS | (use as SMARTS) |

3. **Convert between formats when needed.** SMARTS → SMILES returns
   `''` — patterns aren't real molecules (`DG-FACT-253`).
   ```typescript
   const mol = grok.chem.convert(
     'CCO', DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
   const smi = grok.chem.convert(
     mol, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles);
   ```

4. **Gate long strings before RDKit calls** (`DG-FACT-254`). SMILES
   over 5000 chars typically smuggle in macromolecules.
   ```typescript
   const MAX_SMILES_LENGTH = 5000;
   if (!grok.chem.isMolBlock(molecule) &&
       molecule?.length > MAX_SMILES_LENGTH)
     return;
   ```

5. **Choose the API namespace (Chem vs Bio)** by `semType`
   (`DG-FACT-252`). `grok.chem.*` operates only on `MOLECULE`. Use
   `Bio` for `MACROMOLECULE`; for empty `semType`, run
   `await grok.data.detectSemanticTypes(t)` first (`DG-FACT-249`).
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

6. **Verify the toolkit is reachable.**
   ```typescript
   grok.shell.info(
     `Notations: ${Object.values(DG.chem.Notation).length}, ` +
     `RDKit ok: ${grok.chem.checkSmiles('CCO')}`);
   ```

## Common failure modes

- **`checkSmiles is not a function` / `Function 'isSmiles' not found`.**
  Chem package not installed (`DG-FACT-253`).
- **SDF loaded as a single SMILES.** Load via
  `grok.data.loadTable(path)` (`DG-FACT-251`).
- **`convert(smarts, Smarts, Smiles)` returns `''`.** SMARTS are
  queries, not molecules. Keep them as SMARTS (`DG-FACT-253`).
- **Substructure search on peptides finds nothing.** Peptides carry
  `MACROMOLECULE` — use `Bio` (`DG-FACT-252`).
- **RDKit "Cannot kekulize" / hangs.** Apply the 5000-char gate
  before calling (`DG-FACT-254`).

## See also

- Source: `help/develop/domains/chem/intro.md`.
- Knowledge: `DG-FACT-235`, `236`, `249`, `251`–`254`.
- Related skills: `cheminformatics`, `docking`, `register-identifiers`.
