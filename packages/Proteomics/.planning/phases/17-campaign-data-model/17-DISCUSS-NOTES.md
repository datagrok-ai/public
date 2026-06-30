# Phase 17 — Campaign Data Model — DISCUSS NOTES (in progress)

**Status:** `/gsd:discuss-phase 17` was started 2026-06-09 and INTERRUPTED after 3 of 4
gray areas. The session then pivoted to a long string of v1.3 viewer-UX fixes (volcano,
enrichment, publishing) — see git log `5b68427d0a..9fe7571e28`. CONTEXT.md was never
written. This file captures the decisions made so the discussion can resume.

**To resume:** re-run `/gsd:discuss-phase 17` (it will detect this dir). Confirm the
pending save-eligibility question (D3 below), then it writes 17-CONTEXT.md.

## Locked carry-forward (do NOT re-ask — from REQUIREMENTS/ROADMAP/research)
- Tag-based AppData model (NO Postgres / NO cruddy) — REQUIREMENTS "Out of Scope".
- 4 campaign tags: `proteomics.campaign`, `campaign_id`, `campaign_run_id`,
  `campaign_compound`. 2 new SEMTYPEs: `Proteomics-CompoundVid`, `Proteomics-SpcStatus`
  (register in `src/utils/proteomics-types.ts` + mirror in `detectors.js`).
- Chem via `grok.functions.call('Chem:...')` (no `@datagrok/chem` npm dep).
- Compound canonicalization = explicit `compound_id` OR Chem SMILES canonicalization
  (`Chem:convertMolNotation`); fuzzy-name matching REJECTED (CAMP-04).
- Belt-and-braces tag+column encoding for serialization survival (Phase 15 D-05).
- Audience: biologist-in-the-room demo readiness (Phase 15/16 pin).
- Phase 16 `System:AppData/Proteomics/spc/runs.csv` exists with a stable 13-col shape +
  slugify + upsert-by-(instrument,acquisition_datetime) — reusable patterns.

## DECISIONS MADE

### D1 — Storage layout + per-run persistence (CONFIRMED by user)
`campaigns/<id>/` is the CANONICAL campaign store:
- `campaigns/<id>/campaign.json` — metadata (target, owner, notes) + a run-summary list
  (one entry per run: compound, acquisition date, sample count, protein count, spc_status,
  top hit). The campaign-index DataFrame is built FROM this summary list.
- `campaigns/<id>/runs/<run_id>.csv` — each run's TRIMMED DE table, reusing Phase 15's
  `trimForPublish` allowlist (Protein ID, Gene, log2FC, p-value, adj.p-value, significant,
  direction). Phase 18 re-opens two of these to draw side-by-side volcanos + diff table.
- Phase 16's `spc/runs.csv` stays the SPC dashboard's own store, DECOUPLED. Phase 17
  OPTIONALLY back-stamps a `campaign_id` column onto matching runs.csv rows (honoring
  Phase 16's additivity promise) but does NOT depend on it.
- Tags don't travel in CSVs, so per-run identity lives in `campaign.json`'s summary entry;
  the `proteomics.campaign*` tags are belt-and-braces only.

### D2 — Compound identity / VID (CONFIRMED by user)
- Save dialog compound input auto-prefills from a detected `compound_id` or SMILES column
  (`findColumn` / `Chem:detectSmiles`), editable manually.
- `vid` = explicit `compound_id` if present, else canonical SMILES (`Chem:convertMolNotation`).
- `campaign_compound` JSON = `{vid, smiles (canonical|null), name (display), week?}`.
- **No-structure case: ALLOW name-only with EXACT match + a warning.** vid = the exact name
  string, smiles = null; index compound column shows plain text (no Molecule renderer).
  Cross-week identity then relies on typing the exact same name (exact match, never fuzzy).

### D3 — Save eligibility + SPC dependency (PENDING — user has NOT confirmed)
Claude's PROPOSAL (the AskUserQuestion for this was the one that got cut off):
- Gate ONLY on `proteomics.de_complete`. SPC status never blocks a save — the index shows
  the run's `spc_status` tag if present, else "not computed".
- Spectronaut **Candidates** runs ARE eligible (they carry the canonical DE shape) — this is
  the client-deliverable case; their index SPC column shows "n/a", sample_count shows "—".
- (Alternatives offered: require SPC too / offer inline SPC compute. Resume to decide.)

## Deferred ideas surfaced during the discussion
- True three-way top-level tabs for the analysis view (Table / Volcano / Enrichment as
  separate platform tabs, or one composite view with a `ui.tabControl`). Came up repeatedly
  during the UX-fix session; user chose the light-touch declutter for now. Candidate for its
  own small v1.4.x UX phase.
- SDF/CSV per-campaign compound upload (heavier; deferred — column-detect + manual is the
  v1.4 path).
