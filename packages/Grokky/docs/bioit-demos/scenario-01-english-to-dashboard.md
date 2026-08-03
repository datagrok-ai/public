# Scenario 1 — English → Dashboard

**Theme:** Data access democratization
**Length target:** 2–3 min
**Audience lean:** Both (IT sees curated queries + sharing; scientists see no-SQL flow)
**Recording risk:** Low

## Prompt

> Find ChEMBL compounds with binding data against thrombin where pChEMBL ≥ 7. Show them as a chemical space colored by activity, and save the result as a project called "Thrombin Hits — BioIT Demo" shared with the Demo group.

## Setup

- `Chembl` package (`@datagrok/chembl`) installed and registered — the local Postgres ChEMBL mirror with the RDKit cartridge. Bundles its own connections (`Chembl:Chembl`, `Chembl:ChemblSql`, `Chembl:Unichem`) pointing at the public Datagrok demo DB. No secrets to set.
- A curated query in `packages/Chembl/queries/` that answers "binding data for target X with pChEMBL ≥ Y" — the package already ships `activityDetailsForTarget`, `StructuresByOrganism`, etc. Confirm one matches the prompt deterministically (or add a sibling). The `Query` AI search provider routes to it — see `packages/Grokky/src/ai/search/`.
- `Chem` package installed (provides `chemicalSpace`).
- A user group named `Demo` exists; the demo user has share rights.
- Demo user is signed in to a fresh, empty workspace.

## Expected sequence

1. **Intent routing.** Grokky's `Query` provider matches "thrombin / binding / pChEMBL" to the curated Chembl query. (`CombinedAISearchAssistant` in `packages/Grokky/src/ai/search/combined-search.ts`.)
2. **Parameter inference.** Target → `"thrombin"`. Threshold → `7`. The query runs; result lands in a new view as a dataframe with SMILES + activity columns.
3. **Viewer add.** `Execute` provider emits a `datagrok-exec` block that calls `Chem:Chem Space` (registered name; impl `chemSpaceTopMenu` at `packages/Chem/src/package.ts:878`, top-menu `Chem | Analyze | Chemical Space...`). The function runs a 2D embedding of the molecule column; the resulting scatter is configured to color by the activity column.
4. **Save & share.** Agent calls MCP `create_project` (`packages/Grokky/dockerfiles/mcp-server/src/index.ts:150`), then the gated `share_entity` MCP tool (T1.6). The server returns `requires_confirmation` with the share summary; Grokky pauses; demo user confirms; share executes. A shareable URL appears in the chat.

   **Dependency on T1.6:** the sharing step assumes T1.6 L1 (gated MCP `share_entity`) is shipped. **If T1.6 isn't shipped, drop the sharing step from S1's recording** — without the gate, the "no auto-share" trust point collides with S5's claim and the rep is left explaining the inconsistency on stage.

## Wow moment

The dashboard URL pasted in chat at the end. Click it: same dashboard opens in an incognito tab as the share recipient.

## Talking points

- **IT lead:** "The agent doesn't author SQL — it routes to a query template you already approved. Permissions, audit, and parameter validation come for free. Sharing is your existing group system, no new ACL layer."
- **Scientist:** "Two sentences of English - and you get a chemistry visualization you can hand to a teammate."

## Works today

- Three AI providers (`Help`, `Execute`, `Query`) route prompts. NL→SQL matching works via Sonnet, cached daily.
- MCP `create_project` and sharing exist as tools.
- `Chem:chemicalSpace` is a real registered function.

## Needs polish (see POLISH.md)

- **Reliable viewer configuration** (T1.1) — the LLM currently sets viewer properties via raw JS; brittle for "color by X" and similar one-shot configs. Add structured `configure-viewer` skill.
- **Project creation + group sharing as a single skill** (T1.5) — currently two MCP calls the LLM must sequence; wrap as one.
- **Curated Chembl query with a pChEMBL threshold parameter** — `packages/Chembl/queries/activityDetailsForTarget.sql` returns activity details for a target but **does not filter by pChEMBL** (verified: no `pchembl_value` column in the SELECT or WHERE clause). The demo prompt's "pChEMBL ≥ 7" requires either a new sibling query with a `pChEMBLMin` parameter, or post-filtering via `Execute`. Pre-stage the parameterized version before recording.
- **Gated `share_entity` MCP tool** (T1.6 L1) — verified absent in `packages/Grokky/dockerfiles/mcp-server/src/index.ts`. T1.6 owns the gated version of this tool; once shipped, S1's sharing step works through it. Until then, sharing in S1 should be cut (see Expected sequence step 4 dependency note).

## Backup plan

- If query routing misses, the prompt has explicit terms ("thrombin", "pChEMBL ≥ 7") that the LLM can fall back to as raw SQL via `Execute`. The recording will look fine; the IT-talking-point about routing weakens.
- If `chemicalSpace` is slow on a cold start, pre-warm by running it once before the take.
- If sharing fails, demo the project URL directly — no group share, just "anyone with link."
