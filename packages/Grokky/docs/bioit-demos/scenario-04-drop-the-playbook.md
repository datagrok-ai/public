# Scenario 4 — Drop the Playbook

**Theme:** Agentic orchestration (stretch / showstopper)
**Length target:** 3–5 min
**Audience lean:** Both; especially IT/platform leads who care about workflow capture
**Recording risk:** High — depends on drag-drop UX and chained-tool reliability

This is the "agent moment" — a single drag-and-drop produces a full triage dashboard.

## What the user does

1. Drags a markdown file (`cro-hit-triage.md`, content below) onto the Grokky chat panel.
2. Says: **"Run this playbook on the attached CSV."**
3. Watches.

## The playbook file (`cro-hit-triage.md`)

This file is dragged in. It IS the prompt. Authoring this file deliberately is part of the demo.

```markdown
# CRO hit-list triage playbook

You are processing a new CRO hit list. Execute the following steps in order. After each step, post a one-line status update to chat.

1. Load the attached CSV (columns: Compound_ID, SMILES, Assay_pIC50).
2. Standardize structures — canonicalize SMILES and strip salts.
3. Deduplicate against the corporate compound catalog by canonical SMILES. Tag each row as `Known` or `Novel`.
4. Flag structural alerts on all rows: PAINS, BMS, SureChEMBL. Add a `Flags` column.
5. For Novel + unflagged rows, score toxicity risks and drug-likeness (Lipinski).
6. Create a project named "CRO Triage {today}" with the current table view.
7. Share the project with the `toxicology-review` group at View access.
8. Post the project URL to chat.
```

## Setup

- The **corporate compound catalog connection** is registered and queryable by canonical SMILES. (In a customer deployment this is real; for the demo server, stage a Postgres connection with a `compound_catalog` table keyed on canonical SMILES.)
- `Chem` package installed. Verified functions used: `Chem:canonicalize` (`packages/Chem/src/package.ts:2513`), `Chem:removeWaterAndSalts` (`:3010`), `Chem:Structural Alerts` (`:1340`-area), `Chem:Toxicity Risks` (`widgets/toxicity.ts:58`), `Chem:Drug Likeness` (`widgets/drug-likeness.ts:17`).
- Structural-alerts rule sets supported (verified in `packages/Chem/src/widgets/structural-alerts.ts:6-8`): `PAINS`, `BMS`, `SureChEMBL`, `MLSMR`, `Dundee`, `Inpharmatica`, `LINT`, `Glaxo`. **Brenk is not in the built-in list** — the playbook uses SureChEMBL in its place (catches a similar class of patent-flag liabilities).
- `toxicology-review` group exists.
- The `agents/` folder file-sync is verified working (`packages/Grokky/src/package.ts:62-75`; spec in `docs/IMPLEMENTATION.md`). **Drag-drop onto chat panel is unverified UX** — fall back to paste-the-md if drag-drop isn't wired (see T3.3 in POLISH.md).
- A demo `CRO_2026Q2.csv` file is staged in the user's My Files.

## Expected sequence

1. **File arrives.** Drag-drop or paste lands the playbook content into chat. The `agents/` folder file-sync (`packages/Grokky/src/package.ts:62-75` plus `dockerfiles/claude-runtime/src/user-files.ts`) makes the file visible to the Claude runtime.
2. **Steps execute** as a mix of `datagrok-exec` blocks and MCP tool calls. Step-by-step mechanics:
   - **Standardize** → `datagrok-exec` calling `Chem:canonicalize` and `Chem:removeWaterAndSalts` (registered names; verified in `packages/Chem/src/package.ts:2513` and `:3010`).
   - **Catalog dedup** → today this is a multi-line `datagrok-exec` block that queries the catalog connection via `grok.dapi.queries` and joins on canonical SMILES. There is **no first-class dedup function**; this is the largest reliability gap in the playbook (T1.3 `joinTables({normalizeKey: 'canonicalSmiles'})` covers this once shipped).
   - **Structural alerts** → `Chem:Structural Alerts` top-menu function for PAINS / BMS / SureChEMBL filtering (`packages/Chem/src/widgets/structural-alerts.ts:6-8`).
   - **Toxicity & drug-likeness** → `Chem:Toxicity Risks` (OpenChemLib `ToxicityPredictor` — mutagenicity, tumorigenicity, irritant, reproductive effects; `widgets/toxicity.ts:58`) and `Chem:Drug Likeness` (Lipinski via `DruglikenessPredictor`; `widgets/drug-likeness.ts:17`). **Full ADMET (solubility, clearance) is not shipped** — a pre-trained ChemProp model would be needed; flagged but not in scope for this demo.
   - **Project create** → MCP `create_project` tool (`packages/Grokky/dockerfiles/mcp-server/src/index.ts:150`).
   - **Sharing** → gated `share_entity` MCP tool (T1.6). Phase-1 call returns `requires_confirmation` with the share summary; the user (Alice or presenter) confirms in chat; phase-2 executes. **Dependency on T1.6:** if T1.6's L1 isn't shipped, drop step 7 from the playbook for recording — the same trust argument as S5 applies.
3. **Status updates** appear in chat after each step: e.g. "✓ Standardized 247 structures", "✓ Found 18 catalog duplicates", "✓ 3 PAINS hits flagged", "✓ Toxicity & drug-likeness scored for 226 novel compounds", "✓ Project created", "✓ Shared with toxicology-review".
4. **Final message:** a clickable project URL.

## Wow moment

The status updates streaming in real time, ending with the URL. The audience sees the agent doing a half-day of manual work in two minutes.

## Talking points

- **IT lead:** "The playbook is a markdown file. The agent reads it as instructions and calls a registered platform function for each step."
- **Scientist:** "Standardize, dedupe, flag alerts, score toxicity and drug-likeness, build the project — six steps from the markdown, each one a platform function call."
- **Drug-discovery lead:** "This is the team's CRO-triage procedure. The agent runs it from the same markdown you'd hand a new analyst."

## Works today (verified)

- Claude Opus + MCP can sequence multi-step plans through the Claude Agent SDK runtime.
- `agents/` folder file-sync exists and is verified (`packages/Grokky/src/package.ts:62-75`; spec in `docs/IMPLEMENTATION.md`). Files in `agents/` appear in the Claude container at `/users/{userId}/agents/`.
- **Chem standardization:** `Chem:canonicalize` (`packages/Chem/src/package.ts:2513`) and `Chem:removeWaterAndSalts` (`:3010`).
- **Structural alerts:** `Chem:Structural Alerts` registered function; rule sets verified at `packages/Chem/src/widgets/structural-alerts.ts:6-8` (PAINS, BMS, SureChEMBL, MLSMR, Dundee, Inpharmatica, LINT, Glaxo).
- **Toxicity & drug-likeness panels:** `Chem:Toxicity Risks` (`widgets/toxicity.ts:58`) and `Chem:Drug Likeness` (`widgets/drug-likeness.ts:17`).
- **Project creation via MCP:** `create_project` tool (`packages/Grokky/dockerfiles/mcp-server/src/index.ts:150`).
- **Permissions API:** `grok.dapi.permissions` in `js-api/src/dapi.ts:138` (callable from `datagrok-exec` blocks).

## Needs polish (see POLISH.md) — substantial

- **Catalog dedup via `joinTables` with normalize hook** (T1.3) — today the playbook's "dedupe against catalog" step is a multi-line `datagrok-exec` block: query the catalog connection, normalize SMILES on both sides, join, tag rows. T1.3 `joinTables(df, catalogQuery, {leftKey: smilesCol, rightKey: 'canonical_smiles', normalizeKey: 'canonicalSmiles'})` collapses it to one call. Catalog dedup is not a separate skill — it's `joinTables` with the right options.
- **Drag-drop into chat UX** (T3.3) — file-sync into `agents/` is verified working; the drag-onto-chat-panel affordance is **not verified** (no DOM drop handler found in `packages/Grokky/src/ai/`). If unverified, switch to paste-the-md as the live action. Capability survives; optics change.
- **Gated `share_entity` MCP tool** (T1.6 L1) — does not exist in the MCP server today. Owned by T1.6 with the two-phase confirmation protocol. T3.4 covers the remaining MCP additions (`save_project_layout` etc.) that T1.5 needs.
- **Multi-step orchestration reliability** — six chained tool calls raise partial-failure risk. Each step should be retryable; the agent should be able to say "step 4 failed, retrying" without rerunning steps 1–3.
- **ADMET scope** — the playbook today scores only **toxicity risks + Lipinski drug-likeness** (real, OpenChemLib-based). Full ADMET (solubility, clearance) requires a pre-trained ChemProp model that is not shipped. Either keep the narrower scope (current choice) or commission a model — *out of scope for this demo*.
- **Brenk SMARTS** (data fix, optional) — to honor the original "PAINS, Brenk, BMS" framing, add Brenk SMARTS to `packages/Chem/files/alert-collection.csv` and register it in `widgets/structural-alerts.ts:6-8`. The playbook currently uses SureChEMBL instead, which is shipped.

## Backup plan

- **Tier 1 fallback (most likely):** drag-drop works but executes slowly. Pre-warm by running the playbook once cold. Recording shows real wall-clock; live demo shows pre-warmed wall-clock.
- **Tier 2 fallback:** drag-drop UI doesn't trigger reliably. The presenter pastes the playbook text into chat and says "run this." Same outcome, slightly less magic.
- **Tier 3 fallback (if multi-step orchestration breaks):** demo the first three steps live, then play a pre-recorded video of steps 4–8. Frame it as "and here's what happens next" — still effective.
- **Tier 4 (nuclear):** pull S4 from the booth rotation; lead with S1 + S3 + S5. S4 stays in the highlight reel only because that take is controlled.

## Why this scenario matters

S1, S2, S3 show *features*. S4 shows *the agent*. If you only have time to polish one for the keynote, this is the one that defines the BioIT story. But it's also the one most likely to misbehave on stage — hence the recording.
