# BioIT Demo Scenarios — Grokky

Five demo scenarios for the BioIT conference. Designed for a mixed audience of IT/platform leads and bench/computational scientists.

## How to use these files

Each scenario is a single markdown file with this structure:

1. **Prompt** — the literal user input. Copy-paste into Grokky chat, OR drag the file into the chat window.
2. **Setup** — what must be true before recording (data loaded, connections registered, packages installed, etc.).
3. **Expected sequence** — what Grokky should do, step by step. If a step is brittle, it's flagged.
4. **Wow moment** — the single visual beat the recording should land on.
5. **Talking points** — for live voiceover or title slides/text overlay. Two flavors: IT-lead and scientist.
6. **Works today / Needs polish** — calibrates expectations. The polish items map to [`POLISH.md`](./POLISH.md).
7. **Backup plan** — if a step breaks live, what to do or fewer features delivered by Day 1.

## The five scenarios

| # | Title | Theme | Audience lean | Recording risk |
|---|-------|-------|---------------|----------------|
| 1 | [English → dashboard](./scenario-01-english-to-dashboard.md) | Data access democratization | Both | Low |
| 2 | [Tabular transform](./scenario-02-tabular-transform.md) | Working with tabular data | Scientist | Low |
| 3 | [Substructure SAR](./scenario-03-substructure-sar.md) | Tabular + chem analysis | Scientist (chem) | Medium |
| 4 | [Drop the playbook](./scenario-04-drop-the-playbook.md) | Agentic orchestration (stretch) | Both | High |
| 5 | [Governance in motion](./scenario-05-governance-in-motion.md) | Democratization + governance | IT/platform | Medium |

## Combining for the live demo

Scenarios 2 and 3 can be chained into a single ~8-min "kitchen sink" story for the booth demo: load → calculate → filter by substructure → similarity → SAR → save. The split is for recording focus; in person, run them together.

Scenarios 1 and 5 also pair well — same data democratization story, but S5 adds the governance proof point that IT buyers need before they sign.

## Sequencing tip

For the keynote/highlight reel, lead with **S1** (instant gratification) → cut to **S3** (domain depth) → close with **S5** (enterprise trust). The recordings of S2 and S4 belong in deeper-dive collateral and on the booth screen rotation.

## Related

- [`EXAMPLES.md`](../EXAMPLES.md) — ~40 single-shot prompts already documented for Grokky. Pull from this when designing new scenarios.
- [`POLISH.md`](./POLISH.md) — the missing skills/artefacts that block reliability of these scenarios. Read this before you start polish work.
