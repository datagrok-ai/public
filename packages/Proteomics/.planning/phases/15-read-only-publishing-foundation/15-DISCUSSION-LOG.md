# Phase 15: Read-Only Publishing Foundation - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-07
**Phase:** 15-read-only-publishing-foundation
**Areas discussed:** Target identifier shape, Reviewer-group sourcing & UX, Publish destination + inheritance enforcement, Republish UX flow + enrichment carry shape

---

## Pre-discussion: Areas auto-excluded as already-decided

The following gray areas were NOT presented because REQUIREMENTS.md + research had already locked them. Listed here so the planner / researcher can see what was treated as fixed input:

- Trim allowlist (PUB-02): Protein ID, Gene, log2FC, p-value, adj.p-value, sig, direction — locked
- Deep-clone-first publish primitive (Pitfall 1) — locked
- `assertPublishedShape` round-trip test as phase gate (Pitfall 3, success criterion 3) — locked
- Belt-and-braces tag-AND-column metadata encoding (PUB-11) — locked
- Post-grant `permissions.get` verification (Pitfall 2) — locked
- Republish creates new Project + `proteomics.superseded_by` on prior (PUB-10) — locked
- Project name pattern `Proteomics-Review-<target>-v<N>-<date>` (PUB-09) — locked
- Audit context panel reads `proteomics.*` tags (PUB-07) — locked
- Zero new npm dependencies (v1.4 milestone constraint) — locked

---

## Target identifier shape

| Option | Description | Selected |
|--------|-------------|----------|
| Freeform string + slug sanitization | User types any string; helper slugifies for project name, raw label in tag + metadata column. | ✓ |
| Freeform + per-project history dropdown | Free typing + "Recent targets" computed from prior shares to reduce drift. | |
| Constrained to a Datagrok-platform target entity | First-class Target entity from a registry. Would pull v1.5 scope into v1.4. | |
| Skip — use Claude's default | (Would default to option 1.) | |

**User's choice:** Freeform string + slug sanitization (option 1, recommended).
**Notes:** Aligns with REQUIREMENTS.md "Future Requirements" deferral of first-class target taxonomy to v1.5+. Zero new infra; raw label preserved for human readability while the slug keeps project names well-formed.

---

## Reviewer-group sourcing & UX

| Option | Description | Selected |
|--------|-------------|----------|
| `DG.ChoiceInput` from `dapi.groups.list` | Standard platform widget populated from real groups. | ✓ |
| Freeform name + lookup-or-create hint | User types a group name; helper looks up. Invites typos. | |
| ChoiceInput + "recent targets → group" memory | Same as option 1 + pre-selects per-target last group. | |
| Skip — use Claude's default | (Would default to option 1.) | |

**User's choice:** ChoiceInput from `dapi.groups.list` (option 1, recommended).
**Notes:** Matches platform muscle memory; no typo path to non-existent groups; membership management stays upstream of this package per REQUIREMENTS.md "Out of Scope".

---

## Publish destination + inheritance enforcement

This area required a clarification round. The user asked, in essence, *"what do we mean by 'proteomics reviewer'?"* — a fellow practitioner or a consuming biologist? Claude pulled evidence from PROJECT.md ("biologist consumers"), ROADMAP.md, PITFALLS.md Pitfall 14 ("Cytokinetics demo audience contains biologists"), and the source todo `2026-05-11-share-analysis-read-only-with-biologics-team-filed-by-target.md` to confirm **reviewer = consuming biologist** (the primary persona). The destination question was then re-issued with that framing baked in.

| Option | Description | Selected |
|--------|-------------|----------|
| Publishing user's personal space + verify-and-rollback | Project owned by expert; biologist gets a link. | |
| Single shared 'Proteomics-Reviews' Space + verify-and-rollback | All publishes in one Space; admin audit aggregation. | |
| Per-target Space (Proteomics-Review-<slug>) + verify-and-rollback | Auto-create Space per target; strongest target-isolation. | |
| Skip — use Claude's default | (Would default to option 1.) | |

**User's choice:** Custom — "single general space with sub-space per target" (no canned option).

Follow-up clarification on the sub-space mechanism:

| Option | Description | Selected |
|--------|-------------|----------|
| Project nesting (parent Project per target, version Projects nested inside) | Flat Space + parent/child Project hierarchy. | |
| **True nested Spaces** (Space-per-target under the Reviews Space) | Umbrella Space contains child Spaces. | ✓ |
| Naming-convention-only grouping inside one flat Space | All Projects flat; grouping by name prefix only. | |
| Something else — I'll describe | (Open option.) | |

**User's choice:** True nested Spaces — umbrella `Proteomics-Reviews` Space contains per-target child Spaces, each holding versioned Projects.
**Notes:** Verify-and-rollback gate (Pitfall 2 mitigation) is non-negotiable on top of the nested-Space shape — applies regardless of which destination model wins. **Plan-time research flag (load-bearing):** Datagrok's first-class support for nested Spaces is not 100% verified; researcher must confirm against the live platform before the planner commits, and the fallback shape (Project nesting under one flat Space) is documented in CONTEXT.md D-03.

---

## Republish UX flow

| Option | Description | Selected |
|--------|-------------|----------|
| Same dialog, pre-filled, with "v2 supersedes v1" banner | One menu entry; detect-and-prefill behavior. | ✓ |
| Distinct 'Update Share' menu entry | Two menu items; explicit but doubles surface. | |
| Always-fresh dialog — no pre-fill, no banner | Simplest helper, but accidental republish risk. | |
| Skip — use Claude's default | (Would default to option 1.) | |

**User's choice:** Same dialog, pre-filled, with "v2 supersedes v1" banner (option 1, recommended).
**Notes:** Keeps menu surface flat (no other Proteomics submenu pairs create/update entries). Banner mitigates Pitfall 4 (versioning ambiguity) by making the supersede chain explicit at OK time. Bidirectional pointers (`superseded_by` on prior + `supersedes` on new) for full audit trail.

---

## Enrichment carry (PUB-12, P2)

| Option | Description | Selected |
|--------|-------------|----------|
| Two DataFrames in the Project (protein DF + enrichment DF) | Reuses v1.2 `onCurrentRowChanged` cross-DF wiring. | ✓ |
| One DataFrame: enrichment folded into protein DF columns | Loses term-level p-values and term-side click affordance. | |
| Defer PUB-12 to a v1.4 hotfix — ship Phase 15 protein-only | Cut P2 to keep round-trip-test gate tight. | |
| Skip — use Claude's default | (Would default to option 1.) | |

**User's choice:** Two DataFrames in the Project (option 1, recommended).
**Notes:** Reuses the v1.2 enrichment subscription pattern — explicitly NOT a new `CampaignSelectionBus` (that abstraction is Phase 18's territory per Pitfall 13). When source DF has no enrichment, publish ships protein DF only — opportunistic carry per P2 status.

---

## Claude's Discretion

The following implementation details were captured in CONTEXT.md `### Claude's Discretion` without user input, per the philosophy that the user is visionary and Claude is builder:

- Direction column representation (string vs numeric)
- Sharer's friendly name source (`grok.user.current().friendlyName` or equivalent)
- Single-row metadata column shape (multiple typed one-row columns vs JSON string)
- Audit context panel triggering mechanism (panel decorator + autostart auto-dock)
- Round-trip test fixture choice (synthetic demo + Spectronaut Candidates)
- Slug sanitization rules (charset, length cap)
- Mailto body shape for PUB-13
- `isPublished(df)` helper signature (single boolean tag check on `proteomics.published === 'true'`)

---

## Deferred Ideas

Captured in CONTEXT.md `<deferred>` — preserved here in summary for retrospective audit:

- First-class target taxonomy (v1.5+)
- Publish-block-on-QC-fail policy (v1.5+)
- PDF report export (v1.5+, separate engineering project)
- Reviewer-side comment threads (v1.5+, use platform comments if needed)
- Per-target Space cleanup story (admin concern, not engineering)
- Email / Slack / webhook alerting on publish (out of scope, REQUIREMENTS.md)
- N-way compound comparison for published projects (v1.5+)
- Republish overwrite mode (explicitly rejected — supersede chain only)
- Publishing arbitrary non-DE-complete DFs (explicit success criterion 1 requires `proteomics.de_complete === 'true'`)

---

## Audience clarification (load-bearing for the whole phase)

Mid-discussion the user asked Claude to define "reviewer". Evidence pulled from PROJECT.md + ROADMAP.md + PITFALLS.md + the source todo confirms: **reviewer = consuming biologist** (primary persona). Drives the jargon audit (banned words: DataFrame, tag, semType, ACL, viewer factory), the trim allowlist's rationale (no raw intensities), the audit panel framing, the mailto wording (PUB-13), and the publish-destination model (biologists land via direct link, not by browsing Spaces). This pin lives at the top of CONTEXT.md `<domain>` as the single most decision-shaping fact in the phase.
