# Publishing (Share for Review) — design rationale

Why the Share-for-review subsystem is built the way it is, and **which parts are permanent
architecture vs. workaround weight that the base platform could absorb.**

This is the reference behind the [personas & capabilities](./personas-and-capabilities.md) doc's
enforcement claims. That doc says *what* the read-only boundary is; this says *why the code takes
the shape it does to make that boundary real* — and where we'd rather the platform did the work.

Audience: anyone questioning the complexity of `src/publishing/`, and anyone planning what
platform capabilities to grow so this package can shed code. Per our standing preference, a
workaround defended as *temporary, with a named platform successor* is a stronger position than
one defended as permanent design.

## The one-sentence defense

Most of `src/publishing/` (~1,900 lines across 9 files) is not "sharing." It is the machinery to
**guarantee three properties** a compliance-sensitive review needs, none of which the platform
gave us for free at build time:

1. **Provably read-only** — the reviewer group *cannot* alter the interpretation.
2. **Provably minimized** — only the interpretation ships; raw/again-identifying data does not.
3. **Provably reload-survivable + versioned** — what the reviewer opens is what was published,
   every version is immutable, and supersession is auditable.

The naming scheme and the use of Spaces (the two most visible choices) are *downstream* of these
guarantees, not the point of them.

## Design choices, grouped by the guarantee they serve

Legend for **Verdict**: 🟢 keep in package (genuinely Proteomics-specific) · 🟡 push to platform
(generic capability we had to hand-roll) · 🔵 hybrid (platform does the mechanism, package keeps the policy)

### A. Read-only is enforced, not assumed

| # | Choice | Where | Why | Verdict |
|---|---|---|---|---|
| A1 | Grant **View directly on the Project**, not via Space inheritance | `publish-project.ts:323` | Spike A2: `permissions.get(project)` does **not** surface Space-inherited grants, so an inherited grant is unverifiable. | 🟡 |
| A2 | **Verify-and-rollback ACL gate** — after granting, re-read permissions across three rings (project, child Space, umbrella Space); any `edit`/`share`/`delete` for the reviewer group ⇒ delete the project and abort | `publish-project.ts:331-370` | The feature refuses to ship what it cannot prove is read-only. This is the single most important — and most expensive — thing the code does. | 🔵 |
| A3 | Share to a **group/team, never individuals**; eligible groups exclude hidden, personal, and "All users" | `reviewer-groups.ts:9-19` | Review is a team activity; "All users" would defeat the point of a scoped share. | 🔵 |

### B. Data minimization — a trimmed snapshot, not a live link

| # | Choice | Where | Why | Verdict |
|---|---|---|---|---|
| B1 | **Deep clone; source never mutated** — the publish path operates only on the clone | `trim-dataframe.ts:75-120` | Reviewers get a frozen point-in-time copy; later edits to the working analysis cannot leak. | 🟢 |
| B2 | **Column allowlist** (Protein ID, Gene, log2FC, p, adj.p, significant, [direction, quantities]); everything else dropped | `trim-dataframe.ts:108-120` | The generalized form of the 1.2.0 "client-name removal." A *whitelist* means a new, unexpected column is excluded by default — safer than a blocklist. | 🟢 |

### C. Provenance is embedded and redundant

| # | Choice | Where | Why | Verdict |
|---|---|---|---|---|
| C1 | Rich metadata baked into the snapshot (who/email/DE method/thresholds/version/UUID/date) | `publish-state.ts:71-92`, `trim-dataframe.ts:158-174` | The snapshot is self-describing — provenance travels with the data, not a side channel. | 🟢 |
| C2 | **Dual-write belt-and-braces** — every field written as *both* a DataFrame tag *and* a hidden `_meta_*` column; read column-first, tag-second | `publish-state.ts:24-60,137-225` | Defends against serialization paths that strip tags but keep columns (or vice versa) — an observed failure mode. | 🟡 |
| C3 | Numeric metadata stored **as strings** | `trim-dataframe.ts:153-167` | `addNewFloat` is Float32; string storage preserves exact double thresholds. | 🔵 |

### D. Versioning is first-class and non-destructive

| # | Choice | Where | Why | Verdict |
|---|---|---|---|---|
| D1 | **Opinionated naming** `Proteomics-Review-<slug>-v<N>-<date>`; structural suffix fixed (not configurable) | `publish-project.ts:244-250`, `publish-settings.ts` | The `-v<N>-` segment is parsed to detect prior versions; letting it vary would break version detection. Prefix/Space *are* configurable. | 🔵 |
| D2 | **Republish detection** by name pattern → bump to v(N+1); prior version stays available (never overwritten) | `publish-state.ts:279-328`, `publish-project.ts:138-140` | Immutable, auditable version history. | 🔵 |
| D3 | **Bidirectional supersede pointers** (`supersedes`/`supersededBy`) on both projects, dual-written to project options *and* the DataFrame | `publish-project.ts:436-478` | Either version can be navigated to the other; survives whichever persistence path holds. | 🟡 |

### E. Round-trip fidelity is verified

| # | Choice | Where | Why | Verdict |
|---|---|---|---|---|
| E1 | **Reopen-and-assert gate** (default on) — `closeAll` → reopen → assert published shape survives a reload; rollback on failure | `publish-project.ts:372-434`, `assert-published-shape.ts` | Proves the reviewer sees exactly what was published, before telling the analyst it worked. | 🟡 |
| E2 | **Formula-line self-heal** — a serializer strips volcano threshold lines; re-apply from tags and re-check once before rollback | `publish-project.ts:401-433`, `post-open-recovery.ts` | Works around look-config loss on project round-trip (Phase 13 evidence). | 🟡 |
| E3 | **Re-render viewers on the trimmed snapshot** so the shared dashboard is self-contained | `publish-project.ts:265-295` | The published project must stand alone — no dependency on the analyst's session tables. | 🔵 |
| E4 | **Post-open recovery autostart** — re-applies formula lines and re-wires enrichment↔volcano cross-highlight on reopen | `post-open-recovery.ts`, `package.ts:1010-1037` | Interaction subscriptions and look config don't survive serialization; rebuild them on open. | 🟡 |

### F. Where content lands + idempotency

| # | Choice | Where | Why | Verdict |
|---|---|---|---|---|
| F1 | Publish to **named Spaces, not a personal workspace** | `publish-project.ts:178-239` | A shared, governed location the reviewer team can browse — not tied to one user's private space. | 🔵 |
| F2 | **Two-level Space hierarchy** — umbrella `Proteomics-Reviews` + per-project child Space | `publish-project.ts:178-239` | Groups all reviews under one root, one child per project — a navigable tree instead of a flat pile. | 🔵 |
| F3 | **Create-first-then-enumerate idempotency** for Spaces; fall back to `spaces.list()` on "already exists" | `publish-project.ts:184-238` | The smart-filter API does not reliably surface Space-flagged Projects, so enumeration is the reliable lookup. | 🟡 |

### G. Governance + workflow UX

| # | Choice | Where | Why | Verdict |
|---|---|---|---|---|
| G1 | **Controlled project vocabulary** — analyst picks from an admin-maintained list, never free-types | `share-dialog.ts:57-73`, `publish-settings.ts` | Keeps shared results filed under sanctioned project names; a governance control. | 🟢 |
| G2 | **Precondition gate** — share only available after DE completes | `share-dialog.ts:31-34` | There is no interpretation to review until DE has run. | 🟢 |
| G3 | **Reactive confirmation summary** — exact resulting name/team/supersede status shown before OK | `share-dialog.ts:109-151` | No surprises; the analyst sees the outcome before committing. | 🟢 |
| G4 | **Workspace restoration** after publish — tear down verification artifacts, restore the analyst's original view | `publish-project.ts:480-506` | The verify/supersede steps `closeAll` repeatedly; put the analyst back where they were. | 🔵 |
| G5 | **Full package-settings configurability** (Space, prefix, default team, vocabulary, verify default) | `publish-settings.ts` | Retarget a deployment without a code change. | 🟢 |

## Mapping to platform capabilities

Reading the 🟡/🔵 rows top-down, the workaround weight collapses to a small number of platform
gaps. If the platform closed these, the package could delete or thin most of `src/publishing/`.

| Platform gap (observed) | Evidence in code | What a platform capability would give us | Package code it retires |
|---|---|---|---|
| **P1. No first-class "publish read-only snapshot to a group" primitive** | The entire orchestrator exists to compose this from lower-level parts | A single `publishReadOnly(entity, group, {snapshot, immutable})` that grants view-only, freezes, and returns a stable handle | Most of `publish-project.ts` (A1, A2, F1–F3) |
| **P2. `permissions.get` doesn't reflect Space-inherited grants; may not report share/delete** | Spike A1/A2 notes at `publish-project.ts:113-118,331-370` | A permissions read that returns the *effective* ACL (own + inherited) across all verbs | A2 three-ring check collapses to one `get`; A1 direct-grant workaround unneeded |
| **P3. Serialization strips look config (formula lines) and interaction subscriptions** | E2/E4; `post-open-recovery.ts`, `assert-published-shape.ts:11-17` | Viewer look + formula lines + declared subscriptions persist across `save→open` | E2, E4 self-heal and recovery hook; part of E1's assertions |
| **P4. Tag map / metadata not reliably durable across `clone`/`save`** | C2 dual-write; `publish-state.ts:227-255`, `trim-dataframe.ts:20` | Durable, typed entity/DataFrame metadata that survives clone and round-trip | C2 `_meta_*` columns; the column-first/tag-second reader in `getPublishedMetadata` |
| **P5. Smart-filter can't reliably find Space-flagged Projects** | F3 note at `publish-project.ts:184-203` | `dapi.projects.filter(isSpace=…).first()` that actually works | F3 create-first-then-enumerate dance |
| **P6. No immutable/versioned entity semantics** | D1–D3 name-parsing + manual supersede pointers | Native entity versioning with supersede links | D2/D3 name-pattern versioning and dual-written pointers; D1's suffix could relax |
| **P7. No round-trip integrity check** | E1; `assert-published-shape.ts` | If P3+P4 land, publish *is* faithful by construction and E1 becomes optional | E1 gate downgrades to a test-only assertion |

### What stays in the package regardless

The 🟢 rows are the genuine Proteomics policy and don't belong in the platform:

- **B1/B2** — *which* columns constitute "the interpretation" (the allowlist) is domain knowledge.
- **C1** — *which* provenance fields matter (DE method, FC/p thresholds) is domain knowledge.
- **G1/G2/G3/G5** — the vocabulary, the DE precondition, the confirmation copy, the settings surface.

Even under a fully-featured platform, the package would still call a `publishReadOnly` primitive,
hand it a domain-trimmed snapshot with domain provenance, and gate the UI on DE completion. That
residue is small — roughly `trim-dataframe.ts` + `share-dialog.ts` + `publish-settings.ts`.

## Recommended platform asks, prioritized

1. **P1 + P2 together** — the read-only-to-group primitive with a trustworthy effective-permissions
   read. Highest leverage: retires the bulk of the orchestrator *and* the code we're least confident
   in (the ACL gate can only be as correct as `permissions.get`).
2. **P3 + P4** — durable look/subscription/metadata across round-trip. Retires the self-heal and
   recovery machinery, which is the most brittle code in the subsystem.
3. **P6** — native immutable versioning. Retires the name-parsing version scheme and lets the naming
   suffix relax.
4. **P5, P7** — fall out naturally once the above land.

## Honest caveats to raise proactively

- **The ACL gate is only as strong as `permissions.get`.** If P2 is real (the API may not report
  `share`/`delete`), the three-ring check at `publish-project.ts:358-361` cannot detect those leaks
  today. Worth verifying against the current platform release before claiming the guarantee is
  airtight. *(This claim is inferred from the A1 spike note and should be re-checked live.)*
- **Several verify/rollback steps `closeAll` the user's session.** G4 restoration mitigates it, but
  it's inherent to doing round-trip verification in the same client session; P7 landing would remove
  the need entirely.
- **This is workaround-heavy by design, and that's defensible** *because* each workaround has a named
  platform successor above — not because the complexity is intrinsic to sharing.
