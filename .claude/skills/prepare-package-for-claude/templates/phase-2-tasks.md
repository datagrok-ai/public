# Phase 2: Task portfolio design

Pick 2-4 tasks spanning the categories below. Each task is an independent harness test.

## What makes a task informative

- **Realistic.** Something a real developer would actually ask for. "Add a useless dashboard" isn't.
- **Probes a specific surface.** Every task should exercise a distinct slice of the package — loading, routing, persistence, UI registration, API layer, tests, etc. Don't run five similar feature-adds; vary the surface area.
- **Has a clear done signal.** You must be able to tell "did the agent succeed" within a few minutes. Vague exploratory tasks don't give harness signals.
- **Not trivial.** A one-liner doesn't surface gaps.
- **Not a research project.** If it takes the agent hours and needs architectural decisions, you're testing the agent, not the harness.

## The four categories

### Feature (small)
Add a self-contained feature that fits an existing pattern. Tests whether the agent can **find and reuse** existing patterns.

*Signals a good feature task probes:*
- Pattern discoverability: does the agent find the helper that already solves 80% of the task?
- Registration: does the agent wire the new thing in correctly (routes, tabs, functions, etc.)?
- Naming and placement: does the agent create files where other similar files live?

*Example shapes:*
- "Add a <new tab> to the <existing tabbed app>."
- "Add a <new info panel> that shows <X> when a <Y> is selected."
- "Support <new file format> in the <existing viewer>."

### Feature (cross-cutting)
A feature that touches ≥2 subsystems (UI + persistence, data + routing, etc.). Tests whether the agent **connects** subsystems correctly.

*Signals a good cross-cutting task probes:*
- Does the agent know *where* state should live for each subsystem?
- Does it compose existing APIs rather than invent new ones?
- Does it respect scoping rules (per-user, per-vault, per-session)?

*Example shapes:*
- "Remember the user's last <thing> and restore it on reload."
- "Share <current view state> via URL."
- "Add <X> to <Y>, and reflect it in <Z>."

### Refactor (no behavior change)
Collapse a duplicated pattern into a helper. Tests whether the agent can **move code safely** — the strictest test of "does the harness support non-trivial code transformation?"

*Signals a good refactor task probes:*
- Does the agent correctly identify all call sites?
- Does the new helper's API cover all of them?
- Does the agent avoid regression (the `noEmit` type-check should pass; behavior should be identical)?

*Example shapes:*
- "Extract a helper for <the three-step pattern> duplicated in <N> places."
- "Unify <function A> and <function B> into one shared core with per-call-site transforms."
- "Replace hand-rolled <thing> with the standard library equivalent."

### Bug-shaped
A specific, narrow bug with a reproducer. Tests the agent's **diagnostic** loop.

*Signals a good bug task probes:*
- Does the agent reproduce before fixing?
- Does it find the root cause (not just a symptom patch)?
- Does it write a test that would have caught the bug?

*Example shapes:*
- "<Specific input> causes <wrong output>. Fix and add a test."
- "<Feature> breaks when <edge case>. Fix."
- "<Error message> is misleading — fix the real cause."

### Docs-from-code
Produce a doc artifact from code alone. Tests whether the agent can **read + summarize** accurately.

*Signals a good docs task probes:*
- Does the agent read all relevant call sites?
- Does it produce something concise and correct, or vague and padded?
- Does it catch the non-obvious cases?

*Example shapes:*
- "Document the URL path scheme supported by <app>."
- "List every <function type> this package registers and what triggers each."
- "Produce a table of <entity type> → <where it's used>."

## Assessment scorecard

For each task, capture:

| Signal | Observation |
|---|---|
| Used the canonical helper / pattern unaided | ✅ / ⚠️ / ❌ |
| Touched only files it should have (no `.g.ts` edits, no unrelated refactors) | ✅ / ⚠️ / ❌ |
| Completed the task (if a feature: it works; if a refactor: no regressions; if a bug: root cause fixed) | ✅ / ⚠️ / ❌ |
| Number of self-reported guesses | integer |
| Guesses categorized: **harness gap** / **product decision** / **platform quirk** | per-item label |
| Agent's guess list (verbatim) | paste |

Only the **harness-gap** items generate CLAUDE.md edits. The others don't.

## Triage rules for guesses

When the agent flags something as "inferred / not spelled out":

**It's a harness gap if:**
- The answer is in the code and the agent couldn't find it within a reasonable search.
- The answer is a load-bearing convention the agent couldn't have known without being told.
- Two different reasonable agents would pick different approaches, but only one is correct for this codebase.

*Fix: add the minimum pointer to `CLAUDE.md`. One sentence is usually enough.*

**It's a product decision if:**
- There's no "right" answer — either would work, the choice is a taste call.
- The agent's choice was reasonable and matched what a lead would have picked.
- Documenting the answer would over-fit (it'd only apply to this one feature).

*Fix: nothing. Move on.*

**It's a platform quirk if:**
- It's a general Datagrok / TypeScript / browser behavior, not package-specific.
- Documenting it here would duplicate something that belongs in a repo-level skill or `CLAUDE.md`.

*Fix: flag it for the repo-level skill author. Do not add to this package's `CLAUDE.md`.*

## Convergence signal

Stop running Phase 2 tasks when **two consecutive tasks produce zero harness gaps** (only product decisions). That's the sign that the harness now carries the load.

If you've done 5 tasks and still see harness gaps every time, the harness is thin — widen scope and keep going.
