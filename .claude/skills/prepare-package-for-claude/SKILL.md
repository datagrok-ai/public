---
name: prepare-package-for-claude
description: Prepare a Datagrok package so that Claude (or any coding agent) can work on it effectively from a naive prompt
when-to-use: When the user wants to make a package "Claude-ready" — so that a developer can describe a feature in plain language and Claude implements it correctly without hand-holding
context: fork
effort: high
---

# Prepare a Datagrok Package for Claude

Help the user turn a Datagrok package into a codebase Claude can work on from simple prompts. The output is **not** exhaustive documentation — it's a minimal harness (per-package `CLAUDE.md` + possibly small supporting files) proven to carry real tasks.

## Philosophy

Three rules govern every decision:

1. **Real-task tests drive everything.** Write `CLAUDE.md` content only to address gaps a real agent actually hit on a real task. Don't speculate. Documentation that isn't load-bearing rots fast.
2. **Naive-dev prompts, not expert prompts.** The test is whether a developer who doesn't know the package internals can describe what they want and Claude does it. Prompts that leak internals are cheating.
3. **Keep the harness small.** Every line in `CLAUDE.md` must pull weight. Terse pointers beat thorough prose. If an agent can derive it from code in under 30 seconds, don't document it.

## The 4 phases

### Phase 1 — Baseline CLAUDE.md

Write a per-package `CLAUDE.md` that covers, at minimum:

- **Purpose** — what the package does, 2-3 sentences.
- **Architecture** — the 3-6 files an agent would need to reason about, with one line each describing their role and what lives inside them.
- **Glossary** — domain concepts → code types. Critical for scientific / domain-heavy packages.
- **Conventions** — anything that deviates from repo-wide Datagrok practice.

See [`templates/CLAUDE.md.template`](templates/CLAUDE.md.template) for a ready-to-fill skeleton.

**Do not** at this stage:
- Catalog every file. The agent can `ls`.
- Document conventions the repo-level `CLAUDE.md` already covers.
- Add a "common tasks" section unless you have real tasks to describe.
- Write a glossary entry for anything whose meaning is obvious from its type name.

### Phase 2 — Iterative harness loop (user-driven)

The core of the process. Phase 2 is **interactive**: the user runs tasks in fresh sessions manually and returns feedback. Do **not** simulate naive-dev runs inside this session — context contamination defeats the test.

**Per-task loop:**

1. **Propose 2–4 candidate task prompts.** Cover a mix of categories (feature, refactor, bug, docs) — see [`templates/phase-2-tasks.md`](templates/phase-2-tasks.md) for what makes a task informative and [`templates/naive-dev-prompt.md`](templates/naive-dev-prompt.md) for prompt discipline (no jargon, no file paths, no "use helper X"). Present them as a numbered list and stop.

   Every proposed prompt must end with an explicit feedback request so the fresh agent self-reports what it had to guess. Append this block verbatim to each task:

   ```
   When done, report:
   - Files changed, helpers reused.
   - Anything you had to guess or infer that the codebase / CLAUDE.md did not spell out.
   ```
2. **User selects one** and runs it manually in a **fresh Claude Code session** (no carry-over context). The user collects the agent's feedback — especially anything the agent reports it had to guess, infer, or search around for.
3. **User returns the collected feedback.** Review it and **triage each item**:
   - **Harness gap** — something the agent should have been able to find but couldn't. → Propose a minimum pointer to add to `CLAUDE.md`.
   - **Genuine product decision** — a judgment call with no "right" answer. → Do *not* document.
   - **Platform-wide quirk** — Datagrok / TS / Node, not package-specific. → Flag for repo-level docs or a UI skill; not this package's `CLAUDE.md`.
   Present the proposed `CLAUDE.md` edits to the user and **stop for approval**.
4. **User approves (or edits) the proposed changes.** Apply them to `CLAUDE.md`. Then **roll back the code changes** the fresh agent produced (the task was a harness test — the code artifact is secondary) so the next run starts from the same baseline plus the new harness.
5. **User re-runs the same task** in another fresh session with the updated `CLAUDE.md`. Return to step 3.
6. **Repeat the loop** until the task finishes with minimal or acceptable guesses. At that point, decide with the user whether to keep, revert, or modify the final code artifact.
7. **Continue with the next task** from the list in step 1 (or propose a new batch).

**When to stop:** the **user** decides when Phase 2 is done and Phase 3 can begin. Do not auto-advance. Good convergence signals to surface for the user: two consecutive tasks end with zero harness gaps, or the remaining guesses are all genuine product decisions.

Signals a task was well-designed:
- The agent self-reports 2–5 guesses. Zero guesses = task was too easy; >8 guesses = harness is still thin, or the task was too ambitious for this round.
- Each run surfaces *new* guesses, not re-runs of the same gap.

**Do not** in Phase 2:
- Spawn sub-agents to role-play the naive dev. Fresh sessions must be launched by the user.
- Edit `CLAUDE.md` before the user approves the proposed changes.
- Advance to the next task, or to Phase 3, without the user's explicit go-ahead.

### Phase 3 — Simplify (carefully)

With the harness proven, review the package as a whole and clean up what a reviewer would call out.

**Do NOT:**
- Refactor for style (indent, quotes, semicolons, brace placement).
- Change public/registered function signatures without explicit authorization.
- Touch code that works just because it looks ugly — have a concrete reason (bug, duplication, naming lie, dead code, measurable complexity reduction).

**General reporting rules (both steps):**
- Report findings **grouped by file**, with `file:line` and a one-line description per finding.
- Skip pure style noise (indent, quotes, semicolons, spacing, max-len, trailing whitespace) — call out semantic issues only.
- Filter out `package.g.ts` / `package-api.ts` — regenerated by `grok api`.
- Present a numbered list and ask which to apply (e.g. "all", "none", "1, 3, 5", or "skip"). Never fix before the user selects.

---

**Step 1 — Reviewer read**

Run the prompt below in **two independent sub-agent sessions** (same prompt, separate runs), then union the reports. Empirically, two Opus-4.6 runs of the same prompt overlap only ~50% of findings — running once loses half the achievable coverage for no real saving.

Prompt to run in each sub-agent:

> Analyze the source code in `public/packages/<PackageName>/src`. Find bugs, typos, code smells, and dead code. Suggest improvements and refactorings to make the code shorter, simpler, and easier to maintain. Group findings by file and include `file:line` for each. Skip generated files (`*.g.ts`, `*-api.ts`).

After both runs return:
- Merge reports. Dedupe by `file:line` (same line, same finding = one entry).
- For typos in **exported names** (class / function / type), flag separately — they're breaking changes; list call sites and ask before renaming. Typos in local names, decorator strings, and user-visible literals are safe to fix.
- Present the merged, numbered list to the user. Ask which to apply. Apply only selected ones.

---

**Step 2 — ESLint (optional, ask first)**

After Step 1 is applied, ask the user: *"Run eslint for anything Step 1 missed?"* If no, skip.

If yes, and `package.json` has a `lint` script, run:

```bash
npm run lint
```

If lint output is dominated by style noise, re-run with style rules disabled so only semantic findings remain:

```bash
npx eslint "src/**/*.ts" --rule '{"indent":"off","quotes":"off","semi":"off","max-len":"off","no-trailing-spaces":"off","object-curly-spacing":"off","curly":"off","comma-dangle":"off","eol-last":"off","padded-blocks":"off","space-before-function-paren":"off","keyword-spacing":"off","space-infix-ops":"off","no-multi-spaces":"off","no-multiple-empty-lines":"off","brace-style":"off"}'
```

Report findings grouped by file; present a numbered list; ask the user which to fix; apply selected; re-run to confirm.

If there is no `lint` script, skip. Do not install eslint or hunt for a bundled binary.

---

#### Applying findings

Each step above is self-contained: report → user selects → fix → next step.

Final validation: `npm run build` (catches type errors via ts-loader and regenerates `.g.ts`).

### Phase 4 — Cross-check

Acceptance test. Hand the package + new `CLAUDE.md` to a **fresh** Claude Code session (or a different model entirely) and give it **a new task** — not one from the Phase 2 portfolio. Re-running a Phase 2 task whose code you kept tests nothing; the agent would see the existing implementation and mimic it.

The Phase 4 task should:
- Be something the package has never seen — genuinely new feature, bug, or refactor.
- Ideally probe a surface Phase 2 didn't exercise (new category, new subsystem).
- Follow the same naive-dev prompt discipline as Phase 2.

**Reading the result:**
- Zero to one harness gap, all on genuinely new surfaces → the harness is real.
- Several gaps, or a gap on something Phase 2 *should* have covered → you converged too early. Go back into the Phase 2 loop.
- The agent couldn't even orient itself → the Phase 1 baseline is weak; rewrite the architecture / glossary.

This step is the one most likely to be skipped. Don't skip it. A harness that only works with the author in the room isn't a harness.

## Common pitfalls

- **Over-documenting on the first pass.** Write the Phase 1 baseline minimally. Every real gap can be added in Phase 2.
- **Leaking internals into prompts.** If your task prompt says "use `createCDDTableView`," you're not testing the harness, you're testing the prompt. Strip hints.
- **Treating every guess as a harness gap.** Many guesses are genuine product decisions. Ask: would a rule that eliminates this guess help the *next* task, or just this one?
- **Conflating simplification with harness work.** Phase 3 is after Phase 2 for a reason. Simplifying prematurely can mask gaps the harness would have surfaced.
- **Skipping Phase 4.** "It worked when I tried it" is not the bar.

## When the package is ready

- A real task from a new developer, stated naively, runs to completion in a fresh session.
- The `CLAUDE.md` has no speculative sections — every line was added because something broke without it.
- Two consecutive Phase 2 tasks end with "no harness gaps, only product decisions."
- The `v.next` changelog reflects what changed.

## Output artifacts

When you're done:
- `CLAUDE.md` at the package root.
- Possibly a few supporting `.md` files under `.claude/` if the package has a large feature area that deserves its own doc (rare — only if the Phase 2 loop demanded it).
- `CHANGELOG.md` `v.next` entries for any code changes.
- No leftover task code you don't intend to ship.
