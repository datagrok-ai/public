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

### Phase 2 — Iterative harness loop

The core of the process. Run **2–4 real tasks** against the package. After each task, close the loop: verify, record, update the harness, repeat.

**Per-task loop:**

1. **Pick the next task.** Cover a mix of categories (feature, refactor, bug, docs) — see [`templates/phase-2-tasks.md`](templates/phase-2-tasks.md) for what makes a task informative and examples from each category.
2. **Write a naive-dev prompt.** A developer who just joined and doesn't know the package internals. No jargon, no file paths, no "use helper X." See [`templates/naive-dev-prompt.md`](templates/naive-dev-prompt.md) for the principles and anti-patterns.
3. **Run the task in a fresh Claude Code session.** No carry-over context. Do not course-correct during the run — the point is to see what the harness alone supports.
4. **Assess the result.** Use the scorecard in `templates/phase-2-tasks.md`. What was right? What did the agent report it had to guess or infer?
5. **Triage each guess.** For each item the agent flagged:
   - **Harness gap** — something the agent should have been able to find but couldn't. → Update `CLAUDE.md` with the minimum pointer that would have prevented the guess.
   - **Genuine product decision** — a judgment call that has no "right" answer. → Do *not* document. The agent made a reasonable call; noting it in `CLAUDE.md` would over-fit.
   - **Platform-wide quirk** — a Datagrok / TypeScript / Node thing, not package-specific. → Flag for repo-level documentation or a UI skill; do not put in this package's `CLAUDE.md`.
6. **Apply the update.** Edit `CLAUDE.md` with the minimum pointer. Commit.
7. **If the run produced new code**, decide: keep, revert, or modify. The task was a harness test — the code artifact is secondary. If the agent did a good job and harness gaps are acceptable, keep the work. If not, revert and repeat with the updated CLAUDE.md.

**When to stop:** after 2–5 tasks, or when two consecutive tasks produce **zero** harness gaps (only product decisions remain). That's the convergence signal.

Signals a task was well-designed:
- The agent self-reports 2–5 guesses. Zero guesses = task was too easy; >8 guesses = harness is still thin, or the task was too ambitious for this round.
- Each run surfaces *new* guesses, not re-runs of the same gap.

### Phase 3 — Simplify (carefully)

With the harness proven, do a targeted simplification pass. Focus on:

- **Duplicated patterns** → extract a helper with a discoverable name (future Claude will find it).
- **Unsafe primitives** → replace with safer built-ins (e.g. hand-rolled URL building → `URLSearchParams`).
- **Dead code** → remove and note the removal in CLAUDE.md if the function was named anywhere non-obvious.
- **Latent bugs surfaced by static analysis** (see below).

**Do NOT:**
- Refactor for style.
- Change public/registered function signatures without explicit authorization.
- Touch code that works just because it looks ugly.
- Run `/simplify` over everything without a specific target.

Every change should be one of:
- A move the Phase 2 tasks revealed as useful,
- A fix for a latent bug (typo, dead branch, wrong variable passed, missing `await`, unreachable code, etc.),
- A named helper that replaces a ≥3-site duplicated pattern.

If none of those apply, skip it.

#### Static analysis — run each step separately, report, then ask

Run each step below as a **distinct, interactive loop**: run the check → report results to the user → let the user pick which findings to fix → apply only those → move on. Do **not** batch all steps together, and do **not** auto-fix without the user's selection.

**General reporting rules (all steps):**
- Report findings **grouped by file**, with `file:line` and a one-line description per finding.
- Skip pure style noise (indent, quotes, semicolons, spacing, max-len, trailing whitespace) — call out semantic issues only.
- Filter out `package.g.ts` / `package-api.ts` — regenerated by `grok api`.
- After reporting, present a numbered list of fixable findings and ask the user which to apply (e.g. "all", "none", "1, 3, 5", or "skip"). Never fix before the user selects.
- After applying fixes, re-run the same check to confirm it's clean (or that only user-skipped findings remain), then move to the next step.

---

**Step 1 — TypeScript type-check**

```bash
npx tsc --noEmit -p tsconfig.json
```

- Report errors grouped by file.
- If zero errors: say so in one line and move on.
- If errors exist: list them, ask the user which to fix, apply selected fixes, re-run to confirm.

---

**Step 2 — ESLint (semantic findings only)**

If `package.json` defines a `lint` script, run it first:

```bash
npm run lint
```

If there is no `lint` script, skip this step entirely. Do not install eslint or hunt for a bundled binary.

When lint output is dominated by style noise, re-run eslint with style rules disabled so only semantic findings remain:

```bash
npx eslint "src/**/*.ts" --rule '{"indent":"off","quotes":"off","semi":"off","max-len":"off","no-trailing-spaces":"off","object-curly-spacing":"off","curly":"off","comma-dangle":"off","eol-last":"off","padded-blocks":"off","space-before-function-paren":"off","keyword-spacing":"off","space-infix-ops":"off","no-multi-spaces":"off","no-multiple-empty-lines":"off","brace-style":"off"}'
```

Focus on rules like `no-unused-vars`, `no-throw-literal`, `no-unreachable`, `@typescript-eslint/require-await`, `no-useless-catch`, `no-constant-condition`, etc.

Report semantic findings grouped by file (`file:line: rule — message`). Then:
1. Identify which are safely fixable (e.g. `no-throw-literal` → wrap in `new Error(...)`; unused local imports → remove).
2. Flag ones that need judgment (an unused *exported* symbol may be called externally; ask before removing).
3. Present a numbered list, ask the user which to fix, apply only the selected ones.
4. Re-run the semantic-only eslint command to confirm.

---

**Step 3 — Manual semantic checks**

Run these **every time**, even when eslint passes — rules are often disabled and eslint misses intent-level bugs. Grep + read; do not guess.

Checks:
- `async` functions with no `await` inside (drop `async` or actually await something).
- Empty `try { ... } catch (e) { throw e; }` — pure ceremony, delete.
- `continue` / `break` / `return` as the *last* statement of a loop or function — dead.
- Unreachable code after `return` / `throw`.
- `if (x)` immediately followed by code that would run for falsy `x` anyway (missing `return` / `else`).
- Wrong variable in params — e.g. `paramsStringFromObj({async: true})` where it should have been `paramsStringFromObj(params)`. Audit every call site where a function takes `params` and also has a literal object in the body.
- Mutations on `arr` but reads from `arr.filter(...)` (or vice versa) in the same block.

Collect findings across all checks, then **report a single summary** grouped by file with `file:line` + one-line description + severity guess (bug / dead code / stylistic). Present a numbered list, ask the user which to fix, apply only the selected ones, and validate with `npx tsc --noEmit`.

---

**Step 4 — Typo pass**

CDD Vault taught us these matter because they reach users. Grep for obvious misspellings in:
- **Exported names** (class / function / type names).
- **Decorator strings** (`category`, `description`, `name`).
- **User-visible literals** (`grok.shell.error/warning/info`, button labels, tooltips).
- **Comments** only when they name a symbol.

Patterns: doubled letters (`reorderColummns`, `Filelds`), flipped pairs (`Seach`/`Search`, `Serach`/`Search`, `Stastics`/`Statistics`).

Report findings grouped by file, each annotated as:
- **Local** (param / private method / decorator string) — safe to rename.
- **Exported** — breaking change; list call sites, flag for the user to decide, do not auto-rename.

Present a numbered list, ask the user which to fix, apply only the selected ones.

---

#### Applying findings

Each step above is self-contained: report → user selects → fix → validate → next step. At the end of all four steps, give a brief roll-up:
- **Bugs** fixed (wrong behavior today).
- **Dead code** removed (unused, unreachable). Note: registered `@grok.decorators.func` may be called externally — always flag, never auto-remove.
- **Typos** fixed (local) vs flagged (exported).

Final validation: `npx tsc --noEmit -p tsconfig.json` (filter out `package.g.ts` / `package-api.ts`).

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
