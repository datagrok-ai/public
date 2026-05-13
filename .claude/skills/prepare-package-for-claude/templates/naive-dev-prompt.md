# Writing naive-dev prompts

Phase 2 task prompts must be what a real developer (who does not know the package internals) would write. The whole point is to test whether the **harness** supports the agent — not whether a carefully-written expert prompt does.

## Principles

1. **State the outcome, not the implementation.** The prompt says what should happen from the user's perspective, not which function to call.
2. **Use domain vocabulary, not code names.** "The last search the user ran" — not "persist `SeachEditor.getParams()` result."
3. **Short is honest.** A real ticket is 2–5 lines. If you need a paragraph to explain the task, the task itself is unclear or the prompt is leaking internals.
4. **Include only non-negotiable constraints.** "Must persist across reloads" is a requirement. "Use `grok.userSettings`" is cheating.
5. **Ask the agent to self-report guesses.** Every prompt should end with a request for what the agent had to infer — that's the harness signal.

## Anti-patterns

- ❌ Mentioning specific functions, files, or helpers by name: `"use createCDDTableView"`.
- ❌ Describing the implementation approach: `"create a sync wrapper that calls queryBatches with page_size=100"`.
- ❌ Pre-solving the design: `"store this in grok.userSettings keyed by vaultId"`.
- ❌ Listing edge cases the agent should handle *before seeing the code*: `"make sure to handle the case where protocols haven't loaded yet"`.
- ❌ Long preambles explaining why the task matters. The agent doesn't need that; the task description stands on its own.

## Good prompt examples

### Feature (small) — good
> Add a Batches tab to the app, alongside the existing Protocols, Collections, Saved Searches, and Molecules tabs.
>
> Each batch row should link back to its parent molecule. Routing to the new tab via URL must work the same way the other tabs do.
>
> When done, report: files changed, helpers reused, and anything you had to guess or infer that the codebase / CLAUDE.md did not spell out.

Why it's good: names the outcome and the UI placement, specifies the linking requirement (a real product constraint), mentions routing (a real non-negotiable), doesn't name any functions.

### Feature (cross-cutting) — good
> Remember the last search the user ran and restore it next time they come back. Per-vault: searches in vault A should not show up in vault B. Must persist across reloads.

Why it's good: three lines. States the outcome and two scoping constraints. Says nothing about where state lives, which API to use, or how to hook in.

### Refactor — good
> There's a repeated pattern across several async-export calls: kick off the job, get the export id, poll for results. Collapse that into a helper. No behavior change.

Why it's good: describes the pattern in abstract terms, names the constraint ("no behavior change") that pins down success, doesn't say where to put the helper or how to name it.

## The reporting coda

Every prompt should end with:

> When done, report: files changed, helpers reused, and anything you had to guess or infer that the codebase / CLAUDE.md did not spell out.

That last clause — **"guess or infer"** — is what turns the run into a harness test. Without it, the agent will quietly fill in blanks, you'll never see what it had to invent, and the harness stays broken.

## When to relax these rules

Once, at the end, in Phase 4 (cross-check). That's where you *want* to see whether the harness is complete. If even a realistic-but-minimal prompt fails, the harness isn't done.

During iterations (Phase 2 loop), always keep prompts minimal. A harness that only works with elaborate prompts isn't a harness — it's a carefully-staged demo.
