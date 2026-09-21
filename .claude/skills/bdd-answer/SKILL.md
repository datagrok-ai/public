---
name: bdd-answer
description: Answer a "how do I …" question about the platform with a Gherkin scenario on @datagrok-libraries/bdd that is filmed into a how-to video (grok-bdd guide), so the answer is demonstrated, sent as a video plus numbered steps, and kept as a regression test
when-to-use: When a user or customer asks how to do something in Datagrok, when support wants a walkthrough video or GIF for a question, or when a help page needs its walkthrough regenerated
context: fork
effort: medium
argument-hint: "<the question, or a features/guides/*.feature path> [<package>]"
---

# Answering a question with a filmed scenario

The answer to "how do I open a CSV from my computer?" is a scenario: the steps a person takes,
in the vocabulary the bdd library already has, run against the stand and filmed. What goes back
to the person is the video (or GIF) and the numbered steps with pictures; what stays is the
feature, which the suite runs from then on — an answer that stops being true fails a test.

Rules: **the vocabulary first** (`npx grok-bdd list-steps` in the package's `bdd/`); a phrase the
library lacks is added as a binding or as a name/signal in the core, exactly as for a test
(`public/libraries/bdd/CLAUDE.md`); demo data only (`demog`, the package's `fixtures/`), never a
customer's file; nothing is committed without the lead's order. A question that needs a feature
the platform does not have is not a guide: say so, with what the closest scenario shows, and stop.

## 1. Find or write the scenario

1. Pick the package that owns the area (Browse, spaces, users → `packages/UsageAnalysis/bdd`;
   a viewer → the same, under `features/viewers/`; Bio, Chem, DiffStudio → their own `bdd/`).
2. Grep `features/` for a scenario that already does it — a test scenario often is the answer,
   only written for a checker rather than a reader. Prefer filming it as it is.
3. Otherwise write `features/guides/<slug>.feature`: tag `@guide` (and `@help:<page dir>` when it
   illustrates a help page under `public/help/`), a description that quotes the question, one
   scenario of 5–12 steps. Every `When` names the element a person would click, in the order they
   would; a `Then` after the decisive step shows the result. Step text is the caption of the video
   (`user clicks on X inside Y` → "Click on X in Y"), so it is written for the person asking.
4. `npx grok-bdd compile` then `npx grok-bdd run generated/guides/<slug>.test.ts` until green; a
   step that needs a wait or a pixel is a missing signal in the core, not a `sleep`.

## 2. Film it

```
cd public/packages/<Pkg>/bdd
MSYS_NO_PATHCONV=1 npx grok-bdd guide features/guides/<slug>.feature [--gif]
```

Watch `guides/<feature slug>/<scenario slug>/guide.mp4` once. What to fix and where:

- the pointer goes to the wrong place → the step located a scope, not the element: make the
  phrase name the element (`"Open local file" icon inside browse toolbar`), or add the name to
  the core;
- a dialog or balloon is missing from the "after" picture → `--settle 1000`;
- a step reads badly in the caption → reword the step; a verb the caption leaves in the third
  person goes into the `VERBS` map of `src/runtime/guide.ts`;
- a setup step shows up that should not (or the reverse) → a `Given` that touches nothing is
  skipped; make it a `When` to show it.

## 3. Send it, and keep it

- Reply with `guide.mp4` (or `guide.gif` where a video does not embed) and the body of
  `steps.md`, whose pictures are the `step-NN.png` next to it.
- Add a line to the package's `bdd/features/guides/INDEX.md`: the question, the feature, the date,
  the help page it illustrates if any. Create the file with a one-line header if absent.
- Leave the feature and its generated spec uncommitted for the lead, with the list of core or
  library changes it needed.

For help pages: `npx grok-bdd guide --help-pages` re-films every `@help:` tagged feature of the
package and copies the GIFs into `public/help/<page dir>/img/`; the page embeds them as
`![…](img/<scenario slug>.gif)`.
