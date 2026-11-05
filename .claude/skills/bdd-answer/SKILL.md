---
name: bdd-answer
description: Answer a "how do I …" question about the platform with a Gherkin scenario on @datagrok-libraries/bdd that is filmed into a how-to video (grok-bdd guide), so the answer is demonstrated, sent as a video plus numbered steps, and kept as a regression test; the same scenarios, filmed in help mode, are the GIFs of the help pages
when-to-use: When a user or customer asks how to do something in Datagrok, when support wants a walkthrough video or GIF for a question, or when a help page needs a GIF recorded or replaced
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

What a guide shows, and what the runtime guarantees (`grok-bdd guide` sets it up, nothing to pass):

- **The full shell.** A guide runs with simple mode off — the menus, view tabs and panels as a
  person has them. A test page runs in simple mode; a guide never does, filmed or in a plain
  `grok-bdd run`: its second step is `And simple mode is off`, right after `Given user is logged
  in` (the compiler refuses a `@guide` feature without it). The step is silent, not in the video.
- **Every step a person would take is a UI step**, filmed as a gesture: a `When` names the element
  to click, hover, drag or type into, and the pointer goes there. An API step is only for what
  the person already has when they ask (the open tables); a view switch is a click on the view's
  tab (`user clicks on spgi-100 tab`), which the full shell shows. Opening the tables IS in the
  video: each `Given user opens X dataset` shows the table it opened under its caption. Only the
  login and a step that changed nothing on the page are left out (a help GIF leaves out every
  `Given`: section 4).
- **A menu path is shown stop by stop**: the group in the bar lights and opens, then each item on
  the way, then the leaf; the same for a context menu. A walk that lights only the bar means a
  runtime path did not report its stops: it calls `guide.hop` for each, as `pickTopMenu` and
  `pickMenuPath` do — fix the runtime, not the feature.
- **A choice is shown being made**: a native select opens its list, the option is typed so the
  list highlights it, Enter takes it; a column selector opens its picker, the name is typed into
  its search short of its last letter (a complete unique name is taken on the spot) and the row
  it leaves is clicked. Tests take the same steps through `selectOption` and Enter; only a guide
  walks the list (`selectNative`, `typeInColumnGrid` in `src/runtime/gestures.ts`).
- **A control that shows only under the pointer is hovered in a step of its own.** A filter panel's
  +, a filter card's search icon and a viewer's title-bar icons appear only while the pointer is over
  their panel (a guide films without the test styles that keep them shown). A step's first picture is
  taken before its gestures, so a hover inside the clicking step never shows the control and the
  press lands on empty space: write `When user hovers over filter panel` before `user adds a card for
  …`, `user hovers over scatter plot viewer` before `user clicks on settings icon of scatter plot
  viewer`. A click on such an icon with no hover before it fails ("present, none visible").
- **Text goes into a box the reader has seen clicked.** When another gesture left the focus in the
  box (a card's search icon opens and focuses its search), click into the box in a step of its own
  before typing or pasting, so the video shows where the text goes.
- **The pointer rests on the lit element before every click**; only an icon-sized target (28 px
  or less each way) is zoomed into. The caption sits above the page, clear of a player's timeline.
- **Every press is marked where it landed, under the pointer's tip**: a yellow dot and ring for
  the left button, a green one for the right, twice for a double-click. The place is the one the
  page received — every real press, release and move is logged in the page, whatever sent it (a
  locator's own click, the page's mouse, a drag) — never the centre of the element the step named.
- **The pointer never skips: each movement starts where the previous action ended.** The renderer
  carries the pointer from frame to frame, through drags and menu stops, and refuses to make a
  video in which it skips ("the pointer skipped N px") — a skip is a renderer bug to fix, never a
  frame to accept. Since every move shows, a gesture that sends the pointer across the screen for
  a test's sake (to the page's corner, to clear a hover) steps aside nearby in a guide instead
  (`besidePicker` in `gestures.ts`).
- **Only what a person would look for is in the video.** A `Then` shows when it names something
  on the page — a dialog, a column, a row count, a value, a legend item's color. The checks a
  test needs and a person does not are filmed out by the `HIDDEN_CHECKS` patterns in
  `src/runtime/guide.ts`: error and balloon floors, server state, viewer readings and pixels,
  "than before" claims, property bags, widget counts, task-bar and command bookkeeping, values
  matched against a regular expression (a guide names a new column exactly, not by a pattern). Keep
  them in the feature (they are the test); a new kind of bookkeeping check goes into that list.
- **Set-up a person never does is not filmed**: a pinned setting (`the molecule sketcher is …`) and
  a package's readiness wait are silent (`silent` from `@datagrok-libraries/bdd/runtime` for a
  package's own step), and a wait is written as a `Then` (`… should have finished updating`), not
  an `And` under a `When`, or it is filmed as an action.

## 1. Find or write the scenario

1. Pick the package that owns the area (Browse, spaces, users → `packages/UsageAnalysis/bdd`;
   a viewer → the same, under `features/viewers/`; Bio, Chem, DiffStudio → their own `bdd/`).
2. Grep `features/` for a scenario that already does it — a test scenario often is the answer,
   only written for a checker rather than a reader. Prefer filming it as it is.
3. Otherwise write `features/guides/<slug>.feature`: tag `@guide` (and `@help:<page dir>` when it
   illustrates a help page under `public/help/`), a description that quotes the question, one
   scenario of 5–12 steps. Every `When` names the element a person would click, in the order they
   would; a `Then` after the decisive step shows the result. Step text is the caption of the video
   (`user clicks on X inside Y` → "Click on X in Y"), so it is written for the person asking; where a
   step must stay precise for the test (an atom's index, a drag in pixels, a reading), a `# caption: <text>`
   comment on the line right above it is what the video says instead ("Click the acid's carbon"), and a
   check captioned so is shown even when `HIDDEN_CHECKS` would leave it out. A guide about one widget
   (a sketcher) starts with it open from a silent set-up step that frames the video to it (`guideFrame`).
4. `npx grok-bdd compile` then `npx grok-bdd run generated/guides/<slug>.test.ts` until green; a
   step that needs a wait or a pixel is a missing signal in the core, not a `sleep`.

## 2. Film it

```
cd public/packages/<Pkg>/bdd
MSYS_NO_PATHCONV=1 npx grok-bdd guide features/guides/<slug>.feature [--gif]
```

Watch `guides/<feature slug>/<scenario slug>/guide.mp4/gif` once, and open `audit.png` beside it:
every press as the video shows it, next to the same picture with magenta ticks aimed at where the
press landed. The dot, the ring and the pointer's tip sit between the ticks, on the element the
step lit; the renderer measures both (`audit.json`: `mark_off`, `tip_off` in video pixels, and
prints any press more than 2 px off or outside the lit element). What to fix and where:

- the pointer goes to the wrong place → the step located a scope, not the element: make the
  phrase name the element (`"Open local file" icon inside browse toolbar`), or add the name to
  the core;
- a press is `NOT ON THE LIT ELEMENT` in `audit.png` → the step lit one element and pressed
  another (a stale box, a scope): the same fix — the element the gesture presses is the one it
  names;
- a dialog or balloon is missing from the "after" picture → `--settle 1000`;
- a step reads badly in the caption → reword the step; a verb the caption leaves in the third
  person goes into the `VERBS` map of `src/runtime/guide.ts`;
- a step is missing from the video → it changed nothing on the page (a listener, a server-side
  check) and touched nothing; a step that opens a table or a panel is shown by itself;
- a menu step jumps from the closed menu to the result → the runtime path it went through does
  not report its stops (`guide.hop`); a step that went through the API where a person clicks
  → find or add the named element and make it a `When` (a view switch is `user clicks on X tab`).
- a test feature filmed as a guide fails on a context-panel section or a viewer reading → the
  full shell is not the simple one the test was written in; film a guide feature instead.

Typing shows as it goes in, a drag in several pictures, a wheel turn notch by notch, and a step ends on
a cut to its result. A window narrower than the 1600 px video (`BDD_GUIDE_VIEWPORT`) is filmed at twice
its density. With `--gif` the run also writes `guide-thumb.png` (the result) and `guide-sheet.png`
(twelve frames of the GIF, for a review at a glance); `--fps`, `--colors`, `--size` and `--hold` set
the GIF's rate, palette, width and last-frame hold. The command fails on a press drawn away from the
pointer — a rendering error to fix, not a take to send.

## 3. Send it, and keep it

- Reply with `guide.mp4` (or `guide.gif` where a video does not embed) and the body of
  `steps.md`, whose pictures are the `step-NN.png` next to it.
- Add a line to the package's `bdd/features/guides/INDEX.md`: the question, the feature, the date,
  the help page it illustrates if any. Create the file with a one-line header if absent.
- Leave the feature and its generated spec uncommitted for the lead, with the list of core or
  library changes it needed.

## 4. Help GIFs

A help page shows a GIF about 800 px wide, beside the steps it illustrates (the `write-help` skill
decides where one goes). The GIF is a guide feature filmed in help mode: an 800×500 window at twice
its density (a 1600×1000 GIF, 256 colours, 12 fps), the bare shell (no sidebar, Browse, Toolbox,
Context Panel or Help unless a step opens one), the page as a person sees it, no title card, caption
strip or dimming, no zoom, and no set-up — every `Given` is left out, and the GIF starts a second
before the first action and holds its result two seconds.

1. **Write the feature** in `features/guides/` of the package that owns the area, tagged
   `@guide @help:<page dir>` (the dir under `public/help/`), named by what it shows:
   - one action and its result; set, then undone where the flow allows (a filter applied and reset,
     a column hidden and shown), ending where it started so the loop does not jump;
   - cause before effect: nothing changes on screen before the gesture that changes it;
   - demo data the reader can open (`demog`, a package's sample table) — no fixture with `{run}` or
     `{time}` in a name the GIF shows;
   - the whole viewer in the window, the grid only when the result is in the data (a selection, a
     filter) with the status bar's Selected and Filtered in view;
   - settings changed from the viewer's own context menu where it has one; a step that looks broken
     on video is left out of the feature rather than filmed.
   A replaced GIF repeats the scenario of the one it replaces: watch the old one frame by frame first.
2. **Check it on the local stand**: `npx grok-bdd run generated/guides/<slug>.test.ts`, then
   `npx grok-bdd guide features/guides/<slug>.feature --help-mode` and review `guide-sheet.png`,
   `guide-thumb.png` and single frames at 800 px: nothing cut off, labels readable, cause before effect.
3. **Film it on the stand the help is reviewed on**: `npx grok-bdd guide --help-pages` films every
   `@help:` feature of the package and puts each `<scenario slug>.gif` and `-thumb.png` into
   `public/help/<page dir>/img/`. A GIF the page already has is not overwritten: the new one is
   `<slug>.new.gif`, swapped in after the reviewer approves it; the old file goes only when nothing
   else references it. `public/help` is the datagrok-ai/help submodule: the GIFs are committed there,
   the feature to `public`.
4. **Leave the stand as it was.** What a feature makes goes at feature end (the library's rule), and a
   killed run's leftovers are swept by the fixture's own name family — so on a shared stand a fixture
   carries a name no person would use (the `BDD-` prefix). A step never saves a shared setting (stop
   before **Save and apply**).
