---
name: write-help
description: Write, update, and verify Datagrok help pages (help/, Docusaurus), including GIFs and release updates. Use for user documentation of a platform feature or a Jira ticket, for bringing the help up to date with a release, and for recording or replacing help GIFs.
when-to-use: A feature or a ticket (GROK-12345) needs user documentation; a release shipped and the help must describe it; a help GIF is missing, outdated, or unreadable.
argument-hint: "<topic, page, ticket key (GROK-12345), or release version>"
effort: medium
---

# Writing Datagrok help pages

Mechanics (front matter, headings, links, images, admonitions, build) are in `help/CLAUDE.md`;
language, capitalization, and punctuation in `help/develop/help-pages/writing-style.md` and
`word-list.md`. This skill covers what to write, how to prove it, how to shape it, and how to hand
it over. A GIF on a help page follows `gifs.md` even when it is filmed with `grok-bdd guide`: the
`bdd-answer` skill's full shell and captions are for answers sent to users.

Companion files in this skill's folder:

| File | Read when |
|---|---|
| `release.md` | Updating the help for a release: finding what shipped, help vs community post |
| `gifs.md` | Recording or replacing a GIF: requirements, tools, process, pitfalls |
| `tools/linkcheck.js` | Before every handover: relative links and anchors |
| `tools/rec-lib.mjs`, `tools/examples/` | Recording a GIF when `grok-bdd guide` can't meet `gifs.md` |

## 1. Process

1. **Find what exists before writing.** Grep the help for the topic and list the pages and
   sections that cover it. Extend an existing page rather than create one. Don't create overview
   or workflow pages that restate other pages: the sidebar already links them.
2. **Propose the structure first**: which pages and sections change, and where each fact will
   live. Wait for approval.
3. **A question from the reviewer is a request for your opinion, not a command.** Say whether you
   agree and why, propose the exact text, and wait for the decision.
4. **Work in small, reviewable edits** on a branch cut from the current master. `help/` is a
   submodule of the `datagrok-ai/help` repository: branch, commit, and open the PR there, never
   commit the `help` pointer in `public` (a bot bumps it). Commit and push only when asked, with
   plain commit messages.
5. **Finish with a reading list**: for each topic, the page and anchor that now covers it.

### Starting from a ticket

When the argument is a ticket key, learn the feature from the ticket and the code before deciding
what the help needs:

1. **Read the ticket**: summary, description, comments, status, fix version, linked tickets, and
   PRs (Jira, through the Atlassian tools if they are connected).
2. **Find the change**: `git log origin/master --grep=<KEY>` in core, `public`, and `help`. Read
   the diffs, then the current code on master: commits after the ticket may have changed it.
3. **Decide whether it belongs in the help**: user-facing, finished, and released (see the table in
   `release.md`). Internal changes, unfinished work, and unreleased plugin versions get no help.
   Report the decision with the evidence.
4. **Find the place**: grep the help for the area and its UI labels. Prefer updating the section
   that already covers the area. A help page that now says something wrong comes first.
5. **Propose, verify, write** as in sections 1–6. Put the ticket key at the start of the commit
   message.
6. **Suggest a community post** when the change is new and noticeable to users, changes existing
   results, or requires an action: one line in the reader's words, with the help page and anchor
   if there is one.

## 2. Verify every claim

- **UI labels, menu paths, options, and defaults come from the code**, quoted exactly as rendered.
  Read the code from the current master, not from a local working tree that may be behind.
  - A `@Prop` camelCase field renders as words (`showMovingAverageLine` → **Show Moving Average
    Line**). An explicit `name:` or `caption:` wins. The generated property tables use the field
    name and ignore `name:`, so don't copy labels from them.
  - Menu paths: the package function's `top-menu`, the viewer's context-menu code. When a page
    mentions one menu of a package, check all its paths on the page: stale paths travel together.
- **Behaviour is checked in the running product** before it is written: what a click does, what a
  dialog opens, whether a setting persists, what the result looks like. Permissions, sharing, and
  saving are checked with two accounts (author and recipient). Use a small headless script (see
  `gifs.md`) or the browser, and record what you observed.
- **Numbers come from the data**, not from assumptions.
- **Released, finished, user-facing.** Plugin changes under `## v.next` in a `CHANGELOG.md` are
  not released. Ask before documenting a feature whose ticket is open; its developer may own the
  docs.
- **Never write an unverified claim as fact.** Keep a list of claims checked only in the code or
  not checked at all, and show it. If a check is impossible, drop the claim.
- **Bugs found on the way** go to the reviewer as reproduction steps, not into the help. File a
  ticket only when the reviewer asks.
- **Don't document legacy functionality** or link to it. Ask when unsure.
- **Clean up the stand**: delete every entity you created, by id, and verify it's gone. Never save
  shared settings to try something out.

## 3. Terms and style

- Reuse the terms the help already uses (`word-list.md`): grep before introducing one.
- Bold is for UI elements the reader clicks or reads: buttons, menu items, fields, options.
  Concepts get plain text, or italics on first definition.
- Headings are 1–3 words in sentence case. Actions use gerunds ("Saving a dashboard"),
  consequences use "When …" ("When the source changes").
- Active voice, second person, one idea per sentence. Split sentences longer than about 25
  words. No semicolons. No motivational sentences, no procedures for flows the UI makes obvious,
  no "best practices" lists unless every item is verified and actionable. When asked to rewrite,
  shorten: remove what the reader already knows.

## 4. Shaping the content

- **Visible text states the rule and the recommendation.** Supporting detail (comparison tables,
  long option lists, step-by-step procedures, worked examples) goes under
  `<details><summary>How to use</summary>…</details>` or another descriptive summary.
- **Short sections are prose, not bullets.** Two links in a list become two sentences under a
  heading that says what they are.
- **Every section has a heading.** Text after an image without one reads as an orphan.
- **Admonitions:**
  - `:::note` for one exception to the rule just stated or a fact the reader would miss. At most
    two per section.
  - `:::tip` for advice the reader may skip.
  - `:::warning` and `:::caution` only for irreversible actions or data loss.
  - `:::note developers` for links to `develop/` pages, JS API, or samples, phrased "You can
    [do X](link)".
- **Link, don't repeat.** If a page already explains something, write one sentence and link to
  it. Shared behaviour lives on one page (for example, regression and moving average lines on the
  scatterplot page), and other pages say what is theirs and link there.
- **One link per concept, at the right place.** Link the first mention and the next page the
  reader needs. Two concepts get two links to their own sections, never one link for both. No
  "see also" to pages one click away in the sidebar.
- **Tables or prose.** A table when the reader scans: options and meanings, formats, comparisons
  by several criteria, more than three rows people look up. Prose when it has two or three rows,
  full sentences in every cell, or only "what" and "why" columns. On concept and how-to pages a
  table is supporting detail under `<details>`; on reference pages it stays open.
- **Keep the order of related sections the same across sibling pages** (for viewers: trend lines,
  then formula lines and annotations). Don't restructure page templates unless asked.
- **Alternative GIFs of one section go into tabs** (for example, regression line and moving
  average in `visualize/viewers/line-chart.md`). `help/CLAUDE.md` reserves tabs for OS or language
  variants: when publishing this skill, add this case there. A page with tabs needs
  `mdx:\n  format: mdx` in the front matter and the `Tabs` / `TabItem` import block.
- **No troubleshooting sections in reference pages.** Non-obvious questions go to
  `datagrok/resources/faq.md`.
- **Renaming a heading changes its anchor.** Grep the whole help (release notes included) and
  `js-api/src` for the old anchor and fix every link.
- **Headings the reviewer gives are used as given**, after checking them against the code. Say
  so if one is inaccurate.

## 5. Images and GIFs

- **Only where motion or layout matters**: a new interaction or a multi-step flow. Static
  settings get none. One image per section at most. If the text says it all, remove the image.
- **Reuse** an existing image by relative path instead of recording a duplicate.
- Files go next to the page in `img/`, named by what they show. Every GIF has a `-thumb.png` of
  a meaningful frame next to it: in-app Context Help shows the thumbnail. Alt text describes the
  action. Delete images you replaced.
- How to record, and what a help GIF must look like: `gifs.md`.

## 6. Check and hand over

1. **Links**: from `help/`, run `node <skill dir>/tools/linkcheck.js` (pages changed on the branch
   plus new ones) or pass the pages. Fix every link you broke. Report pre-existing breakage, and
   fix it only if the reviewer agrees.
2. **Build**: the local server and build commands are in `help/CLAUDE.md`. If the API docs fail to
   build locally (they need built libraries), save this as `docusaurus/docusaurus.nodocs.config.js`
   and delete it before committing:

   ```js
   const config = require('./docusaurus.config.js');
   module.exports = {...config, plugins: config.plugins.filter((p) =>
     p[0] !== 'docusaurus-plugin-typedoc' && !(p[1] && p[1].id === 'api'))};
   ```

   Then, from `docusaurus/`: `npx docusaurus start --port <free port> --config
   docusaurus.nodocs.config.js` to review, and `npx docusaurus build --config
   docusaurus.nodocs.config.js` before merging. The build fails on broken links: with this
   config, links to `/api` are expected; everything else is real. Check which process holds a
   port before stopping anything, and never stop a process you didn't start.
3. **Review table**: page → sections changed → link to the anchor on localhost. Keep a running
   list of approved and pending pages. Re-open a page for review when you change it again.
4. **Secrets**: examples use placeholders (`<password>`, `<token>`). The repository runs
   gitleaks on every push, and `.gitleaksignore` entries are pinned to line numbers: adding lines
   above an ignored example breaks the ignore, so replace such values with placeholders instead.
5. **Report**: what was documented (page → section), what was left out and why, what was verified
   only in the code, bugs found (steps), and what is pending.
