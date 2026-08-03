# Datagrok Browse — project context (handoff)

This document accompanies `browse_manual_tests.md`. It contains all the context needed to continue the work in a new chat: the goal, the environment, scope decisions, sources, a summary of how Browse is structured, and the plan for the next steps.

---

## 1. Goal and current stage

- We are writing **manual test cases for Datagrok Browse**, covering it through the UI.
- The end goal is **Playwright automated tests**. The manual test suite is deliberately structured for easy porting to code.
- Current stage: the manual suite is ready (`browse_manual_tests.md`, 89 cases, 17 sections). Next — capture selectors on dev and start porting to Playwright.

## 2. Environment

- **Target instance:** https://dev.datagrok.ai/
- Tests run as a logged-in user with access rights to demo data.
- On dev there are **demo files** (Files > Demo), a **demo DB** (Databases), and **demo dashboards** (Dashboards) — the preconditions rely on them.

## 3. Test case format

- Markdown, prose descriptions (not tables).
- Each case: **Goal → Type → Preconditions → Steps (action + expected reaction) → Final check → Postconditions/cleanup**.
- **ID:** `Browse-<Area>-<NN>` — the format is compatible with Datagrok Test Track (cf. `Scripts-Edit-6`).
- **Type:** `smoke` / `functional` / `regression` / `negative` — maps to Playwright tags (`@smoke`, etc.).
- Steps are atomic: action → `await` call, expectation → `expect`. Preconditions → `beforeEach`/fixture, postconditions → `afterEach`/teardown.

## 4. Scope decisions (agreed)

- **We cover the whole Browse tree**: My stuff, Spaces, Apps, Files, Dashboards, Databases, Platform + cross-cutting concerns (navigation, tree, view modes, filters/search, routing).
- **Spaces** — in Browse only **smoke** (clicking a node without console errors). Full Spaces coverage (rename, moving, permission inheritance, data loss) goes into a **separate Spaces test suite**, not here.
- **Model Hub / Compute** — only **smoke** in Browse (historical failure points).
- **Global search, Context Panel controls, Split view / viewing modes** — added as base sections (15–17).

## 5. Where test cases live in Datagrok

- Datagrok has a built-in tool, **Test Track** (DevTools → Test Track), for manual test cases. ref: GROK-14402 (DevTools: Test Track: Manual testing).
- The team is moving toward replacing tickets with full-fledged manual test cases. ref: GROK-15857 (Replace tickets in the test track with manual test cases) — our work fits into this approach.
- Cases in Test Track are named using the pattern `<Area>-<Action>-<number>` — we use a compatible format.

## 6. Sources

- Documentation: `https://datagrok.ai/help/datagrok/navigation/` (Navigation / Browse section).
- Jira: instance `reddata.atlassian.net`, project **GROK** (Datagrok). cloudId: `42e2928b-7dc3-47f5-8ab7-52961685645a`.
  - Jira access: JQL search (`searchJiraIssuesUsingJql`) and reading issues work; full-text Rovo `search` returned 403 ("app is not installed") — use JQL.

## 7. Source bugs (regression-case linkage)

A list of Browse bugs referenced by the cases (the `ref:` field). They serve both as regression targets and as a map of "where things historically broke."

- **GROK-16261** — the browse tree should remember its expand/collapse state (Tree-05).
- **GROK-19802** — nested nodes "stick" expanded when the side menu is collapsed (Tree-06).
- **GROK-17922** — clicking an object without access rights: should show a placeholder, not crash (Tree-09, negative).
- **GROK-19848** — Add to Favorites from MyFiles (MyStuff-05, Fav).
- **GROK-20032** — the Apps list took ~40 sec to load (Apps-02).
- **GROK-19638** — errors in the app tooltip/details (Apps-04).
- **GROK-19562** — opening a .h5 file (Files-04).
- **GROK-19844** — a new file in S3/Spaces appeared only after a manual refresh (Files-05).
- **GROK-19847** — the share-folder icon (Files-06).
- **GROK-19662** — a new dashboard appeared only after a refresh (Dash-03).
- **GROK-19934** — duplicated content in the Context Panel when switching dashboards (Dash-04).
- **GROK-16857** — Postgres > CHEMBL > Browse > Summary crashed (DB-06).
- **GROK-19691** — empty/irrelevant properties in the filter panel (Filter-01).
- **GROK-19689** — the CreatedRecently filter for Files (Filter-02).
- **GROK-19688** — the UsedByMe filter for Apps/Dockers (Filter-03).
- **GROK-19692** — the Users filter by group, errors/"balloons" (Filter-04).
- **GROK-19690** — filter and search working together (Filter-05).
- **GROK-19703** — the Models filter by Namespace (Filter-06).
- **Model Hub:** GROK-19965 (double-click crashed, right-click→Run worked), GROK-19740 (click/hover on a model), GROK-19628 (expanding Uncategorized), GROK-17896 / GROK-17664 (opening Model Catalog from Apps), GROK-17770 (opening after saving).

> The list is not exhaustive. To extend it, search Jira: `project = GROK AND summary ~ "Browse" ORDER BY created DESC`. Also look at active Spaces bugs (GROK-19843/54/57/58, 19845) — but those go into the separate Spaces test suite.

## 8. Summary of how Browse is structured (from the docs)

**What Browse is:** a single entry point to everything — projects, queries, files, connections, settings, apps. Opened via the icon on the Sidebar.

**The Browse tree (top-level sections):**
- **My stuff** — Recent, Favorites, Shared with me, My Dashboards, My Files, etc. An object without an explicit project lands here.
- **Spaces** — shared spaces (by permissions).
- **Apps** — installed applications.
- **Files** — My files / Spaces / App Data / Demo.
- **Dashboards** — owned + shared.
- **Databases** — connections → tables → columns → queries; includes a Schema Browser.
- **Platform** — instance configuration (plugins, functions, access control); by permissions.

**Browse Top Menu panel:** Home, Import file, Import text, Refresh tree, Collapse all, Locate current object.

**Browsing mode vs Persistent view:**
- Single click → **unpinned** view (replaced by the next click).
- Double click / any editing / pin → **persistent** view (not replaced).
- On the Sidebar, views of the same type are grouped with a counter badge.

**Tree management:**
- Click / double click / right click (context menu).
- Keyboard: ↑/↓ — navigation, →/← — expand/collapse, Enter/Space — open the item.
- Drag-and-drop — reorganizing the tree and dragging entities into panels.

**Context Panel** (right screen): toggle **F4**; in the header — Back/Forward, Clone and detach, Collapse all / Expand all, Favorites (star). The content follows the selected object.

**Split & viewing modes:** Hamburger → Split right / Split down; modes with a Status Bar — Tabbed / Simple / Presentation.

**Global search:** a search bar on the Home Page; searches by metadata, tags (`#tag`), identifiers, with NLP support.

**Favorites:** available only for **entities** (not for individual files and not for cell values). Adding: context menu, the star in the Context Panel, drag-and-drop into the Favorites Panel.

**Routing:** the URL reflects the state — dashboards, files (`/files/{namespace}.{share}/{path}`), parameterized queries (`/q/{ns}.{conn}.{query}?param=value`).

## 9. Notes for automation (summary)

- **baseURL** = `https://dev.datagrok.ai/`; authentication via `storageState`/global setup.
- **assertNoErrors(page)** — a shared helper: `page.on('pageerror')` + `page.on('console')` (error filter) + a check that no Datagrok error balloon is present. Call it at the end of all "no errors" cases.
- **Creating data for preconditions** (Nav-06, Files-05, MyStuff-04, ModelHub-05) — via the API/a second context, not manually in the UI.
- **Cleanup** — names with a `qa_autotest_<timestamp>` prefix + a shared teardown that deletes everything matching the prefix (a safeguard against failed tests).
- **Asynchronicity** — wait for nodes/elements to appear, not `sleep`.
- **Timings** — move the Apps-02 threshold (5 s) into a config constant.
- **Selectors not captured yet.** On the first dev run, capture stable `data-*`/roles for: the Browse/Favorites icons on the Sidebar; tree nodes and the expand/collapse arrows; the Top Menu; the context menu; the Favorites star and Context Panel controls; tabs and the pin button; the Hamburger (Split); the mode switcher on the Status Bar; the global search bar; the filter panel and its search field.

## 10. Plan for the next steps

1. Walk through the **P1 smoke** cases manually on dev → capture the real selectors → add a "Selectors" block to each case.
2. Align the suite with the team; enter it into Test Track if needed.
3. Port to Playwright in order: first `@smoke` (P1), then `@regression` (bug cases), then the rest.

## 11. Open questions for later

- The exact names of demo entities on dev (file/DB/dashboard) that are stably present — pin them down in fixtures.
- The time threshold for Apps-02 (currently 5 s as an assumption) — confirm with the team.
- The exact composition of the context menu for different entity types (Tree-07) — capture on the first run.
- Whether a whole separate Spaces section is needed (currently agreed — a separate test suite outside this file).
