---
name: js-api
description: Navigate the Datagrok JS API — pick the right entry point (`grok` / `ui` / `DG`) for a task, register functions, and subscribe to events
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# js-api

## When to use

You're writing TypeScript inside a Datagrok package (or an ad-hoc
`Functions | Scripts | New JavaScript Script`) and need to know which of
the three JS API entry points exposes the thing you want, what the import
line looks like, and how to verify the surface is reachable. Triggers:
"how do I import the JS API", "where does `grok.shell` live", "register a
function from JS", "subscribe to a platform event", "dock a custom div".

## Prerequisites

- A package scaffold (`grok create <Name>`); `datagrok-api` already in
 `dependencies` (the `grok create` template wires this).
- A logged-in platform tab to verify with the Console (`Tools → Console`)
 and Inspector (`Alt+I`).

## Steps

1. **Add the three canonical imports.**
 Every package file that touches the API uses the same three lines.
 `grok` is the discoverability root; `ui` builds the DOM; `DG` is the
 raw class namespace (knowledge `DG-FACT-224`). Confirm with
 `packages/Tutorials/src/package.ts`, `packages/PowerPack/src/package.ts`,
 `packages/Chem/src/package.ts`.
 ```typescript
 import * as grok from 'datagrok-api/grok';
 import * as ui from 'datagrok-api/ui';
 import * as DG from 'datagrok-api/dg';
 ```
 Expected: TypeScript resolves all three; `webpack` build succeeds with
 no `Cannot find module 'datagrok-api/*'` errors.

2. **Pick the entry point from the task.**
 Match the intent against the table; the namespace inventory is
 exhaustive (knowledge `DG-FACT-225`). When in doubt, the article's
 "Grok" section (lines 31-47) recommends starting at `grok.` and
 following IntelliSense; drop to `DG.<Class>` only when you need to
 instantiate a class directly.

 | Intent | Entry point | Knowledge / skill |
 |---|---|---|
 | Load / build a DataFrame, run a query | `grok.data.*` | `access-data`, `DG-FACT-008` |
 | Open / dock a view, set current view | `grok.shell.*` | `routing`, `DG-FACT-018/019` |
 | Register a function (panel, viewer, widget…) | `grok.functions.register({...})` | `DG-FACT-227` |
 | Server entities (projects, files, packages, users) | `grok.dapi.*` | `access-data`, `DG-FACT-026/032` |
 | Per-user persistent KV (`<5000` chars/value) | `grok.userSettings.*` | `user-settings-storage`, `DG-FACT-046..049` |
 | Subscribe to platform / DataFrame events | `grok.events.*`, `frame.on*` | `DG-FACT-226` |
 | Cheminformatics (similarity, descriptors, SAR) | `grok.chem.*` | (chem how-tos) |
 | Machine learning (clustering, PCA, predict) | `grok.ml.*` | (ml how-tos) |
 | Build a dialog / button / div / accordion | `ui.*` | `custom-views`, `DG-FACT-213..217` |
 | Construct a class (`DataFrame`, `Widget`, `ViewLayout`, …) | `DG.<Class>` | per-class refs |
 | Custom semantic-type renderer | `DG.ObjectHandler` (NOT `EntityMeta`) | `DG-FACT-059/060` |
 | Dock pose constants (`LEFT`, `RIGHT`, `TOP`, `DOWN`, `FILL`) | `DG.DOCK_TYPE.*` | `js-api/src/const.ts:775-782` |

 Expected: one entry point picked. If two seem to fit, prefer the one
 that exposes an instance helper over a raw class instantiation.

3. **Register a JS function from the package.**
 `grok.functions.register({signature, run, tags?, isAsync?, namespace?,
 options?})` — the signature uses Datagrok's type vocabulary (`String`,
 `int`, `List<String>`, `Widget`, `DataFrame`), NOT raw TypeScript types
 (knowledge `DG-FACT-227`). Set `isAsync: true` whenever `run` is async.
 ```typescript
 grok.functions.register({
 signature: 'String jsConcat(int foo, int bar)',
 run: (foo: number, bar: number) => `${foo}_${bar}`,
 });

 grok.functions.register({
 signature: 'Widget jsClock',
 tags: 'Widgets',
 run: => {
 const e = document.createElement('div');
 setInterval( => { e.innerText = new Date.toTimeString; }, 1000);
 return new DG.Widget(e);
 },
 });
 ```
 Expected: typing `jsConcat(42, 33)` in the platform Console returns
 `'42_33'`; "Widgets" dashboard exposes `jsClock`.

4. **Subscribe to a platform or DataFrame event.**
 Event getters return RxJS `Observable<T>` (knowledge `DG-FACT-226`);
 keep the returned `Subscription` so you can `unsubscribe` from
 `init` / `autostart` cleanup. Per-instance event names live on the
 instance (`frame.onValuesChanged`, `view.onClosed`, etc.).
 ```typescript
 const sub = grok.events.onCurrentProjectChanged.subscribe((_) =>
 grok.shell.info(`Project: ${grok.shell.project.name}`));
 // …later:
 sub.unsubscribe;
 ```
 Expected: switching the current project pops a balloon. To find an
 event name, open Inspector (`Alt+I`) → "Client Log" tab → perform the
 action; the panel shows the event id and an auto-generated handler
 snippet (knowledge `DG-FACT-228`).

5. **Dock an arbitrary visual element.**
 `grok.shell.dockManager.dock(element, dockType, refNode?, title?)`
 accepts either a `DG.DOCK_TYPE` enum value or its string equivalent
 (`'left' | 'right' | 'up' | 'down' | 'fill'` — note `TOP === 'up'` and
 `DOWN === 'down'`).
 ```typescript
 const panel = ui.divText('hello from JS');
 grok.shell.dockManager.dock(panel, DG.DOCK_TYPE.RIGHT, null, 'JS');
 ```
 Expected: a docked pane labelled "JS" appears at the right edge.

6. **Smoke-test the API surface from a fresh package.**
 In `src/package.ts`, add an `init` function and `webpack && grok publish
 <alias>`. Reload the platform; `grok.shell.info(...)` confirms the
 bundle loaded and all three entry points resolved.
 ```typescript
 //name: <Pkg>Init
 //tags: init
 export async function init {
 grok.shell.info(
 `JS API ok — grok=${typeof grok}, ui=${typeof ui}, DG=${typeof DG}`);
 }
 ```
 Expected: a balloon reading `JS API ok — grok=object, ui=object,
 DG=object`. Anything else (red error toast, `Cannot read properties of
 undefined`) means an import or version drift.

## Common failure modes

- **`grok.testData(...)` throws "is not a function".** The article's
 line-81 snippet drops `data.`; canonical form is `grok.data.testData(name,
 rows?, cols?)`. Fix: prefix with `grok.data.`.
- **`extends DG.EntityMeta` fails to compile.** That class does not
 exist — the article's "User-defined types" section (lines 210-217) is
 out of date. The actual base is `DG.ObjectHandler` registered via
 `DG.ObjectHandler.register(handler)`. Fix: subclass `ObjectHandler`; see the
 `register-identifiers` skill for the full pattern.
- **Article's IntelliSense list `[shell, chem, data, data, ml]` is wrong.**
 `data` is listed twice and the list is incomplete; the real inventory
 is `functions`, `events`, `dapi`, `shell`, `settings`, `data`,
 `userSettings`, `ai`, `log`, `chem`, `ml`, `decorators` (knowledge
 `DG-FACT-225`). Fix: trust the inventory, not the article.
- **`grok.dapi.fetchProxy(...)` example URL `restcountries.eu` returns
 an empty array.** That domain has been off-air for years; the project
 moved to `restcountries.com`. Fix: swap the
 URL when copy-pasting line 167 of the article.
- **Event subscription leaks across hot reloads.** `subscribe(...)`
 returns a `Subscription`; if you discard it, repeated `init` calls
 during dev mode multiply handlers. Fix: store the subscription on the
 package object (or in a module-level array) and `unsubscribe` when
 the package re-inits.

## See also

- Source articles:
 - `help/develop/packages/js-api.md` (this overview).
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
 `DG-FACT-008` (testData enum), `DG-FACT-018/019` (shell view setters),
 `DG-FACT-026` (fetchProxy), `DG-FACT-046..049` (userSettings),
 `DG-FACT-059/060` (ObjectHandler), `DG-FACT-213..217` (custom views),
 `DG-FACT-224` (three entry points), `DG-FACT-225` (`grok.*` inventory),
 `DG-FACT-226` (events as RxJS Observables), `DG-FACT-227` (function
 signature shape), `DG-FACT-228` (Inspector for event discovery),
 (testData prefix), (dead URL).
- Related skills: `extensions` (function-role tags),
 `access-data` (`grok.data.query`, `grok.dapi.*`),
 `user-settings-storage` (`grok.userSettings.*`),
 `register-identifiers` (`DG.ObjectHandler`),
 `custom-views` (`grok.shell.dockManager`, `ui.*`),
 `manipulate-viewers` (`DG.Viewer` instantiation).
