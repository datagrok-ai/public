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

1. **Add the three canonical imports** (`DG-FACT-224` — `grok` is the
 discoverability root, `ui` builds DOM, `DG` is the raw class namespace).
 ```typescript
 import * as grok from 'datagrok-api/grok';
 import * as ui from 'datagrok-api/ui';
 import * as DG from 'datagrok-api/dg';
 ```

2. **Pick the entry point from the task** (full inventory at
 `DG-FACT-225`). Start at `grok.` and follow IntelliSense; drop to
 `DG.<Class>` only when you need to instantiate a class directly.

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

3. **Register a JS function from the package** (`DG-FACT-227` —
 signature uses Datagrok types like `String`/`int`/`Widget`, not raw
 TypeScript; set `isAsync: true` for async `run`).
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

4. **Subscribe to a platform or DataFrame event** (`DG-FACT-226` —
 event getters return RxJS `Observable<T>`; keep the `Subscription`
 for `unsubscribe`). Discover event names via Inspector (`Alt+I`) →
 "Client Log" (`DG-FACT-228`).
 ```typescript
 const sub = grok.events.onCurrentProjectChanged.subscribe((_) =>
 grok.shell.info(`Project: ${grok.shell.project.name}`));
 // …later:
 sub.unsubscribe;
 ```

5. **Dock an arbitrary visual element.**
 `grok.shell.dockManager.dock(element, dockType, refNode?, title?)`
 accepts either `DG.DOCK_TYPE` enum or string (`'left'|'right'|'up'|
 'down'|'fill'` — note `TOP === 'up'`).
 ```typescript
 const panel = ui.divText('hello from JS');
 grok.shell.dockManager.dock(panel, DG.DOCK_TYPE.RIGHT, null, 'JS');
 ```

6. **Smoke-test the API surface** — `webpack && grok publish <alias>`,
 then reload and check the balloon.
 ```typescript
 //name: <Pkg>Init
 //tags: init
 export async function init {
 grok.shell.info(
 `JS API ok — grok=${typeof grok}, ui=${typeof ui}, DG=${typeof DG}`);
 }
 ```

## Common failure modes

- `grok.testData(...)` not a function — canonical form is `grok.data.testData(name, rows?, cols?)`.
- `extends DG.EntityMeta` fails to compile — class doesn't exist; subclass `DG.ObjectHandler` (`DG-FACT-059`).
- Article's IntelliSense list is wrong — see `DG-FACT-225` for real inventory.
- `restcountries.eu` returns empty — domain off-air, use `restcountries.com`.
- Event subscription leaks across hot reloads — store the `Subscription` and `unsubscribe` on re-init.

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
