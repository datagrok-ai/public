---
name: home-page-widgets
version: 0.1.0
description: |
  Register a function in a Datagrok package that the platform renders
  as a tile on the welcome screen the user lands on after logging in —
  alongside the built-in "Recent projects" and "Usage" panels. For
  plugin authors who want a domain dashboard (recent SDF files, top
  hits, KPI badge) visible before the user navigates anywhere.
  Produces a `static` method on `PackageFunctions` decorated with
  `@grok.decorators.dashboard({...})` plus a `DG.Widget` subclass for
  the tile body.
  Use when asked to "show a tile on the welcome screen when users log in",
  "put a custom panel on the start page next to Recent projects",
  or "add a landing-page dashboard tile for our plugin".
triggers:
  - tile on the welcome screen
  - custom panel on the start page
  - landing screen dashboard tile
  - widget next to recent projects
  - panel before users navigate anywhere
  - home page tile for plugin
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# home-page-widgets

## When to use

Your package should surface a small tile on the Datagrok welcome screen
— the view that opens on login, alongside "Recent projects" and "Usage".
Typical phrasings: "show the three most-recently-opened SDF files when
a chemist logs in", "put a KPI counter for active screening campaigns
on the landing page".

## Prerequisites

- Familiarity with `DG.Widget` — every tile extends it (`DG-FACT-141`).

## Steps

1. **Author the widget class.** Subclass `DG.Widget` and call
   `super(rootElement)` (see `DG-FACT-141`). Always override
   `get type(): string` — the default is `'Unknown'` (see `DG-FACT-143`,
   `DG-FACT-DRIFT-HPW-001`). Use `super.addProperty(name, type, default, options?)`
   for tunable settings — it registers AND returns the default (see
   `DG-FACT-142`).

   ```typescript
   // src/widgets/recent-files-widget.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   import * as ui from 'datagrok-api/ui';

   export class RecentFilesWidget extends DG.Widget {
     get type(): string { return 'RecentFilesWidget'; }

     caption: string;
     limit: number;

     constructor() {
       super(ui.panel([], 'recent-files-widget'));
       this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Recent files');
       this.limit = super.addProperty('limit', DG.TYPE.INT, 5);
       // fetch + append to this.root
     }
   }
   ```

   Expected: file type-checks; `DG.Widget`, `DG.TYPE`, and `ui.panel`
   resolve from the imports above.

2. **Register as a dashboard function.** Declare a zero-arg `static`
   method on `PackageFunctions` returning `DG.Widget`, decorated with
   `@grok.decorators.dashboard({...})` — a specialized decorator with
   `order?` / `test?` on top of `FunctionOptions` (see `DG-FACT-137`,
   `DG-FACT-138`, `DG-FACT-139`).

   Placement is the decorator-level `order` string — lower sorts first;
   PowerPack uses `'-1'` (Spotlight) / `'6'` (Community) (see
   `DG-FACT-140`). Do NOT use `super.addProperty('order', ...)` for this.

   ```typescript
   // src/package.ts
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';
   import {RecentFilesWidget} from './widgets/recent-files-widget';

   export const _package = new DG.Package();

   export class PackageFunctions {
     @grok.decorators.dashboard({order: '5', name: 'Recent SDFs'})
     static recentFilesWidget(): DG.Widget {
       return new RecentFilesWidget();
     }
   }
   ```

   Expected: file type-checks; `@grok.decorators.dashboard` resolves
   from `datagrok-api/grok`.

3. **Build and verify the auto-emitted wrapper.**
   ```bash
   npm install && npm run build
   ```
   The regenerated `src/package.g.ts` must include both
   `//meta.role: dashboard` and `//output: widget result` (see
   `DG-FACT-138`). Don't hand-edit `package.g.ts` — it is overwritten.
   Compare `packages/PowerPack/src/package.g.ts:17-24` and `26-32`.

4. **Publish.** The `dashboard` role IS the registration — no
   explicit `grok.functions.register(...)` is needed.

   ```bash
   grok publish <host>   # add --release once stable
   ```
   After reload the tile appears in the slot dictated by `order`.

## Common failure modes

- Tile never appears — `package.g.ts` is missing `//meta.role: dashboard`
  (lowercase) or `//output: widget result` (`DG-FACT-137`, `DG-FACT-138`).
- Tile in wrong slot — `order` set via `addProperty` instead of the
  decorator (`DG-FACT-140`); verify `//meta.order: 5` in `package.g.ts`.
- Widget serializes as `type: 'Unknown'` — missing `get type()` override
  (`DG-FACT-143`).
- Settings absent from gear panel — `super.addProperty(...)` called
  outside the constructor or its return not assigned (`DG-FACT-142`).

## See also

- Source article: `help/develop/how-to/packages/home-page-widgets.md`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-137/138/139/140/141/142/143` and drift
  `DG-FACT-DRIFT-HPW-001`.
- Reference packages: `packages/PowerPack/src/package.ts:138-155` +
  `widgets/community-widget.ts:6-28` +
  `widgets/recent-projects-widget.ts:8-19` (Spotlight, Community,
  Recent projects); `packages/UsageAnalysis/src/package.ts:397-415`
  (Usage / Reports, gated via `meta.canView`).
- Related skills: `create-package`; `column-tooltip` (same
  `DG.Widget` shape, different surface).
