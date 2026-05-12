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
harness-authored: true
---

# home-page-widgets

## When to use

Your package should surface a small tile on the Datagrok welcome screen
— the view that opens on login, alongside "Recent projects" and "Usage".
Typical phrasings: "show the three most-recently-opened SDF files when
a chemist logs in", "put a KPI counter for active screening campaigns
on the landing page".

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from package root.
- `datagrok-api` imported — `grok`, `DG`, `ui`.
- Familiarity with `DG.Widget` — every tile extends it (`DG-FACT-141`).

## Steps

1. **Author the widget class.** Subclass `DG.Widget`, call
   `super(rootElement)` with the tile's root `HTMLElement`, append
   content to `this.root` (`DG-FACT-141`). Add a `get type(): string`
   override so the widget registry and layout serialization see a
   stable identifier instead of the default `'Unknown'` (`DG-FACT-143`,
   `DG-FACT-DRIFT-HPW-001` — the article example omits this; every
   PowerPack widget has it). Tunable settings (caption, limit) are
   declared with `super.addProperty(name, type, default, options?)`,
   which returns the default so each field initializes in one line
   (`DG-FACT-142`).

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
   method on `PackageFunctions` returning `DG.Widget` and decorate it
   with `@grok.decorators.dashboard({...})` — a specialized decorator
   (not the generic `func`) whose `DashboardOptions` adds `order?: string`
   and `test?: string` on top of `FunctionOptions` (`DG-FACT-139`). The
   role string is `dashboard` lowercase (`DG-FACT-137`); inputs are
   empty, output is `DG.Widget` (`DG-FACT-138`).

   Placement is controlled by `order` — lower numbers (including
   negatives) sort first. PowerPack uses `order: '-1'` for Spotlight,
   `order: '6'` for Community (`DG-FACT-140`). Do NOT declare `order`
   via `super.addProperty('order', ...)` inside the widget — that is a
   per-widget runtime field, NOT the placement control (older article
   draft conflated the two; the rewrite removed that example).

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

3. **Build, and verify the auto-emitted wrapper in `package.g.ts`.**
   `npm run build` runs the function-generator plugin, which scans
   `PackageFunctions` for decorated methods and regenerates
   `src/package.g.ts` with one wrapper per function. Do NOT hand-edit
   `package.g.ts`; it is overwritten on every build.

   ```bash
   npm install && npm run build
   ```

   Expected: `npm run build` exits 0; the regenerated `src/package.g.ts`
   gains a wrapper of the shape

   ```
   //name: Recent SDFs
   //output: widget result
   //meta.role: dashboard
   //meta.order: 5
   export function recentFilesWidget() : any {
     return PackageFunctions.recentFilesWidget();
   }
   ```

   The MUST-have lines are `//meta.role: dashboard` and
   `//output: widget result` (`DG-FACT-138`). Compare
   `packages/PowerPack/src/package.g.ts:26-32` (Community) or
   `:17-24` (Spotlight, with `order: -1`).

4. **Publish.** The `dashboard` role IS the registration — the platform
   picks the function up from the `package.g.ts` wrapper at load; no
   explicit `grok.functions.register(...)` is needed.

   ```bash
   grok publish <host>   # add --release once stable
   ```

   Expected: exits 0 and reports the package version pushed. After a
   Datagrok reload, the welcome screen shows the new tile in the slot
   dictated by `order`.

## Common failure modes

- **Tile never appears on the welcome screen.** Role token missing or
  misspelled. `src/package.g.ts` MUST contain `//meta.role: dashboard`
  (lowercase) and `//output: widget result` (`DG-FACT-137`,
  `DG-FACT-138`). `Dashboard`/`dashboards` won't match.
- **Tile appears but in the wrong slot.** `order` set on the wrong
  surface. Placement is the function-level
  `@grok.decorators.dashboard({order: '5'})`, NOT
  `super.addProperty('order', ...)` (`DG-FACT-140`). Verify
  `//meta.order: 5` in `package.g.ts`.
- **Widget serializes as `type: 'Unknown'`.** Class missing the
  `get type(): string` override (`DG-FACT-143`, `DG-FACT-DRIFT-HPW-001`).
  Add `get type(): string { return '<ClassName>'; }` at the top of
  the class body.
- **Settings don't appear in the gear panel.** `super.addProperty(...)`
  called outside the constructor, or its return wasn't assigned
  (`DG-FACT-142`). Canonical:
  `this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Recent files')`.

## Verification

- `npm run build` and `grok publish <host>` both exit 0.
- `src/package.g.ts` carries `//meta.role: dashboard`,
  `//output: widget result`, and `//meta.order: <n>` matching the
  decorator (compare `packages/PowerPack/src/package.g.ts:17-32`).
- In Datagrok, the welcome screen shows the tile in the slot set by
  `order`; right-click → *Settings* exposes `caption` and other
  `addProperty` fields.

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
