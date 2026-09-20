# Recipe: a domain app as a spec (configuration)

The second tier of a domain (EMS) app: the same page the [zero-code tier](crud-app.md) renders,
laid out as a `dg-ui/1` spec — a `u2-domain-source` and the `u2-domain-*` controls bound to it.
No TypeScript beyond mounting the spec; the designer edits it; an LLM can emit it. Reference:
**Stockroom**'s `src/app.spec.json` (`packages/Stockroom`), abridged here.

## The shape

```json
{"$schema": "dg-ui/1",
 "root": {"tag": "u2-panel", "name": "app", "children": [
   {"tag": "u2-domain-source", "name": "substances", "props": {"table": "stockroom.substance", "pageSize": 50}},
   {"tag": "u2-splitter", "props": {"direction": "horizontal", "sizes": [0.35, 0.65]}, "children": [
     {"tag": "u2-panel", "children": [
       {"tag": "u2-domain-search",  "bind": {"source": "$.substances.source"}, "props": {"placeholder": "Name or CAS"}},
       {"tag": "u2-domain-filters", "bind": {"source": "$.substances.source"}},
       {"tag": "u2-domain-list",    "bind": {"source": "$.substances.source"}, "props": {"mode": "cards"}}]},
     {"tag": "u2-panel", "children": [
       {"tag": "u2-domain-form",     "bind": {"source": "$.substances.source"}},
       {"tag": "u2-domain-children", "bind": {"source": "$.substances.source"},
        "props": {"tables": ["container", "sds_document"], "mode": "list"}},
       {"tag": "u2-domain-history",  "bind": {"source": "$.substances.source"}},
       {"tag": "u2-button", "props": {"text": "New"},                    "on": {"click": "cmd:substances.newRow"}},
       {"tag": "u2-button", "props": {"text": "Save", "primary": true}, "on": {"click": "cmd:substances.save"}},
       {"tag": "u2-button", "props": {"text": "Discard"},                "on": {"click": "cmd:substances.discard"}}]}]}]}}
```

Mounting it is the [shell-integration](shell-integration.md) shape — the spec instance is a
control, so it is the view's content; the designer opens the same file for editing:

```ts
import {renderSpec} from '@datagrok-libraries/u2';
import {appView, designerView, registerPlatformComponents} from '@datagrok-libraries/u2/src/dg/index.js';
import spec from './app.spec.json';                              // tsconfig: resolveJsonModule

registerPlatformComponents();                                    // once per package: every u2-* tag
export function stockroomApp(): DG.ViewBase {
  return appView({name: 'Stockroom', content: renderSpec(spec)});
}
export function stockroomDesigner(): DG.ViewBase {
  return designerView(spec, {name: 'Stockroom (design)'});
}
```

## The pieces

| Tag | What it is |
|---|---|
| `u2-domain-source` | the data: `table`, `query` and `search` (bindable), `pageSize`, `defaults` (a parent's id on a child table), `empty`, `draft` (a create form's source). Functions for `cmd:`: `save`, `discard`, `newRow` |
| `bind: {source: "$.substances.source"}` | how every control takes its source — whole, through the `source` step; the control inherits the source's access. `$.substances` alone is the rows, `$.substances.currentRow.<col>` a value a plain input may bind to |
| `u2-domain-list` (`mode`, `itemHeight`, `empty`) | the rows, one line or a card each; selecting a row makes it the source's current row; Open / Delete / the table's actions per row |
| `u2-domain-form` (`include`, `exclude`, `layout`) | the current row's editors under the row's access; readonly is text |
| `u2-domain-grid` | the platform grid over the same rows, editable in place — for many columns or bulk edits |
| `u2-domain-search` (`placeholder`, `debounceMs`) | writes the source's `search`: the schema's searchable columns |
| `u2-domain-filters` (`mode: query \| builder`, `placeholder`) | writes the source's `query`, two-way; a change while dirty goes through the unsaved-changes gate |
| `u2-domain-children` (`tables`, `mode: grid \| list`) | one tab per referring table over the current row, every child source in the same session |
| `u2-domain-history` | the current row's audit trail, newest first |
| `u2-domain-pick` (`table`, `filter`, `value`) | a picker over another table — a generated form uses it by itself for a `ref` column |
| `cmd:substances.save` / `.discard` / `.newRow` | a button's `on.click`: the named source's functions. Save goes through the session, so it lands every source of the spec |

Every tag's registry entry carries its props, defaults and a `usage` line — what the designer's
palette and an LLM read. `u2-domain-source` is in `manifest.json` (the platform-free registry);
the `u2-domain-*` controls are platform-side and arrive with `registerPlatformComponents`
(`src/dg/domain/registrations.ts`).

## The session

Every `u2-domain-source` of one spec instance joins the instance's `SharedSession` (the ambient
session `renderSpec` builds under). Two sources on one page — a `projects` source and an `issues`
source with `defaults: {project_id: ...}` — are one unit of work: one Save button, wired to either
source's `save`, writes both as one transaction; the parent may still be a draft (the child holds
its `~new:` id, the server orders the inserts); Discard drops both; the summary counts both. A
`u2-domain-children` adds its child sources to the parent's session by itself.

## What the designer edits

- Every prop in the tables above, live, from the property panel; layout (`u2-splitter` direction
  and sizes, the `u2-panel`s) and appearance (the shared CSS props) the same way.
- Composition: drag a `u2-domain-*` tag from the palette onto the canvas, bind its `source` to
  the page's source by picking it; a second `u2-domain-source` for a second table.
- Design/Run toggle in the ribbon: the running page is the same spec over the real backend.
- Phase 4 stores this spec on the table itself (designer-edited, schema-admin owned), so an app
  has a default page without a file — nothing written today may make that harder: every domain
  control is spec-registered from day one (GOAL ruling 8).

## When to reach for code

An action on a row, a validator beyond the schema's, a card of your own, presets or shortcuts —
the [code tier](custom-app.md): the spec keeps the layout, the code fills the table's registries.

## Anti-patterns

- Binding `$.substances` or `$.substances.currentRow` to a domain control — the control needs
  the source: `$.substances.source`.
- A Save button per control, or `cmd:substances.save` beside `cmd:containers.save` — one session,
  one Save; the second button is a duplicate of the first.
- Two `u2-domain-source`s over the same table on one page to get "the list" and "the form" — one
  source; the form follows the list's selection through `currentRow`.
- A `u2-text-input` bound to `$.substances.search` with a debounce written in a script — that is
  `u2-domain-search`.
