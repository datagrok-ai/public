# PowerPack

Commonly used Datagrok platform enhancements: dashboard widgets, search providers, formula lines editor, viewers gallery, Excel import, and custom inputs.

## Directory Layout

```
src/
  package.ts             # Entry point, PackageFunctions class with all registrations
  package.g.ts           # AUTO-GENERATED (do not edit)
  widgets/               # Dashboard widgets and custom inputs
    activity-dashboard-widget.ts, kpi-widget.ts, web-widget.ts, html-widget.ts, ...
    cron-input.ts         # CronInput (JsInputBase<string>) - cron schedule editor
    cron-utils.ts         # Pure cron parsing/validation/description utilities
  dialogs/               # Dialog UIs (formula-lines, add-new-column)
  search/                # Search providers (entity-search, power-search, help)
  workers/               # Web workers (exceljs-worker)
  tests/                 # Test files
css/
  power-pack.css         # All styles, loaded via "sources" in package.json
```

## Build

```bash
npm run build    # grok api && grok check --soft && webpack
```

## Key Patterns

### Function Registration
All functions are registered via `@grok.decorators` on the `PackageFunctions` class in `package.ts`.

### Dashboard Widgets
Extend `DG.Widget`, register with `@grok.decorators.dashboard()`. See `KpiWidget`, `CommunityWidget`.

### Custom Inputs (Value Editors)
Extend `DG.JsInputBase<T>`, register with:
```typescript
@grok.decorators.func({
  meta: { propertyType: 'string', semType: 'yourSemType', role: 'valueEditor' },
  outputs: [{type: 'object', name: 'result'}],
})
static yourInput(): DG.InputBase { return new YourInput(); }
```
See `CronInput` for a complete example.

### CSS
All styles use prefixed class names (`pp-kpi-`, `pp-cron-`, `power-pack-`). CSS is loaded by the platform via the `"sources"` field in `package.json`.
