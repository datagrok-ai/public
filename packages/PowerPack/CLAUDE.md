# PowerPack

Commonly used Datagrok platform enhancements: dashboard widgets, search providers, formula lines editor, 
viewers gallery, Excel import, custom inputs, and the spotlight feature.

## Directory Layout

```
src/
  package.ts             # Entry point, PackageFunctions class with all registrations
  package.g.ts           # AUTO-GENERATED (do not edit)
  widgets/               # Dashboard widgets and custom inputs
    spotlight-widget.ts, kpi-widget.ts, web-widget.ts, html-widget.ts, ...
    cron-input.ts         # CronInput (JsInputBase<string>) - cron schedule editor
    cron-utils.ts         # Pure cron parsing/validation/description utilities
  dialogs/               # Dialog UIs (formula-lines, add-new-column)
    formula-lines.ts      # FormulaLinesDialog: wires the parts below
    formula-lines/        # Host (source items, viewer rules), ItemList, Table, Preview (+ viewer-axes,
                          # preview-viewer), Editor (+ item / formula-line / annotation-region forms),
                          # CreationControl (+ item-factory, history), axis-sync, formulas, types
  search/                # Search providers (entity-search, power-search, help)
  workers/               # Web workers (exceljs-worker)
  tests/                 # Test files
css/
  power-pack.css         # All styles, loaded via "sources" in package.json
```
