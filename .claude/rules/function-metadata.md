---
paths:
  - packages/**/src/**
  - packages/**/scripts/**
---

## Function Registration Pattern

Functions are declared via JSDoc-style metadata comments (NOT TypeScript decorators), processed by `grok api`:

```typescript
//name: myFunction
//description: What it does
//input: dataframe df {caption: Input}
//input: column col {type: numerical}
//input: string name {optional: true}
//output: dataframe result
//meta.role: panel
//meta.cache: client
export function myFunction(df: DG.DataFrame, col: DG.Column, name?: string): DG.DataFrame {
}
```

`grok api` generates `package.g.ts` and `package-api.ts` from these comments. `grok check` validates them.

## Function Roles

| Tag | Role |
|-----|------|
| `#app` | Application entry point |
| `#panel` | Info panel (extends context panel) |
| `#init` | Package initialization (runs once before first function call) |
| `#autostart` | Runs at platform startup |
| `#semTypeDetector` | Semantic type detector |
| `#cellRenderer` | Custom cell renderer (`cellRenderer-<type>`) |
| `#fileViewer` | File viewer for specific extensions (`fileViewer-<ext>`) |
| `#fileExporter` | Data export function |
| `#dashboard` | Dashboard |
| `#packageSettingsEditor` | Custom settings editor UI |
| `meta.role: converter` | Type converter (with `meta.inputRegexp` for matching) |

## Script Metadata (Python/R/Julia)

```python
#name: Descriptors
#language: python
#input: dataframe df
#input: column smiles {type: chemical/smiles}
#output: dataframe result
```

## SQL Query Metadata

```sql
--name: ProteinClassification
--connection: Chembl
--input: string target = "CHEMBL301"
SELECT * FROM protein_classification WHERE chembl_id = @target
--end
```
