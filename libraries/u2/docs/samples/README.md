# u2 samples — "Build with u2"

Copy-ready u2 code for the everyday jobs: inputs, forms, dialogs, panels, lists, tables, menus, async data,
reactive state, specs and the platform shell. Each sample is the whole code, written as function components
(`component((self) => …)`), views through `ViewHost`, imports from the entry points only.

| File | What it is |
|---|---|
| `<id>.ts` | the sample: the pattern, then — under a `// Live demo` line — its demo data and `render(host)` for a live one |
| `<id>.md` | its write-up: front matter (`title`, `group`, `order`, `summary`, `live: yes \| no`, `height`, `old`), the description, then `## What it takes` (`{lines}`, `{cssLines}` and `{oldLines}` are filled in from the code where the sample is shown) |
| `<id>.css` | the sample's own CSS (a height, a colour), imported by it as `./<id>.css`; never a rule over u2's classes |
| `<id>.old.ts` | the same job with the platform's `ui.*` today, for the comparison; its tab title is `old:` in the `.md` |
| `tsconfig.json` | compiles the folder against u2's sources (`tsc -p docs/samples` from `libraries/u2`) |

`tests/recipes.test.js` compiles every file here and runs every live sample headless — strict usage, nothing left
live, the operability probe clean. A sample that cannot run without the platform says why in that test's
`HEADLESS` table.
