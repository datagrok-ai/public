# Charts behavioral tests

Gherkin features under `features/`, compiled by `@datagrok-libraries/bdd` (`public/libraries/bdd`)
into the Playwright specs under `generated/` — committed, never edited by hand.

| Feature | What it claims |
|---|---|
| `viewers/viewer-rules-edge-cases` | the rules every viewer shares (no rows, a column removed or renamed, a layout naming a missing column, a click filter left by a closed viewer or a removed column, a value out of its declared bounds, an empty category, the 500-category limit), one Examples row per Charts viewer; the core viewers' rows are in `packages/UsageAnalysis/bdd/features/viewers/viewer-rules-edge-cases.feature` |
| `viewers/radar/radar-edge-cases` | Values columns removed from the table, every row filtered out, every row deleted |
| `viewers/sankey-chord/sankey-chord-edge-cases` | Source equal to Target, more than 1500 Chord pairs, the Source / From column removed |
| `viewers/sunburst/sunburst-edge-cases` | an empty hierarchy, a hierarchy column removed, On Click switched from Filter to Select |
| `viewers/timelines/timelines-edge-cases` | a cleared Color, 200-character lane names, no rows, an event whose End is before its Start |
| `viewers/tree/tree-edge-cases` | a column renamed under a filter, a hierarchy column removed, empty values, On Click switched from Filter to Select |
| `viewers/word-cloud/word-cloud-edge-cases` | the word column removed |

What the features cannot claim yet is in the `MISSING.md` next to each feature.

From a fresh checkout of `public`, against a local stand on `http://localhost:8888` with this
package published:

```bash
cd public/libraries/bdd && npm ci && npm run build   # the library (a path dependency; dist/ is not committed)
npx playwright install chromium                      # its browser, once per machine
cd ../../packages/Charts && npm install              # the package
npx grok-bdd link                                    # ONE Playwright: the library's copy into node_modules (redo after every npm install)
npx webpack && grok publish localhost                # the features read the viewers' own status
npx grok-bdd run --reporter=list                     # compile --check, then Playwright (one page per worker; --workers N)
```
