---
feature: viewers
realizes_atlas: []   # untyped: nx cross-viewer integration chain
priority: p2
target_layer: manual-only
coverage_type: smoke
---

### Legend position
1. Open **SPGI**
2. Apply the old layout attached to [#3203](https://github.com/datagrok-ai/public/issues/3203)
   (download `SPGI_v2 - legend position old.zip` from the issue, unpack it, and drag
   the `.layout` file onto the view)
3. Open **Context Panel**, verify Legend Position = `right` for each viewer
4. Close All
---
{
  "order": 6,
  "datasets": ["System:DemoFiles/chem/SPGI.csv"]
}