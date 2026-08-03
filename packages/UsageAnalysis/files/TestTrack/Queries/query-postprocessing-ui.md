---
feature: queries
target_layer: playwright
coverage_type: edge
priority: p2
realizes: []
realized_as: []
related_bugs: []
---

# Query Post-Processing — manual UI checks

This is the **manual companion** to `query-postprocessing.md`. The autotest
`query-postprocessing.test.ts` covers the parts that map to stable DOM automation;
this file covers the rest, all of which need a human in front of the browser.

The autotest covers: create the query via UI, set Name + SQL, run via Play, Save,
patch `postProcessScript` via dapi (UI fallback — see SCOPE NOTE 2 in the spec
file), and verify on round-trip that the Post-Process tab mounts correctly with the
saved body.

The manual steps below cover three carve-outs the autotest could not exercise:

1. **Layout-tab viewers** — adding a Scatterplot and a Correlation Plot to the
   Layout-tab preview requires drag-drop or a JS handle on the nested TableView,
   neither of which respond to synthesised browser events.
2. **Typing into the Post-Process tab editor** — the CodeMirror inside ScriptView
   does not propagate JS-injected values into `DataQueryView.postProcessScript`,
   so the autotest patches via dapi instead. A real human typing in the editor
   does work, and that path is what we exercise here.
3. **Run + green info balloon `77`** — Right-click → Run on a query whose
   `postProcessScript` was patched in the same session raises a red
   "Handler for null is not registered" error in our automation
   (the client function registry doesn't rebind the new handler). A real fresh
   browser session does NOT see this — the balloon fires correctly. So the
   end-to-end runtime assertion lives here.

## Pre-conditions

* You have a `Test_Postprocessing` query saved (run the autotest once, or follow
  steps 1–6 of `query-postprocessing.md` manually).
* The query's `postProcessScript` already contains
  `grok.shell.info(result.rowCount);` — the autotest sets this for you, or you can
  type it manually inside the Post-Process tab.

## Steps

### Editing the post-process via UI (covers carve-out 2)

1. Right-click the `Test_Postprocessing` query in **Browse > Databases > Postgres
   > NorthwindTest** and select **Edit...**
2. Switch to the **Post-Process** tab — the editor shows a JavaScript template
   with `result` (input dataframe) and `out` (output dataframe).
3. On line 7 (or anywhere inside the script body), add:
   ```javascript
   grok.shell.info(result.rowCount);
   ```
4. Click **Save** in the ribbon — the query updates silently (no dialog).

### Building a layout (covers carve-out 1)

5. Switch to the **Layout** tab.
6. Click **Run query** in the Layout-tab preview — the preview TableView populates.
7. Drag a **Scatterplot** viewer from Toolbox > Viewers onto the layout.
8. Drag a **Correlation Plot** viewer the same way.
9. **Save** the query.

### Running and verifying the post-process (covers carve-out 3)

10. **Close All**.
11. Go to **Browse > Databases > Postgres > NorthwindTest** and click the
    `Test_Postprocessing` query — verify in the preview that:
    * Both viewers (Scatterplot + Correlation Plot) render in the saved layout.
    * A green info balloon with `77` appears (the post-process ran on preview).
12. Right-click the query and select **Edit...**
13. Switch to the **Post-Process** tab and click **Run query** — verify the green
    info balloon with `77` appears.
14. Switch to the **Layout** tab and click **Run query** — verify the green info
    balloon with `77` appears AND both viewers are still displayed.

## What to look for

* The post-process you typed in step 3 persists across Close All + reopen.
* Both viewers persist into the saved layout and show up after Close All + reopen.
* The post-process info balloon (`77`) fires on every Run, regardless of which tab
  you press Run from.
* No red error balloons or console errors during any step.
