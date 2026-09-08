---
feature: queries
realizes: [views.queries, viewers.scatter-plot, viewers.correlation-plot]
priority: p2
target_layer: manual-only
coverage_type: edge
manual_only_reason: |
  Layout-tab viewers are added by drag-and-drop, only real typing exercises
  the Post-Process editor, and the run's green info balloon is verified in a
  real browser session.
related_bugs: []
---

# Query Post-Processing — manual UI checks

This is the **manual companion** to `query-postprocessing.md`. The automated
companion covers the parts that map to stable automation; this file covers the
rest, all of which need a human in front of the browser.

The manual steps below cover three carve-outs the autotest could not exercise:

1. **Layout-tab viewers** — adding a Scatterplot and a Correlation Plot to the
   Layout-tab preview requires drag-drop, which must be done by hand.
2. **Typing into the Post-Process tab editor** — only real typing in the
   editor exercises this path.
3. **Run + green info balloon `77`** — the end-to-end runtime assertion lives
   here: in a real browser session, Right-click → Run fires the green balloon
   correctly.

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

## Cleanup

* Revert the `Test_Postprocessing` query to its pre-run state: right-click →
  **Edit...**; on the **Layout** tab remove the Scatterplot and Correlation Plot
  added in steps 7–8, and on the **Post-Process** tab remove the
  `grok.shell.info(result.rowCount);` line added in step 3; then **Save**.
  Alternatively, delete the whole query (right-click → **Delete...**) — the
  autotest or steps 1–6 of `query-postprocessing.md` re-create it for the next run.
* Close All.

---
{
  "order": 13,
  "datasets": []
}
