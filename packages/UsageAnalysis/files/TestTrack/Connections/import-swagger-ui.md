---
feature: connections
realizes: []
priority: p2
target_layer: manual-only
coverage_type: edge
manual_only_reason: |
  OS-originated drag-and-drop of the swagger file into a browser tab must be
  done by hand.
related_bugs: []
---

# Import Swagger — manual UI checks

This is the **manual companion** to `import-swagger.md`. Covers the parts the
autotest cannot exercise: OS-originated drag-and-drop into a browser tab must
be done by hand.


## Pre-conditions

- The Samples package's `openweathermap.yaml` is downloaded locally — easiest path is to clone the public repo or "Download raw file" from the GitHub mirror

## Steps

1. Open the Datagrok browser tab
2. Drag `openweathermap.yaml` from your local file explorer onto the tab
3. Wait for the file-import dialog (or auto-toast) to recognise it as an OpenAPI swagger
4. Confirm the import — the platform should register an `OpenWeatherMap` connection under **Browse > Platform > Functions > OpenAPI**

## What to look for

- Toast / balloon that the file was recognised as a swagger (no red error balloons)
- The new `OpenWeatherMap` node appears under **OpenAPI** without page reload
- The node carries the standard connection icon (no broken / placeholder icons)
- Right-clicking the new node shows the standard connection context menu (Edit…, Delete…, etc.)

## Cleanup

- Right-click the imported `OpenWeatherMap` node under **Browse > Platform > Functions > OpenAPI**
  → **Delete...** → confirm, so the next import starts from a clean tree.

---
{
  "order": 5,
  "datasets": []
}
