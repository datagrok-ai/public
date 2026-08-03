---
feature: connections
target_layer: playwright
coverage_type: edge
priority: p2
realizes: []
realized_as: []
related_bugs: []
---

# Import Swagger — manual UI checks

This is the **manual companion** to `import-swagger.md`. Covers the parts that
`import-swagger.test.ts` cannot exercise because OS-originated drag-and-drop into
a browser tab cannot be synthesised by Playwright (the browser only accepts a
real `dataTransfer` from the OS shell).

The autotest covers:

- The OpenWeatherMap connection appearing under **Browse > Platform > Functions > OpenAPI** after the YAML is ingested
- Right-click → **Edit...** dialog opening, ApiKey input accepting the real keyboard, **SAVE**
- Right-click on a query → **Run** producing either a result grid or a non-error balloon

The manual steps below cover the part the autotest skips.

## Pre-conditions

- Logged into Datagrok (dev or release)
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
