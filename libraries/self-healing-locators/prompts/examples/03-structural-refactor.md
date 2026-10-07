# Example 03 — Structural Refactor

The viewer toolbar was wrapped in a new container, classes were
re-hashed, and `data-testid` was lost in the migration. The accessible
name and role remain.

## Input

### Fingerprint
```json
{
  "schemaVersion": "1",
  "anchorName": "scatter-toolbar-settings",
  "semantic": {
    "testId": "scatter-toolbar-settings",
    "ariaLabel": "Settings",
    "ariaRole": "button",
    "name": "Settings",
    "dgWidgetType": "Viewer",
    "dgViewerType": "scatter-plot"
  },
  "structural": {
    "tagName": "button",
    "stableClasses": ["toolbar-button", "icon-cog"],
    "parentChain": [
      { "tagName": "div", "stableClasses": ["viewer-toolbar"] },
      { "tagName": "div", "dgWidgetType": "Viewer" }
    ]
  },
  "text": { "title": "Settings" },
  "visual": {
    "boundingBox": { "x": 1180, "y": 60, "width": 28, "height": 28 },
    "viewportSize": { "w": 1440, "h": 900 }
  },
  "hints": {
    "textIsLikelyDynamic": false
  }
}
```

### Failed selector
`[data-testid="scatter-toolbar-settings"]`

### Already tried
- `[data-testid="scatter-toolbar-settings"]`
- `.toolbar-button.icon-cog`
- `[role="button"][aria-label="Settings"]`

(Tier 3 found three "Settings" buttons on the page — the scatter
toolbar, the bar chart toolbar, and the page-level settings menu — so
ambiguity demoted it to LOW. Hence the case is here.)

### Accessibility tree (excerpt)
```
[dg-widget-type="Viewer" dg-viewer-type="scatter-plot" name="Scatter plot"]
  [role="toolbar" name="Scatter plot toolbar"]
    [class="vtb-controls"]
      [role="button" name="Settings"]
      [role="button" name="Properties"]
      [role="button" name="Close"]

[dg-widget-type="Viewer" dg-viewer-type="bar-chart"]
  [role="toolbar"]
    [role="button" name="Settings"]
    ...

[role="navigation"]
  [role="button" name="Settings"]
```

## Expected response

```json
{
  "candidates": [
    {
      "selector": "[dg-viewer-type=\"scatter-plot\"] [role=\"toolbar\"] [role=\"button\"][aria-label=\"Settings\"]",
      "selector_type": "css",
      "matched_attributes": [
        "dgWidgetType",
        "dgViewerType",
        "ariaRole",
        "ariaLabel",
        "parentChain"
      ],
      "rationale": "Scopes by dgViewerType to disambiguate the three 'Settings' buttons on the page. Role and accessible name are preserved.",
      "risk_notes": "If the toolbar role is dropped in a future change, this would break; consider adding a stable testid in the platform."
    },
    {
      "selector": "[dg-viewer-type=\"scatter-plot\"] .vtb-controls [role=\"button\"][aria-label=\"Settings\"]",
      "selector_type": "css",
      "matched_attributes": [
        "dgViewerType",
        "stableClasses",
        "ariaRole",
        "ariaLabel",
        "parentChain"
      ],
      "rationale": "Adds the inner controls class as additional structural scoping. Slightly more specific than the toolbar-role version."
    }
  ],
  "unable_to_match": false
}
```

## What the validator does next

1. First candidate: matches one element. Runtime scorer:
   - `dgWidgetType` + `dgViewerType` match: 0.80
   - `ariaRole` match: combined into role+name 0.85
   - `name` match: contributes via role+name
   - `ariaLabel` match: 0.85
   - `parentChain` match: 1.0 (rooted on viewer)
   - testId absent on both sides: doesn't contribute either way
   - visual: bounding box close enough → 0.25 if visual layer enabled
   - Composite ~0.86 → HIGH band.
2. Decision: heal silently.
3. Note the second candidate's rationale mentions `vtb-controls` — the
   validator would check whether that class survived. If it did and
   the score is similar, the first candidate still wins on simplicity.

This is a typical "structural refactor" case. Without `dgViewerType` as
a Datagrok-specific anchor, the resolver would have failed at LOW. The
platform's semantic layer is what makes this case healable.
