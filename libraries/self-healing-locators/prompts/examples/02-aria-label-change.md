# Example 02 — Accessible Name Change (MID band)

The platform changed the accessible name on the dialog OK button from
"OK" to "Confirm". No other changes. The system should heal but flag.

## Input

### Fingerprint
```json
{
  "schemaVersion": "1",
  "anchorName": "delete-confirm-ok",
  "semantic": {
    "ariaLabel": "OK",
    "ariaRole": "button",
    "name": "OK",
    "dgWidgetType": "Dialog"
  },
  "structural": {
    "tagName": "button",
    "stableClasses": ["dlg-button", "dlg-button-default"],
    "parentChain": [
      { "tagName": "div", "stableClasses": ["dlg-buttons"] },
      { "tagName": "div", "ariaRole": "dialog", "dgWidgetType": "Dialog" }
    ]
  },
  "text": { "visibleText": "OK" },
  "hints": {
    "textIsLikelyDynamic": false
  }
}
```

### Failed selector
`[aria-label="OK"]`

### Already tried
- `[aria-label="OK"]`
- `text="OK"`

### Accessibility tree (excerpt)
```
[role="dialog" dg-widget-type="Dialog" name="Delete project?"]
  [role="heading" level=1 name="Delete project?"]
  [role="button" name="Confirm"]
  [role="button" name="Cancel"]
```

## Expected response

```json
{
  "candidates": [
    {
      "selector": "role=dialog >> role=button[name=\"Confirm\"]",
      "selector_type": "role",
      "matched_attributes": ["ariaRole", "dgWidgetType", "parentChain"],
      "rationale": "Same dialog (dgWidgetType + structural parent), button role preserved, accessible name changed from 'OK' to 'Confirm'. The default action button position in the dialog matches.",
      "risk_notes": "Name changed; flag for review. Verify this is the intended affirmative button and not a different action."
    },
    {
      "selector": ".dlg-buttons > .dlg-button-default",
      "selector_type": "css",
      "matched_attributes": ["stableClasses", "parentChain"],
      "rationale": "Backup using the default-button class which has historically tracked the affirmative action across renamings.",
      "risk_notes": "Class-only fallback — weaker than role-based."
    }
  ],
  "unable_to_match": false
}
```

## What the validator does next

1. The `role=dialog >> role=button[name="Confirm"]` candidate matches
   one element. Runtime scorer:
   - `ariaRole` match: 0.85
   - `dgWidgetType` match: 0.80
   - `parentChain` match: 1.00 (depth 2)
   - `stableClasses` match: 0.66 (`dlg-button` matches, `dlg-button-default`
     matches → 2/2 = 1.0; minus mild penalty if stableClasses semantics
     in the tool consider both equally weighted; let's say full)
   - `ariaLabel` mismatch (OK vs Confirm): 0
   - `name` mismatch: 0
   - `visibleText` mismatch: 0

   Composite score lands around 0.62 — MID band.
2. Decision: heal AND flag.
3. The healing PR row reads `score=0.62, band=MID, flag=review`. The
   reviewer sees the rationale and confirms the rename was intentional.

This is the exact case where the policy doc § 4.1 ("synonym wording
changes") applies: text was decisive, so we flag.
