# Example 01 — `data-testid` Rename

A common, easy case. The platform team renamed `data-testid="submit"`
to `data-testid="submit-btn"` for consistency. Everything else is
unchanged.

## Input

### Fingerprint
```json
{
  "schemaVersion": "1",
  "anchorName": "login-submit",
  "semantic": {
    "testId": "submit",
    "ariaLabel": "Submit",
    "ariaRole": "button",
    "name": "Submit"
  },
  "structural": {
    "tagName": "button",
    "stableClasses": ["primary-action"],
    "parentChain": [
      { "tagName": "div", "stableClasses": ["form-actions"] },
      { "tagName": "form", "ariaRole": "form" }
    ]
  },
  "text": { "visibleText": "Submit" },
  "hints": {
    "domIdLooksStable": false,
    "classesAreCssInJs": false,
    "textIsLikelyDynamic": false
  }
}
```

### Failed selector
`[data-testid="submit"]`

### Already tried
- `[data-testid="submit"]`
- `#submit`
- `[role="button"][aria-label="Submit"]`

(In practice, Tier 1 of the resolver would have found this case via
`getByRole('button', { name: 'Submit' })`. The example assumes that
hit was somehow weak — perhaps another `Submit` exists on the page.)

### Accessibility tree (excerpt)
```
[role="form"]
  [role="group" name="form-actions"]
    [role="button" name="Submit" testid="submit-btn"]
    [role="button" name="Cancel" testid="cancel-btn"]
[role="navigation"]
  [role="button" name="Submit feedback"]
```

## Expected response

```json
{
  "candidates": [
    {
      "selector": "[data-testid=\"submit-btn\"]",
      "selector_type": "testid",
      "matched_attributes": ["testId"],
      "rationale": "Same role and name; testid renamed from 'submit' to 'submit-btn'. Unique under form-actions group.",
      "risk_notes": "Note that 'Submit feedback' button exists in nav; this selector avoids it via the testid."
    },
    {
      "selector": "form [role=\"button\"][aria-label=\"Submit\"]",
      "selector_type": "role",
      "matched_attributes": ["ariaLabel", "ariaRole", "parentChain"],
      "rationale": "Backup using the form scope and accessible name. Does not depend on testid stability."
    }
  ],
  "unable_to_match": false
}
```

## What the validator does next

1. `[data-testid="submit-btn"]` → exactly one match → score ~0.92
   (testId 1.0, ariaLabel 1.0, ariaRole+name 1.0, stableClasses 0.66,
   parentChain 1.0). HIGH band. Heal.
2. The second candidate is also strong but not picked because the first
   one wins.

Outcome: PR opens with `[data-testid="submit"]` →
`[data-testid="submit-btn"]`.
