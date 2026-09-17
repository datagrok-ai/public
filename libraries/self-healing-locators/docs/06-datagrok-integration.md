# 06 — Datagrok Platform Integration

> Status: compact. Captures specifically what we use from the platform
> and what we recommend changing in the platform.

## What we use from the platform today

**`Widget.getWidgetStatus()`.** Returns the widget's runtime structure
intended for automated testing and introspection (per the JS API). We
read `status.widgetType` and, where applicable, `status.viewerType` into
the fingerprint. This survives DOM-level reorganization because it
reflects the platform's own model of the page.

**The accessibility tree.** Datagrok's UI primitives (`ui.button`,
`ui.input.*`, `ui.dialog`, ribbon items, etc.) emit standard ARIA roles
and accessible names. We read these via the test driver's a11y APIs
(Puppeteer `page.accessibility.snapshot()`, Playwright's `getByRole`).

**Inspector data source.** The Inspector tool (`Alt + I`) exposes the
registered widget tree. Our capture utility taps the same data source
where possible so fingerprints reflect the platform's model rather than
just the DOM.

## What we recommend (separate, non-blocking work)

These are not required for the system to work. They make it work better.

1. **`data-testid` convention on `ui.*` primitives.** When the caller
   passes a logical name (`ui.button('Run analysis', ...)`), the
   primitive should set `data-testid="ui-button-run-analysis"` (slugified)
   by default. Skippable via opt-out. This gives us strong Tier 1
   anchors without authors thinking about it.
2. **Stable role/name on dynamic UI.** Dialog headers, viewer toolbars,
   ribbon items — anywhere the visible text doubles as the accessible
   name, ensure both are set explicitly. Our resolver will benefit;
   accessibility users will also benefit.
3. **`getWidgetStatus()` for non-viewer widgets.** Currently most
   reliable on viewers. Extending to dialogs, ribbon panels, and toolbox
   items would give us additional Datagrok-specific anchors.
4. **A `Test ID` mode in the Inspector.** A toggle that overlays
   suggested anchor names on the page would help test authors choose
   good explicit anchors and would let our capture tool offer them
   automatically.

These recommendations should be tracked separately. Any of them landing
makes the deterministic resolver better; none of them is required for
the v1 of this library.

## Constraints we accept from the platform

- **The library is a downstream consumer.** We do not assume the platform
  will change for us. Recommendations above are nice-to-have.
- **Tests run against a live platform instance.** Snapshots are what we
  capture, not what we control. The platform owns the source of truth
  about UI semantics; we just observe it.

## API touchpoints

| Use case | Platform API | Notes |
|---|---|---|
| Read widget type from element | `DG.Widget.fromRoot(el).getWidgetStatus()` | |
| Read accessibility tree | driver-native (Puppeteer/Playwright) | |
| Identify view boundary | `grok.shell.v` (current view) | Used as parent-chain stop |
| Read viewer type | `viewer.type` | When element is inside a viewer |
| Detect dialog containment | climb to nearest `[role="dialog"]` | |

If any of these become unstable across platform versions, the capture
utility records the platform version in the fingerprint
(`runMetadata.platformVersion` in the queue, separately in the
fingerprint via `hints.platformVersion`) so the offline tool can detect
incompatibilities and skip rather than guess.
