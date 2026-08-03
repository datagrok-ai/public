# UsageAnalysis

Datagrok plugin providing the **Usage Analysis** application — dashboards over platform usage,
events, sessions, errors, and operational health.

## Feature docs

- Tabs architecture (app entry, `ViewHandler`/`MultiView`, shared toolbox, viewers, drilldown, routing):
  [src/tabs/TABS.md](src/tabs/TABS.md)
- Metrics tab (operational-health dashboard; storage snapshot job and disk stats):
  [src/tabs/METRICS.md](src/tabs/METRICS.md)
- Release-readiness dashboard (dev-only `/release` app; Tests/Stress/Vulnerabilities/Tickets + overview):
  [src/release/RELEASE.md](src/release/RELEASE.md)
