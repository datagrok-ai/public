# Domain UI changelog

## v.next

* GROK-20298: Initial release — `DomainFrameEditor` (batch editing of domain rows with the
  `~state` / `~changes` / `~errors` service-column state, one-transaction save, standard
  version-conflict dialog, rebuild-on-refresh) and `DomainGrid` (platform-decorated editable
  grid with change highlighting, per-cell validation markers, deleted-row exclusion, and
  capability-driven read-only degradation); re-exports of the core domain surface
