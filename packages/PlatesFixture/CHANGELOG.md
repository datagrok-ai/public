# PlatesFixture changelog

## v.next

* GROK-20298: Added `plate.source_query` (`ref: Core.queries`) — a soft cross-schema ref into a Core table (index, no FK; per-row `fk` write probe; filter/expand/facet hop into the query; a deleted query leaves the id dangling)
* GROK-20590: Added the `filters` manifest section to `plate_well` (default filter panel: categorical, histogram with bins, FK-path column with a label, system-column range)
* GROK-20308: Introduced the plates domain-schema fixture — all three security modes, constraints, referential actions, and audit flags exercised by CI tests over `grok.dapi.domains`
* GROK-20314: Added the depth-2 `well_reading` table (master mode delegating to `plate_well`) with a CI test covering transitive security and cascade through the chain

## 0.0.1 (2026-07-04)
