# API Samples changelog

## v.next

* GROK-20601: Replaced dapi/domains-dapi2.js (the generated dapi2 domains client is removed at typed-surface parity) with dapi/domains-schema.js — schema lifecycle handle: manifest export + whole-schema audit
* GROK-20600: dapi/domains-batch.js cleanup now uses one `deleteWhere` filtered bulk delete instead of the query + per-row delete loop
* GROK-20591: Added dapi/domains-filters.js — batched facets (filter-panel counts) and shareable saved filter presets on a domain table
* GROK-20572: Added data-access/db/query_builder_domain.js — TableQueryBuilder over a domain schema's virtual connection
* GROK-20316: Added dapi/domains-dapi2.js — querying domain-table rows via the generated grok.dapi2 REST client
* GROK-20315: Added dapi/domains-batch.js, domains-transaction.js, domains-aggregate.js, domains-dataframe.js, domains-idempotency.js — batch upsert, multi-entity transactions, aggregation, queryDf → grid, idempotent retries + optimistic concurrency
* GROK-20307: Added dapi/domains.js — domain-table row CRUD via grok.dapi.domains

## 1.2.1 (2-24-11-24)

* New samples, version bump

## 1.0.12 (2024-10-17)

* Add samples before for scripts
* Sample fixes 

## 1.0.11 (2024-07-09)

* Add samples for HelmInput

## 1.0.7 (2023-07-27)

* Dependency: datagarok-api >= 1.16.0*

### Features

* All examples now have auto tests 
