# BenchlingLink

BenchlingLink is a [Datagrok](https://datagrok.ai) package that integrates with
[Benchling](https://www.benchling.com/), letting you browse tenant content —
AA sequences, DNA sequences, DNA oligos, molecules, mixtures, plates, assay results,
assay runs, and projects — and post new entities back to the registry, all from inside
Datagrok.

## Status

> **This package currently runs against mock fixtures, not a live Benchling tenant.**
> Every list/create call in `src/*Api.ts` returns a hard-coded 100-row fixture; the real
> `grok.dapi.fetchProxy(...)` calls are commented out alongside them. The UI, input
> forms, and DataFrame conversion are real — only the HTTP hop is stubbed.
>
> Before switching to live HTTP, see **[Activating real HTTP](#activating-real-http)**
> and the auth caveat below.

## Features

- Browse Benchling registry content in a sidebar tree under **Browse > Apps > Chem > Benchling**
- Typed input form per resource with all of Benchling's documented filter params
- List responses flattened into Datagrok DataFrames (nested `creator.name` becomes a `"creator.name"` column; arrays of objects become comma-joined identifier CSVs)
- `createX` counterparts for posting new entities (AA sequence, DNA sequence, DNA oligo, molecule, mixture, plate, assay result, assay run)
- Auto-registered Datagrok functions (`BenchlingLink:getAASequences`, `BenchlingLink:createMolecule`, …) callable from pipelines, scripts, and the platform search bar

## Resources wrapped

| Resource       | List                              | Create                               |
|----------------|-----------------------------------|--------------------------------------|
| AA Sequences   | `BenchlingLink:getAASequences`    | `BenchlingLink:createAASequence`     |
| DNA Sequences  | `BenchlingLink:getDNASequences`   | `BenchlingLink:createDNASequence`    |
| DNA Oligos     | `BenchlingLink:getDnaOligos`      | `BenchlingLink:createDnaOligo`       |
| Molecules      | `BenchlingLink:getMolecules`      | `BenchlingLink:createMolecule`       |
| Mixtures       | `BenchlingLink:getMixtures`       | `BenchlingLink:createMixture`        |
| Plates         | `BenchlingLink:getPlates`         | `BenchlingLink:createPlate`          |
| Assay Results  | `BenchlingLink:getAssayResults`   | `BenchlingLink:createAssayResult`    |
| Assay Runs     | `BenchlingLink:getAssayRuns`      | `BenchlingLink:createAssayRun`       |
| Projects       | `BenchlingLink:getProjects`       | —                                    |

The upstream Benchling REST API v2 OpenAPI spec lives at `swaggers/benchling_openapi.yaml`; it is the source of truth for field names, required params, and response shapes.

## Usage

Open the app from **Browse > Apps > Chem > Benchling**. The sidebar tree shows one
node per resource. Selecting a node opens the standard platform preview for the matching
`getX` function — an input form above a live grid. Fill in any filters and run the search.

Create operations are exposed as separate `createX` functions. They accept array/object
fields (`aliases`, `annotations`, `fields`, `customFields`, `ingredients`, `wells`, `fieldValidation`, …)
as JSON strings and parse them inside the wrapper.

## Activating real HTTP

1. **Verify authentication first** — see the caveat below.
2. In the target `src/<resource>Api.ts`, uncomment the `// const token = … / url / response / data / df`
   block inside `queryXxx` and `postXxx`.
3. Replace the `'YOUR_BENCHLING_API_TOKEN'` placeholder with `await getApiKey()` and add
   `import { getApiKey } from './credentials';`.
4. Delete (or gate behind a mock flag) the `return mockXxx...` line.
5. Configure the `apiKey` credential on the deployed package (via `grok s` or the package settings panel).

The `buildQueryString` helper in `src/utils.ts` automatically translates TS underscore
suffixes (`_anyOf`, `_gt`, `_lte`, `_mol`, `_smiles`, …) into Benchling's dotted
wire names (`.anyOf`, `.gt`, `.lte`, …) — no per-resource rewrite required.

## ⚠️ Auth mechanism is UNVERIFIED

The current `src/credentials.ts` assumes a single static API key stored as the `apiKey`
package credential and injected as `Authorization: Bearer <apiKey>`. **This has never
been tested against a live Benchling tenant.**

The swagger's `POST /token` endpoint (line ~17482 of `swaggers/benchling_openapi.yaml`)
describes OAuth2 client-credentials auth (`basicClientIdSecretAuth`): a `clientId` +
`clientSecret` pair exchanged for a short-lived access token that must be refreshed on
expiry. Different tenants may enforce different schemes.

Before wiring real HTTP:

1. Confirm with the target Benchling tenant admin which scheme is enabled — static API
   token, OAuth2 client-credentials, or both.
2. If OAuth2 is required, `credentials.ts` must be rewritten end-to-end (credential
   parameters become `clientId` + `clientSecret`; add a `POST /token` call;
   cache the access token with its `expires_in` TTL; refresh on 401).
3. Even in the static-bearer case, revisit the module-level `apiKey` cache — it holds
   the credential for the full page session, so server-side rotations take effect only
   after reload; concurrent first callers race on `_package.getCredentials()`.

Do not treat the current `getApiKey()` as load-bearing until the auth layer has been
verified against a live tenant.

## Package structure

```text
src/
  package.ts           App + tree browser + every @grok.decorators.func wrapper
  <resource>Api.ts     Per-resource interfaces + query/post + mock fixture
  utils.ts             dataFrameFromObjects + buildQueryString + translateParamName + randomDnaSequence
  credentials.ts       getApiKey() (currently unused, awaiting real HTTP activation)
  types.ts             Shared types (UserSummary, Organization, ArchiveRecord)
files/
  benchling_mock_molecules.csv   SMILES source for the Molecules mock fixture
swaggers/
  benchling_openapi.yaml         Upstream Benchling OpenAPI spec
```

For the developer-facing conventions (tree order, per-resource template, schema-matched
create shapes, Benchling glossary), see [CLAUDE.md](./CLAUDE.md).

## See also

- [Benchling Developer Platform documentation](https://docs.benchling.com/docs)
- [Benchling API reference](https://benchling.com/api/reference)
- [Datagrok documentation](https://datagrok.ai/help/)
