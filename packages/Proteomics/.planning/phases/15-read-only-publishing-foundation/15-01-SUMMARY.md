---
phase: 15
plan: 01
status: complete
type: implementation
---

# Plan 15-01 Summary — Publish-State Module

Built `src/publishing/publish-state.ts`, the single source of truth for Phase 15's published-tag namespace. Every downstream plan (02 trim, 03 assert, 04 orchestrator, 05 dialog, 06 panel, 07 wireup, 08 test) imports its constants and helpers from this file — no other file in the package may inline `proteomics.published*` tag string literals.

## Exports

**Constants**
- `DE_COMPLETE_TAG` — `'proteomics.de_complete'` (I-10 fix)
- `PUBLISHED_TAGS` — 13-entry const record, all keys prefixed `proteomics.`
- `META_COLUMNS` — 13-entry const record, parallel keyspace to `PUBLISHED_TAGS`, names prefixed `_meta_published_*`

**Interfaces**
- `PublishedMetadata` — 12 fields (target, publishedAt, publishedBy/Email, deMethod, fcThreshold, pThreshold, version, publishId, includesEnrichment, supersedes, supersededBy)
- `PublishOptions` — 5 fields (target, reviewerGroup, note, priorVersion, optional umbrellaName for test injection)
- `MailtoOptions` — 4 fields (sharerEmail, sharerName, projectName, publishedDateStr)

**Functions**
- `isPublished(df)` — boolean gate, one-line mirror of `proteomics.de_complete === 'true'` idiom
- `getPublishedMetadata(df)` — column-FIRST / tag-SECOND read order (Pitfall 3 belt-and-braces); returns `null` on non-published DataFrames; per-field corruption falls back gracefully
- `setPublishedTags(df, meta)` — writes every required `PUBLISHED_TAGS.*` entry as a string tag; intentionally omits null `supersedes`/`supersededBy` so reader can disambiguate "missing" from "empty"
- `slugifyTarget(raw)` — D-01 sanitization: charset `[A-Za-z0-9._-]`, case PRESERVED, cap 64, empty → `'unnamed'`
- `findPriorShare(target, group)` — async, returns the highest-version prior `DG.Project` matching `Proteomics-Review-<slug>-v*` or `null`; uses `like` smart-filter (PRIMARY) with `list()` + client-side filter as FALLBACK
- `buildMailtoUrl(opts)` — pure `encodeURIComponent`-safe `mailto:` builder; co-located here rather than in `share-dialog.ts` because Plan 06 (wave 2) needs it before Plan 05 (wave 4) lands (B-1 wave-cycle fix)

## Filter syntax used

`findPriorShare` uses `grok.dapi.projects.filter('name like "<pattern>"').list()` as PRIMARY. Spike 15-00 confirmed `like` works against `release/1.27.3` (count 1, firstName matched). The `list()` + client-side `startsWith` fallback is retained for defense in depth in case the platform smart-filter parser changes.

## No string-literal drift (production code)

`grep -rE "proteomics\.published" src --include='*.ts' | grep -v 'publishing/publish-state.ts' | grep -v 'tests/publish-spike.ts'` returns zero matches at this commit. The only inlines are in `src/tests/publish-spike.ts` (the Wave 0 enumeration probe, category `Publishing-Spike`, excluded from CI) — that file intentionally hardcoded the namespace to verify it independently of the helper module. Every subsequent production plan must import from `'./publish-state'` rather than inlining the strings.

## buildMailtoUrl co-location

Per the B-1 fix in 15-01-PLAN.md, `buildMailtoUrl` lives here and not in `share-dialog.ts` because:
1. Plan 06 (reviewer panel, wave 2) needs it for PUB-13
2. share-dialog.ts is wave 4 — importing wave 4 from wave 2 would break the wave graph
3. The function is pure (no DOM, no `grok.shell` calls, no dialog dependency) so it has zero coupling to the dialog code

Plan 05 will import `buildMailtoUrl` from `./publish-state`, not the other way around.

## Verification

- Project-wide `tsc --noEmit` passes with no errors mentioning `publish-state.ts`
- File exports: 6 functions, 3 interfaces, 2 const records, 1 string constant (per the artifacts spec)
- All 13 entries in `PUBLISHED_TAGS` start with `proteomics.`
- `PUBLISHED_TAGS` and `META_COLUMNS` have identical keyspaces (13 entries each, same key names)

## Threat-model status

T-15-03 (slug injection): mitigated by `slugifyTarget` regex + 64-char cap; raw target is preserved in tag/column only and never used for paths/URLs.
T-15-05 (findPriorShare scope): mitigated by platform ACL implicit scoping on `dapi.projects.filter`; fallback `list()` is also platform-scoped.
T-15-SP-06 (mailto header injection): mitigated by wrapping every variable part in `encodeURIComponent`.

## Output

`src/publishing/publish-state.ts` — 232 lines, type-checks clean under strict mode.
