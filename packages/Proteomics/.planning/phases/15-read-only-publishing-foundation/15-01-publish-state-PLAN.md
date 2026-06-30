---
phase: 15-read-only-publishing-foundation
plan: 01
type: execute
wave: 1
depends_on: ["15-00"]
files_modified:
  - src/publishing/publish-state.ts
autonomous: true
requirements: [PUB-04, PUB-07, PUB-10, PUB-11, PUB-13]
must_haves:
  truths:
    - "Any caller can ask `isPublished(df)` and get a deterministic boolean (true iff `proteomics.published === 'true'`)"
    - "The published tag namespace has exactly one source of truth — `publish-state.ts` — and no other file inlines `proteomics.published*` tag string literals"
    - "`findPriorShare(target, group)` returns the most recent published Project matching that target+group, or null if none exist"
    - "`getPublishedMetadata(df)` returns a typed object whose fields mirror the belt-and-braces metadata column names (reader uses column FIRST, tag SECOND per Pitfall 3)"
    - "The slug sanitization helper accepts any freeform target string and produces a safe `[A-Za-z0-9._-]` slug capped at 64 chars (D-01)"
    - "`buildMailtoUrl(opts)` is a pure helper exported here (PUB-13) so the reviewer panel can import it without depending on share-dialog.ts (which is wave 4 — would otherwise create a wave/dependency cycle for the wave-2 panel)"
    - "`DE_COMPLETE_TAG = 'proteomics.de_complete'` constant exported so callers (Plan 05, Plan 07) don't inline the string"
  artifacts:
    - path: "src/publishing/publish-state.ts"
      provides: "Single source of truth for publishing tags + helpers + buildMailtoUrl"
      exports: ["isPublished", "getPublishedMetadata", "setPublishedTags", "findPriorShare", "slugifyTarget", "buildMailtoUrl", "PUBLISHED_TAGS", "META_COLUMNS", "DE_COMPLETE_TAG", "PublishedMetadata", "PublishOptions", "MailtoOptions"]
  key_links:
    - from: "src/publishing/publish-state.ts"
      to: "grok.dapi.projects.filter"
      via: "findPriorShare smart-filter query"
      pattern: "grok\\.dapi\\.projects\\.filter"
    - from: "src/publishing/publish-state.ts"
      to: "src/utils/proteomics-types.ts"
      via: "SEMTYPE import"
      pattern: "import.*SEMTYPE.*proteomics-types"
---

<objective>
Build the single source of truth for Phase 15's published-tag namespace, plus pure helpers (`isPublished`, `getPublishedMetadata`, `setPublishedTags`, `findPriorShare`, `slugifyTarget`, `buildMailtoUrl`). Sibling pattern of `src/analysis/experiment-setup.ts`'s `getGroups`/`setGroups`: tag I/O lives in one module, every consumer imports from it, no other file inlines tag strings.

Purpose: prevent tag-string drift across `trim-dataframe.ts` (writer), `publish-project.ts` (orchestrator), `assert-published-shape.ts` (gate), `share-dialog.ts` (republish-detection consumer), and `published-analysis-panel.ts` (reader). A single typo in a tag string would silently break the panel; centralizing is the antidote.

**Wave-ordering note (revision):** `buildMailtoUrl` was originally placed in `share-dialog.ts` (Plan 05, wave 4), but Plan 06 (the reviewer panel, wave 2) needs it for PUB-13. Importing wave-4 from wave-2 would break the wave graph. Resolution: `buildMailtoUrl` is pure, has no UI / DOM / dialog dependency, and conceptually belongs alongside the other publishing primitives — so it lives here. Plan 05 imports it from this module (or re-exports as a convenience for any external caller).

Output: `src/publishing/publish-state.ts` exporting the 6 helpers + the typed shapes + the `PUBLISHED_TAGS` + `META_COLUMNS` constant records + the `DE_COMPLETE_TAG` constant. Pure module — no DOM, no UI calls.
</objective>

<execution_context>
@$HOME/.claude/get-shit-done/workflows/execute-plan.md
@$HOME/.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md
@.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md
@.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md
@packages/Proteomics/CLAUDE.md
@packages/Proteomics/src/analysis/experiment-setup.ts
@packages/Proteomics/src/utils/proteomics-types.ts
@packages/Proteomics/src/utils/column-detection.ts

<interfaces>
<!-- Pulled from src/analysis/experiment-setup.ts (analog) -->

From src/analysis/experiment-setup.ts:
  export interface GroupAssignment {
    group1: {name: string; columns: string[]};
    group2: {name: string; columns: string[]};
  }
  export function setGroups(df: DG.DataFrame, groups: GroupAssignment): void
  export function getGroups(df: DG.DataFrame): GroupAssignment | null

From datagrok-api:
  df.getTag(key: string): string | null
  df.setTag(key: string, value: string): void
  df.col(name: string): DG.Column | null    // single-row metadata read source
  col.get(rowIdx: number): unknown          // for reading metadata column[0]
  grok.dapi.projects.filter(smartFilter: string): HttpDataSource<DG.Project>
  grok.dapi.projects.list(): Promise<DG.Project[]>  // fallback if filter unsupported
  grok.dapi.groups: GroupsDataSource
  DG.Group: { id: string; friendlyName: string; ... }
</interfaces>
</context>

<tasks>

<task type="auto" tdd="false">
  <name>Task 1: Define typed shapes, PUBLISHED_TAGS, META_COLUMNS, DE_COMPLETE_TAG constants</name>
  <files>src/publishing/publish-state.ts</files>
  <read_first>
    - @packages/Proteomics/src/analysis/experiment-setup.ts (full file — analog for tag helper module)
    - @packages/Proteomics/CLAUDE.md (tag namespace conventions, function-naming prefixes)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (decisions D-01, D-04 — slug rules + supersede chain)
    - @.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md (Section 1 — analog and notes)
    - @packages/Proteomics/src/utils/proteomics-types.ts (SEMTYPE constants for typing)
  </read_first>
  <action>
Create `src/publishing/publish-state.ts`. Follow the exact import + comment + JSDoc convention of `src/analysis/experiment-setup.ts` lines 1-27. Mirror its tag-helper shape.

Declare the typed shapes:

`PublishedMetadata` interface — fields one-to-one with the belt-and-braces metadata columns Plan 02 will write. Required fields (per CONTEXT.md Claude's discretion + RESEARCH §"Pattern 3"):
  - `target: string` (raw user input, unsanitized — never use for paths/URLs)
  - `publishedAt: Date | string` (ISO string at write time; Date when read from `addNewDateTime` column)
  - `publishedBy: string` (`grok.shell.user.friendlyName`)
  - `publishedByEmail: string | null` (`grok.shell.user.email`, may be null per RESEARCH Open Question 5)
  - `deMethod: string` (one of `'limma' | 'deqms' | 't-test' | 'spectronaut'`)
  - `fcThreshold: number`
  - `pThreshold: number`
  - `version: number` (1, 2, 3, ...)
  - `publishId: string` (the published Project's id — set after `projects.save` returns)
  - `includesEnrichment: boolean`
  - `supersedes: string | null` (prior Project id, set by Plan 04 supersede step)
  - `supersededBy: string | null` (next Project id, set retroactively when republished)

`PublishOptions` interface — fields for `share-dialog.ts -> publishAnalysis(df, opts)` input:
  - `target: string` (raw freeform user input)
  - `reviewerGroup: DG.Group`
  - `note: string`
  - `priorVersion: DG.Project | null` (set by `findPriorShare` if a republish; null on fresh share)
  - `umbrellaName?: string` (OPTIONAL override for tests — defaults to `'Proteomics-Reviews'`; Plan 08 Test 5 uses this to inject a throwaway umbrella; production callers never set it)

`MailtoOptions` interface — fields for `buildMailtoUrl(opts)`:
  - `sharerEmail: string | null`
  - `sharerName: string`
  - `projectName: string`
  - `publishedDateStr: string`

Declare `PUBLISHED_TAGS` as a `const` readonly record. Every tag key uses the `proteomics.` prefix. Required entries (resolves PUB-11 + D-04 + Pitfall 3 belt-and-braces mirror, plus the volcano-threshold belt-and-braces from B-2 fix):
  - `PUBLISHED: 'proteomics.published'`
  - `PUBLISHED_AT: 'proteomics.published_at'`
  - `PUBLISHED_BY: 'proteomics.published_by'`
  - `PUBLISHED_BY_EMAIL: 'proteomics.published_by_email'`
  - `PUBLISHED_TARGET: 'proteomics.published_target'`
  - `PUBLISHED_DE_METHOD: 'proteomics.published_de_method'`
  - `PUBLISHED_FC_THRESHOLD: 'proteomics.published_fc_threshold'`
  - `PUBLISHED_P_THRESHOLD: 'proteomics.published_p_threshold'`
  - `PUBLISHED_VERSION: 'proteomics.published_version'`
  - `PUBLISHED_ID: 'proteomics.published_id'`
  - `PUBLISHED_INCLUDES_ENRICHMENT: 'proteomics.published_includes_enrichment'`
  - `SUPERSEDED_BY: 'proteomics.superseded_by'`
  - `SUPERSEDES: 'proteomics.supersedes'`

Declare the matching `META_COLUMNS` record (one-to-one with `PUBLISHED_TAGS`, names prefixed `_meta_published_*` per RESEARCH §"Pattern 3" naming). Used by both Plan 02 (writer) and Plan 06 (reader) to avoid string drift.

Also declare `DE_COMPLETE_TAG = 'proteomics.de_complete' as const;` — exported so Plan 05 dialog precondition gate and Plan 07 menu handler do not inline the string (I-10 fix).

Use `const PUBLISHED_TAGS = { ... } as const;` so TypeScript narrows the values to string literals (downstream `as const` consumers don't widen).

Per CLAUDE.md code style: 2-space indent, single quotes, semicolons, kebab-case filename (already correct), TypeScript strict. JSDoc on every exported member (sibling pattern of `experiment-setup.ts`).
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/publish-state.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "export interface PublishedMetadata" src/publishing/publish-state.ts | grep -q '^1$' &amp;&amp; grep -c "export interface PublishOptions" src/publishing/publish-state.ts | grep -q '^1$' &amp;&amp; grep -c "export interface MailtoOptions" src/publishing/publish-state.ts | grep -q '^1$' &amp;&amp; grep -c "PUBLISHED_TAGS" src/publishing/publish-state.ts | grep -qv '^0$' &amp;&amp; grep -c "DE_COMPLETE_TAG" src/publishing/publish-state.ts | grep -qv '^0$' &amp;&amp; grep -q "as const" src/publishing/publish-state.ts</automated>
  </verify>
  <done>`src/publishing/publish-state.ts` exists; exports `PublishedMetadata`, `PublishOptions`, `MailtoOptions`, `PUBLISHED_TAGS`, `META_COLUMNS`, `DE_COMPLETE_TAG`; type-checks cleanly under strict mode; tag values use `proteomics.` prefix and `as const` narrowing.</done>
</task>

<task type="auto" tdd="false">
  <name>Task 2: Implement isPublished + getPublishedMetadata + setPublishedTags helpers</name>
  <files>src/publishing/publish-state.ts</files>
  <read_first>
    - @packages/Proteomics/src/analysis/experiment-setup.ts (lines 49-70 — getGroups/setGroups parse pattern with try/catch)
    - @packages/Proteomics/src/viewers/heatmap.ts (line 32 — `df.getTag('proteomics.de_complete') === 'true'` idiom)
    - @.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md (Section 1 — column-first, tag-second read order)
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"Pattern 3" — belt-and-braces column shape)
  </read_first>
  <action>
Add to `src/publishing/publish-state.ts` (same file as Task 1):

`isPublished(df: DG.DataFrame): boolean` — single line: `return df.getTag(PUBLISHED_TAGS.PUBLISHED) === 'true';`. Mirrors the `proteomics.de_complete === 'true'` idiom in `heatmap.ts:32`.

`getPublishedMetadata(df: DG.DataFrame): PublishedMetadata | null` — returns `null` if `!isPublished(df)`. Otherwise reads each field. Per CONTEXT.md "Belt-and-braces is the design philosophy" + PATTERNS.md Section 6 "COLUMN FIRST, TAG SECOND" — read order:
  - First try `df.col(META_COLUMNS.PUBLISHED_TARGET)?.get(0)` (cast/coerce to string)
  - If null/undefined, fall back to `df.getTag(PUBLISHED_TAGS.PUBLISHED_TARGET)`
  - If still null, return null (corrupted publish — caller surfaces error)

Apply the same column-first, tag-second order for every field. `publishedAt` reads as `Date` from `addNewDateTime` column or ISO `string` from tag; normalize to `Date` in the returned object. `fcThreshold` / `pThreshold` parse as `parseFloat`. `version` parses as `parseInt(..., 10)`. `includesEnrichment` checks `=== 'true'`. `supersedes` / `supersededBy` may be null and that is non-fatal.

Wrap each field-read in a try/catch so a corrupt single field returns `null` for that field rather than crashing the whole read. This matches the `try/catch` JSON parse in `experiment-setup.ts:54-61`.

`setPublishedTags(df: DG.DataFrame, meta: PublishedMetadata): void` — single-purpose: writes every `PUBLISHED_TAGS.*` key as a string tag on the DataFrame. Used by Plan 02 (`trimForPublish`) AFTER cloning. Per Pitfall 3 mitigation: the clone's tag-preservation is partial, so every required tag is re-set explicitly. Convert non-string values to string: `String(meta.fcThreshold)`, `meta.publishedAt.toISOString()` etc. For `supersedes`/`supersededBy` null values, do NOT call setTag (or call `df.setTag(key, '')` — pick one and document; recommend NOT setting null fields so reader can disambiguate "missing" from "empty").
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/publish-state.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "export function isPublished" src/publishing/publish-state.ts | grep -q '^1$' &amp;&amp; grep -c "export function getPublishedMetadata" src/publishing/publish-state.ts | grep -q '^1$' &amp;&amp; grep -c "export function setPublishedTags" src/publishing/publish-state.ts | grep -q '^1$' &amp;&amp; grep -E "df\.col\(META_COLUMNS" src/publishing/publish-state.ts | grep -v '^#'</automated>
  </verify>
  <done>Three exports `isPublished`, `getPublishedMetadata`, `setPublishedTags` present; `getPublishedMetadata` reads metadata column FIRST and falls back to tag SECOND (per Pitfall 3 mitigation); `setPublishedTags` writes every key in `PUBLISHED_TAGS`; strict-mode TypeScript passes.</done>
</task>

<task type="auto" tdd="false">
  <name>Task 3: Implement slugifyTarget + findPriorShare + buildMailtoUrl</name>
  <files>src/publishing/publish-state.ts</files>
  <read_first>
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (D-01 slug rules — preserve case, replace special chars with `-`, collapse repeats, drop trailing punctuation, cap ~64 chars, charset `[A-Za-z0-9._-]`; Claude's discretion mailto body shape)
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"Don't Hand-Roll" row "Slug sanitization" — single regex chain; §"Don't Hand-Roll" row "mailto via encodeURIComponent"; §"Open Question 1" — smart-filter `like` syntax; §"Open Question 2" — version increment from prior+1; spike Task 1 confirms `like` works or falls back; §"Open Question 5" — sharer email null fallback)
    - @packages/Proteomics/.planning/phases/15-read-only-publishing-foundation/15-00-SUMMARY.md (if present — Plan 00 resolved Assumption A8)
  </read_first>
  <action>
Add to `src/publishing/publish-state.ts`:

`slugifyTarget(raw: string): string` — implements D-01 sanitization rules. Sequence:
  1. Replace every char NOT in `[A-Za-z0-9._-]` with `-` (regex: `/[^A-Za-z0-9._-]+/g`)
  2. Collapse repeats: replace `/-{2,}/g` with single `-`
  3. Drop leading/trailing `-` and `.` punctuation: `replace(/^[-.]+|[-.]+$/g, '')`
  4. Cap at 64 chars: `.slice(0, 64)`
  5. After cap, drop trailing `-` or `.` again (cap may have created one)
  6. Return; if result is empty string after sanitization, return `'unnamed'` (defensive — never produce empty slug that would break Project.name)

CASE PRESERVED (D-01 explicit: "preserve case"). Do NOT `.toLowerCase()`.

`findPriorShare(target: string, group: DG.Group | null): Promise<DG.Project | null>` — searches for the most recent published Project matching the slug+group. T-15-05 mitigation: query is scoped to projects the current user can administer (server-side filter handles this implicitly when querying via `dapi.projects.filter`).

Implementation (PRIMARY — use spike Task 1 confirmation of A8):
```
const slug = slugifyTarget(target);
const namePattern = `Proteomics-Review-${slug}-v%`;
const candidates = await grok.dapi.projects.filter(`name like "${namePattern}"`).list();
```

Fallback (if A8 resolved that `like` is NOT supported per spike output): replace with `grok.dapi.projects.list()` and client-side filter via regex `^Proteomics-Review-${slug}-v(\d+)-\d{4}-\d{2}-\d{2}$`.

After getting candidates, if `group` is non-null, filter further to only those whose Project carries the reviewer-group association. (Resolution path depends on spike A2 output — Plan 00 confirmed how Space-inheritance propagates. If `permissions.get(project)` shows the group in `view`, keep. If we cannot reliably query group from candidates client-side, skip group filter and return latest version — reviewer-group is a hint, not a hard predicate, for prefilling.)

Parse version from project name via regex `/-v(\d+)-/`, return the candidate with `max(version)`. Return `null` if no candidates.

Add JSDoc noting the function is purely informational (used by Plan 05 `share-dialog.ts` for republish-detection banner + Plan 04 `publish-project.ts` for supersede-chain wiring) — NOT a security gate.

`buildMailtoUrl(opts: MailtoOptions): string` — pure helper. Per RESEARCH §"Don't Hand-Roll" and CONTEXT.md PUB-13 Claude's discretion:
  - Subject: `'Re-run request: ' + opts.projectName`
  - Body: `'Hi ' + opts.sharerName + ', could you re-run with [different parameters]? Looking at ' + opts.projectName + ' (shared ' + opts.publishedDateStr + ').'`
  - URL: `'mailto:' + (opts.sharerEmail ? encodeURIComponent(opts.sharerEmail) : '') + '?subject=' + encodeURIComponent(subject) + '&body=' + encodeURIComponent(body)`
  - Return string. Pure function — no DOM, no `grok.shell` calls.

**Co-location rationale (B-1 fix):** `buildMailtoUrl` is exported HERE rather than in `share-dialog.ts` because Plan 06 (reviewer panel, wave 2) needs it for PUB-13. share-dialog.ts is wave 4 — a wave-2 import from wave 4 would break the dependency graph. Plan 05 imports `buildMailtoUrl` from this module (`./publish-state`), not the other way around. The function is pure with zero UI/dialog/DOM dependency, so this is the cleanest home.
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/publish-state.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "export function slugifyTarget" src/publishing/publish-state.ts | grep -q '^1$' &amp;&amp; grep -c "export async function findPriorShare" src/publishing/publish-state.ts | grep -q '^1$' &amp;&amp; grep -c "export function buildMailtoUrl" src/publishing/publish-state.ts | grep -q '^1$' &amp;&amp; grep -c "encodeURIComponent" src/publishing/publish-state.ts | grep -qv '^0$' &amp;&amp; grep -c "mailto:" src/publishing/publish-state.ts | grep -qv '^0$'</automated>
  </verify>
  <done>`slugifyTarget` sanitizes per D-01 rules (preserves case, charset `[A-Za-z0-9._-]`, cap 64, never empty); `findPriorShare(target, group)` returns the most-recent prior Project matching the slug or null; `buildMailtoUrl(opts)` builds an `encodeURIComponent`-safe mailto URL; all four exported; strict-mode TypeScript passes.</done>
</task>

</tasks>

<threat_model>
## Trust Boundaries

| Boundary | Description |
|----------|-------------|
| user input (target string) -> slugifyTarget | Freeform user input; could contain shell/SQL/URL injection chars |
| current user's filter scope -> findPriorShare | Returns projects user can administer; never reveals projects user cannot see |
| mailto subject/body inputs -> buildMailtoUrl | Values flow from server-stored metadata (publishedBy, projectName); `encodeURIComponent` prevents header injection |

## STRIDE Threat Register

| Threat ID | Category | Component | Disposition | Mitigation Plan |
|-----------|----------|-----------|-------------|-----------------|
| T-15-03 | Tampering (freeform target injection) | slugifyTarget | mitigate | Strict charset `[A-Za-z0-9._-]` regex strips every other char; cap at 64 chars prevents long-name attacks; raw target preserved in tag/column only, never used in URL/Space name |
| T-15-05 | Information disclosure (findPriorShare scope) | findPriorShare | mitigate | `dapi.projects.filter` is scoped server-side to projects the user can administer (per platform ACL); helper does NOT use `dapi.projects.list()` unless A8 fallback forces it; in fallback path, server still scopes by user permissions |
| T-15-SP-03 | Information disclosure (tag corruption) | getPublishedMetadata | accept | Corrupt single field returns null for that field; caller (Plan 06 panel) renders "Unavailable" instead of crashing; no security impact |
| T-15-SP-06 | Tampering (mailto header injection) | buildMailtoUrl | mitigate | All variable parts wrapped in `encodeURIComponent`; subject/body templates use only sanitized strings; no raw user input bypasses encoding |
</threat_model>

<verification>
- TypeScript strict-mode compiles cleanly
- `PUBLISHED_TAGS` and `META_COLUMNS` records have parallel keys (Plan 02 + Plan 06 depend on parity)
- `slugifyTarget('!@#$ MYH7 -DMD ')` returns a valid `[A-Za-z0-9._-]` string (smoke check via REPL or test)
- `findPriorShare('nonexistent-target', null)` returns `null` against the live server (smoke check)
- `buildMailtoUrl({sharerEmail: 'a@b.com', sharerName: 'A', projectName: 'P', publishedDateStr: '2026-06-07'})` returns a string starting with `mailto:` and containing `encodeURIComponent`-safe values (smoke check)
</verification>

<success_criteria>
- File `src/publishing/publish-state.ts` exists with all 6 exports + 3 typed shapes + 2 const records + DE_COMPLETE_TAG constant
- All 13 tag keys in `PUBLISHED_TAGS` start with `proteomics.`
- `META_COLUMNS` keys match `PUBLISHED_TAGS` keys one-to-one (parity check)
- Every export has JSDoc
- No tag-string literal appears outside this file (Plan 02-08 import the constants)
- `buildMailtoUrl` is pure (no `grok.shell` calls, no DOM) so it can be safely imported by wave-2 callers (Plan 06)
</success_criteria>

<output>
Create `.planning/phases/15-read-only-publishing-foundation/15-01-SUMMARY.md` when done with:
- Exports summary (6 functions + 3 interfaces + 2 const records + 1 string constant)
- Notes on `findPriorShare` filter syntax actually used (PRIMARY `like` or FALLBACK `list()` per spike output)
- Confirmation that no other file inlines `proteomics.published*` strings
- Confirmation that `buildMailtoUrl` lives here (not in share-dialog.ts) per B-1 wave-cycle fix
</output>
