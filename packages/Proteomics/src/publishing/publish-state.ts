import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {reviewNamePrefix} from './publish-settings';

/**
 * Single source of truth for Phase 15's published-tag namespace + pure helpers.
 *
 * Sibling pattern of `src/analysis/experiment-setup.ts`'s `getGroups`/`setGroups`:
 * tag I/O lives in one module, every consumer imports from it, no other file
 * inlines `proteomics.published*` tag string literals. A single typo in a tag
 * string would silently break the reviewer panel; centralizing is the antidote.
 */

/** DE-pipeline completion gate. Exported here so Plan 05 dialog precondition
 *  and Plan 07 menu handler do not inline the string (I-10 fix). */
export const DE_COMPLETE_TAG = 'proteomics.de_complete' as const;

/**
 * Canonical published-tag keyspace. Every tag uses the `proteomics.` prefix.
 * Mirrored 1:1 by {@link META_COLUMNS} below — Plan 02 writes BOTH; Plan 06
 * reads column-first / tag-second (Pitfall 3 belt-and-braces).
 */
export const PUBLISHED_TAGS = {
  PUBLISHED: 'proteomics.published',
  PUBLISHED_AT: 'proteomics.published_at',
  PUBLISHED_BY: 'proteomics.published_by',
  PUBLISHED_BY_EMAIL: 'proteomics.published_by_email',
  PUBLISHED_TARGET: 'proteomics.published_target',
  PUBLISHED_DE_METHOD: 'proteomics.published_de_method',
  PUBLISHED_FC_THRESHOLD: 'proteomics.published_fc_threshold',
  PUBLISHED_P_THRESHOLD: 'proteomics.published_p_threshold',
  PUBLISHED_VERSION: 'proteomics.published_version',
  PUBLISHED_ID: 'proteomics.published_id',
  PUBLISHED_INCLUDES_ENRICHMENT: 'proteomics.published_includes_enrichment',
  SUPERSEDED_BY: 'proteomics.superseded_by',
  SUPERSEDES: 'proteomics.supersedes',
} as const;

/**
 * Belt-and-braces metadata column names. One-to-one with {@link PUBLISHED_TAGS}.
 * Plan 02 writes a single-row column per entry alongside the tag; Plan 06
 * reads column FIRST, tag SECOND so a serialization path that strips tags but
 * preserves columns (or vice versa) still recovers the metadata.
 */
export const META_COLUMNS = {
  PUBLISHED: '_meta_published',
  PUBLISHED_AT: '_meta_published_at',
  PUBLISHED_BY: '_meta_published_by',
  PUBLISHED_BY_EMAIL: '_meta_published_by_email',
  PUBLISHED_TARGET: '_meta_published_target',
  PUBLISHED_DE_METHOD: '_meta_published_de_method',
  PUBLISHED_FC_THRESHOLD: '_meta_published_fc_threshold',
  PUBLISHED_P_THRESHOLD: '_meta_published_p_threshold',
  PUBLISHED_VERSION: '_meta_published_version',
  PUBLISHED_ID: '_meta_published_id',
  PUBLISHED_INCLUDES_ENRICHMENT: '_meta_published_includes_enrichment',
  SUPERSEDED_BY: '_meta_superseded_by',
  SUPERSEDES: '_meta_supersedes',
} as const;

/** Typed view of one published-analysis DataFrame's metadata. */
export interface PublishedMetadata {
  /** Raw user input — never use directly for paths/URLs (use {@link slugifyTarget}). */
  target: string;
  /** Normalized to Date when read; ISO string accepted on write. */
  publishedAt: Date;
  publishedBy: string;
  publishedByEmail: string | null;
  deMethod: 'limma' | 'deqms' | 't-test' | 'spectronaut' | string;
  fcThreshold: number;
  pThreshold: number;
  /** Monotonically increasing per (target, group); first share is version 1. */
  version: number;
  /** The published Project's id — set after `projects.save` returns. */
  publishId: string;
  includesEnrichment: boolean;
  /** Prior Project id; null on the first share. */
  supersedes: string | null;
  /** Next Project id; set retroactively on the prior version when republishing. */
  supersededBy: string | null;
}

/** Input shape consumed by `share-dialog.ts` → `publishAnalysis(df, opts)`. */
export interface PublishOptions {
  target: string;
  reviewerGroup: DG.Group;
  note: string;
  /** Set by {@link findPriorShare} when this is a republish; null on first share. */
  priorVersion: DG.Project | null;
  /**
   * Test-only override for the umbrella Space name.
   * Defaults to `'Proteomics-Reviews'`. Production callers never set this;
   * Plan 08 Test 5 uses it to inject a throwaway umbrella.
   */
  umbrellaName?: string;
  /**
   * Whether to run the reopen-and-verify round trip after publishing. Set by the
   * share dialog's "Verify published dashboard" checkbox (which defaults to the
   * `verifyPublishedDashboard` package setting). When undefined (non-dialog callers),
   * `publishAnalysis` falls back to the package setting.
   */
  verify?: boolean;
}

/** Input shape for {@link buildMailtoUrl}. Pure data — no DOM. */
export interface MailtoOptions {
  sharerEmail: string | null;
  sharerName: string;
  projectName: string;
  publishedDateStr: string;
}

/** Returns true iff `proteomics.published === 'true'`. Mirror of the
 *  `proteomics.de_complete === 'true'` idiom in `viewers/heatmap.ts`. */
export function isPublished(df: DG.DataFrame): boolean {
  return df.getTag(PUBLISHED_TAGS.PUBLISHED) === 'true';
}

/**
 * Reads the published metadata from a DataFrame, column-FIRST and tag-SECOND
 * (Pitfall 3 belt-and-braces). Returns `null` when the DataFrame is not a
 * published analysis. Returns a populated object when at least the required
 * load-bearing fields (target, publishId) are recoverable; per-field corruption
 * falls back gracefully without crashing the whole read.
 */
export function getPublishedMetadata(df: DG.DataFrame): PublishedMetadata | null {
  if (!isPublished(df)) return null;

  const readString = (colKey: keyof typeof META_COLUMNS, tagKey: keyof typeof PUBLISHED_TAGS): string | null => {
    try {
      const col = df.col(META_COLUMNS[colKey]);
      if (col != null) {
        const v = col.get(0);
        if (v != null && v !== '') return String(v);
      }
    } catch { /* fall through to tag */ }
    try {
      return df.getTag(PUBLISHED_TAGS[tagKey]) ?? null;
    } catch {
      return null;
    }
  };

  const readDate = (): Date | null => {
    try {
      const col = df.col(META_COLUMNS.PUBLISHED_AT);
      if (col != null) {
        const v = col.get(0);
        if (v instanceof Date) return v;
        if (typeof v === 'string' && v) return new Date(v);
        if (typeof v === 'number') return new Date(v);
      }
    } catch { /* fall through */ }
    try {
      const raw = df.getTag(PUBLISHED_TAGS.PUBLISHED_AT);
      if (raw) return new Date(raw);
    } catch { /* return null */ }
    return null;
  };

  const readNumber = (colKey: keyof typeof META_COLUMNS, tagKey: keyof typeof PUBLISHED_TAGS): number => {
    const s = readString(colKey, tagKey);
    if (s == null) return NaN;
    const n = parseFloat(s);
    return Number.isFinite(n) ? n : NaN;
  };

  const readBool = (colKey: keyof typeof META_COLUMNS, tagKey: keyof typeof PUBLISHED_TAGS): boolean => {
    const s = readString(colKey, tagKey);
    return s === 'true';
  };

  const target = readString('PUBLISHED_TARGET', 'PUBLISHED_TARGET');
  const publishId = readString('PUBLISHED_ID', 'PUBLISHED_ID');
  if (target == null || publishId == null) return null;

  const publishedAt = readDate();
  const versionRaw = readString('PUBLISHED_VERSION', 'PUBLISHED_VERSION');
  const versionNum = versionRaw != null ? parseInt(versionRaw, 10) : NaN;

  return {
    target,
    publishedAt: publishedAt ?? new Date(NaN),
    publishedBy: readString('PUBLISHED_BY', 'PUBLISHED_BY') ?? '',
    publishedByEmail: readString('PUBLISHED_BY_EMAIL', 'PUBLISHED_BY_EMAIL'),
    deMethod: readString('PUBLISHED_DE_METHOD', 'PUBLISHED_DE_METHOD') ?? '',
    fcThreshold: readNumber('PUBLISHED_FC_THRESHOLD', 'PUBLISHED_FC_THRESHOLD'),
    pThreshold: readNumber('PUBLISHED_P_THRESHOLD', 'PUBLISHED_P_THRESHOLD'),
    version: Number.isFinite(versionNum) ? versionNum : 1,
    publishId,
    includesEnrichment: readBool('PUBLISHED_INCLUDES_ENRICHMENT', 'PUBLISHED_INCLUDES_ENRICHMENT'),
    supersedes: readString('SUPERSEDES', 'SUPERSEDES'),
    supersededBy: readString('SUPERSEDED_BY', 'SUPERSEDED_BY'),
  };
}

/**
 * Writes every {@link PUBLISHED_TAGS} entry as a string tag on the DataFrame.
 * Pitfall 3 mitigation: even though spike 15-00 showed all 14 tags survive the
 * basic save→find→open path, `df.clone()` may not carry the full tag map;
 * Plan 02 calls this AFTER cloning to guarantee a complete tag namespace.
 *
 * Null `supersedes` / `supersededBy` are intentionally NOT written so reader
 * can disambiguate "missing" from "empty". The PUBLISHED tag is always set to
 * 'true' regardless of the metadata input.
 */
export function setPublishedTags(df: DG.DataFrame, meta: PublishedMetadata): void {
  df.setTag(PUBLISHED_TAGS.PUBLISHED, 'true');
  df.setTag(PUBLISHED_TAGS.PUBLISHED_AT,
    meta.publishedAt instanceof Date ? meta.publishedAt.toISOString() : String(meta.publishedAt));
  df.setTag(PUBLISHED_TAGS.PUBLISHED_BY, meta.publishedBy);
  if (meta.publishedByEmail != null)
    df.setTag(PUBLISHED_TAGS.PUBLISHED_BY_EMAIL, meta.publishedByEmail);
  df.setTag(PUBLISHED_TAGS.PUBLISHED_TARGET, meta.target);
  df.setTag(PUBLISHED_TAGS.PUBLISHED_DE_METHOD, meta.deMethod);
  df.setTag(PUBLISHED_TAGS.PUBLISHED_FC_THRESHOLD, String(meta.fcThreshold));
  df.setTag(PUBLISHED_TAGS.PUBLISHED_P_THRESHOLD, String(meta.pThreshold));
  df.setTag(PUBLISHED_TAGS.PUBLISHED_VERSION, String(meta.version));
  df.setTag(PUBLISHED_TAGS.PUBLISHED_ID, meta.publishId);
  df.setTag(PUBLISHED_TAGS.PUBLISHED_INCLUDES_ENRICHMENT, meta.includesEnrichment ? 'true' : 'false');
  if (meta.supersedes != null)
    df.setTag(PUBLISHED_TAGS.SUPERSEDES, meta.supersedes);
  if (meta.supersededBy != null)
    df.setTag(PUBLISHED_TAGS.SUPERSEDED_BY, meta.supersededBy);
}

/**
 * Sanitizes a freeform target string per D-01 into a safe slug.
 *
 * Rules:
 *  - Charset `[A-Za-z0-9._-]` (case PRESERVED — do NOT lowercase)
 *  - Replace runs of disallowed chars with single `-`
 *  - Drop leading/trailing `-` or `.`
 *  - Cap at 64 chars
 *  - Empty result → `'unnamed'` (never produce empty slug — would break Project.name)
 */
export function slugifyTarget(raw: string): string {
  if (raw == null) return 'unnamed';
  let s = String(raw)
    .replace(/[^A-Za-z0-9._-]+/g, '-')
    .replace(/-{2,}/g, '-')
    .replace(/^[-.]+|[-.]+$/g, '')
    .slice(0, 64)
    .replace(/[-.]+$/g, '');
  if (s.length === 0) s = 'unnamed';
  return s;
}

/**
 * Looks up the most recent prior published Project matching the slugified
 * target. Used by Plan 05 (`share-dialog.ts`) for republish-detection banner
 * and Plan 04 (`publish-project.ts`) for the supersede chain.
 *
 * Smart-filter `like` confirmed working by spike 15-00 (assumption A8); the
 * fallback `list()` + client-side filter path is retained as defense in depth
 * in case the platform smart-filter parser changes.
 *
 * NOT a security gate — purely informational. T-15-05 mitigation relies on
 * the platform ACL implicitly scoping `dapi.projects.filter` to projects the
 * current user can administer.
 */
export async function findPriorShare(target: string, _group: DG.Group | null): Promise<DG.Project | null> {
  const slug = slugifyTarget(target);
  const namePrefix = `${reviewNamePrefix()}-${slug}-v`;
  const versionRe = /-v(\d+)-/;

  const pickLatest = (cands: DG.Project[]): DG.Project | null => {
    let best: DG.Project | null = null;
    let bestV = -1;
    for (const p of cands) {
      const name = (p as any)?.name ?? '';
      const m = versionRe.exec(name);
      if (!m) continue;
      const v = parseInt(m[1], 10);
      if (Number.isFinite(v) && v > bestV) {
        bestV = v;
        best = p;
      }
    }
    return best;
  };

  try {
    const namePattern = `${namePrefix}%`;
    const candidates = await grok.dapi.projects.filter(`name like "${namePattern}"`).list();
    if (Array.isArray(candidates) && candidates.length > 0)
      return pickLatest(candidates);
  } catch { /* fall through to client-side filter */ }

  try {
    const all = await grok.dapi.projects.list();
    const matching = (Array.isArray(all) ? all : [])
      .filter((p) => ((p as any)?.name ?? '').startsWith(namePrefix));
    if (matching.length > 0) return pickLatest(matching);
  } catch { /* return null */ }

  return null;
}

/**
 * Pure helper for the reviewer panel's "Request re-run" button (PUB-13).
 * Builds an `encodeURIComponent`-safe `mailto:` URL with the sharer's email,
 * a project-aware subject, and a body that names the project + share date.
 *
 * Co-located here rather than in `share-dialog.ts` because Plan 06 (wave 2)
 * needs it BEFORE Plan 05 (wave 4) lands — importing wave 4 from wave 2 would
 * break the wave graph. No DOM / UI / dialog dependency.
 */
export function buildMailtoUrl(opts: MailtoOptions): string {
  const subject = `Re-run request: ${opts.projectName}`;
  const body = `Hi ${opts.sharerName}, could you re-run with [different parameters]? ` +
    `Looking at ${opts.projectName} (shared ${opts.publishedDateStr}).`;
  const to = opts.sharerEmail ? encodeURIComponent(opts.sharerEmail) : '';
  return `mailto:${to}?subject=${encodeURIComponent(subject)}&body=${encodeURIComponent(body)}`;
}
