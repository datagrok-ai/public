import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, delay} from '@datagrok-libraries/test/src/test';

import {SEMTYPE} from '../utils/proteomics-types';
import {findColumn} from '../utils/column-detection';
import {
  META_COLUMNS, PUBLISHED_TAGS, PublishedMetadata, getPublishedMetadata,
} from './publish-state';

/**
 * Load-bearing distinguishing marker. Plan 04 step 8's catch block matches on
 * this prefix to identify the "look/filter config stripped by serializer"
 * failure mode and trigger the post-open formula-line recovery hook (B-2 / W-7)
 * BEFORE rolling back the published Project.
 */
export const FORMULA_LINE_ASSERTION_PREFIX = 'FORMULA_LINE_MISSING:' as const;

/** What the publisher claimed at publish time, used as the expected baseline
 *  for {@link assertPublishedShape}. */
export interface PublishedShapeContract {
  /** The trimmed DF's name (e.g. `<source>_published_<YYYY-MM-DD>`). */
  expectedName: string;
  /** The Project's name (e.g. `Proteomics-Review-<slug>-v<N>-<YYYY-MM-DD>`). */
  expectedProjectName: string;
  /** The 7-or-fewer column names the trim chose. */
  expectedAllowlist: string[];
  expectedMeta: PublishedMetadata;
  /** Whether a volcano viewer should be present on reopen. */
  expectVolcano: boolean;
  /** Whether an enrichment DataFrame should be present in `grok.shell.tables`. */
  expectEnrichment: boolean;
  /** Whether FC/p threshold formula lines should be present on the volcano.
   *  PUB-06 / SC-2 — defaults to true when `expectVolcano` is true. */
  expectFormulaLines: boolean;
}

/**
 * Reopens the saved Project in a fresh session and asserts every observable
 * invariant of the publish contract. Throws on first violation with a
 * specific error including which assertion failed and what was observed.
 *
 * Single source of truth: BOTH the Plan 04 orchestrator (step 8 non-negotiable
 * gate) AND the Plan 08 regression test (`src/tests/publish-roundtrip.ts`)
 * call this helper. If the contract changes, it changes here.
 *
 * The volcano-formula-line assertion (8a) throws an error PREFIXED with
 * {@link FORMULA_LINE_ASSERTION_PREFIX}. Plan 04 step 8 catches this specific
 * prefix to drive its self-healing recovery (W-7) — re-add the formula lines
 * from `proteomics.published_fc_threshold` / `proteomics.published_p_threshold`
 * tags and re-run this assertion exactly once before rolling back.
 *
 * Caller is responsible for `grok.shell.closeAll()` after the assertion if
 * additional cleanup is desired.
 */
export async function assertPublishedShape(
  project: DG.Project,
  contract: PublishedShapeContract,
): Promise<void> {
  const projId = project.id;
  grok.shell.closeAll();
  await delay(100);

  const reopened = await grok.dapi.projects.find(projId);
  await reopened.open();

  await awaitCheck(
    () => !!grok.shell.tv && !!grok.shell.tv.dataFrame,
    'assertPublishedShape: TableView never materialized after project.open()',
    5000,
  );

  const reDf = grok.shell.tv!.dataFrame;

  // ASSERTION 1 — df.name matches
  if (reDf.name !== contract.expectedName) {
    throw new Error(
      `assertPublishedShape: df.name: got "${reDf.name}", expected "${contract.expectedName}"`);
  }

  // ASSERTION 2 — Project name matches (PUB-09)
  if ((reopened as any).name !== contract.expectedProjectName) {
    throw new Error(
      `assertPublishedShape: project.name: got "${(reopened as any).name}", expected "${contract.expectedProjectName}"`);
  }

  // ASSERTION 3 — Every allowlist column present
  for (const colName of contract.expectedAllowlist) {
    if (!reDf.col(colName)) {
      throw new Error(
        `assertPublishedShape: allowlist column missing on reopen: "${colName}"`);
    }
  }

  // ASSERTION 4 — Column count matches allowlist (no extras except hidden
  // _meta_* belt-and-braces columns and ~-prefixed technical columns such as
  // the volcano's '~Volcano label' — neither is part of the published shape).
  const visibleCols = reDf.columns.toList()
    .map((c) => c.name)
    .filter((n) => !n.startsWith('_meta_') && !n.startsWith('~'));
  if (visibleCols.length !== contract.expectedAllowlist.length) {
    const unexpected = visibleCols.filter((n) => !contract.expectedAllowlist.includes(n));
    const missing = contract.expectedAllowlist.filter((n) => !visibleCols.includes(n));
    throw new Error(
      `assertPublishedShape: visible column count: got ${visibleCols.length}, expected ${contract.expectedAllowlist.length}. ` +
      `Unexpected: [${unexpected.join(', ')}]. Missing: [${missing.join(', ')}].`);
  }

  // ASSERTION 5 — PROTEIN_ID semType present on protein-id column
  const proteinIdCol = findColumn(reDf, SEMTYPE.PROTEIN_ID, ['protein id', 'majority protein id', 'accession']);
  if (proteinIdCol == null) {
    throw new Error('assertPublishedShape: protein-id column not found via findColumn');
  }
  if (proteinIdCol.semType !== SEMTYPE.PROTEIN_ID) {
    throw new Error(
      `assertPublishedShape: protein-id semType: got "${proteinIdCol.semType ?? 'null'}", ` +
      `expected "${SEMTYPE.PROTEIN_ID}". Detectors may not have re-fired on reopen — ` +
      `re-assign in trim-dataframe.ts Step E.`);
  }

  // ASSERTION 6 — Required proteomics.published* metadata readable (column FIRST, tag SECOND)
  const meta = getPublishedMetadata(reDf);
  if (meta == null) {
    throw new Error(
      'assertPublishedShape: getPublishedMetadata returned null — published flag ' +
      'not set or all critical fields unreadable');
  }
  const expectedMeta = contract.expectedMeta;
  const eqFloat = (a: number, b: number) => !(Math.abs(a - b) > 1e-9);
  if (meta.target !== expectedMeta.target)
    throw new Error(`assertPublishedShape: meta.target: got "${meta.target}", expected "${expectedMeta.target}"`);
  if (meta.deMethod !== expectedMeta.deMethod)
    throw new Error(`assertPublishedShape: meta.deMethod: got "${meta.deMethod}", expected "${expectedMeta.deMethod}"`);
  if (!eqFloat(meta.fcThreshold, expectedMeta.fcThreshold))
    throw new Error(`assertPublishedShape: meta.fcThreshold: got ${meta.fcThreshold}, expected ${expectedMeta.fcThreshold}`);
  if (!eqFloat(meta.pThreshold, expectedMeta.pThreshold))
    throw new Error(`assertPublishedShape: meta.pThreshold: got ${meta.pThreshold}, expected ${expectedMeta.pThreshold}`);
  if (meta.version !== expectedMeta.version)
    throw new Error(`assertPublishedShape: meta.version: got ${meta.version}, expected ${expectedMeta.version}`);
  if (meta.publishId !== expectedMeta.publishId)
    throw new Error(`assertPublishedShape: meta.publishId: got "${meta.publishId}", expected "${expectedMeta.publishId}"`);
  if (meta.includesEnrichment !== expectedMeta.includesEnrichment)
    throw new Error(`assertPublishedShape: meta.includesEnrichment: got ${meta.includesEnrichment}, expected ${expectedMeta.includesEnrichment}`);
  if (meta.publishedBy !== expectedMeta.publishedBy)
    throw new Error(`assertPublishedShape: meta.publishedBy: got "${meta.publishedBy}", expected "${expectedMeta.publishedBy}"`);

  // ASSERTION 7 — Belt-and-braces metadata column present (PUB-11)
  for (const key of Object.keys(META_COLUMNS) as Array<keyof typeof META_COLUMNS>) {
    const colName = META_COLUMNS[key];
    if (!reDf.col(colName)) {
      throw new Error(`assertPublishedShape: metadata column missing: "${colName}"`);
    }
  }

  // ASSERTION 8 — Volcano viewer present (PUB-06)
  let volcanoViewer: DG.Viewer | null = null;
  if (contract.expectVolcano) {
    try {
      await awaitCheck(
        () => {
          const tv = grok.shell.tv;
          if (!tv) return false;
          for (const v of tv.viewers) {
            if (v.type === DG.VIEWER.SCATTER_PLOT) {
              volcanoViewer = v;
              return true;
            }
          }
          return false;
        },
        'assertPublishedShape: volcano scatter plot not present on reopen — consider post-open recomputeVolcano in publish-project.ts',
        5000,
      );
    } catch (e) {
      throw new Error(`assertPublishedShape: volcano viewer assertion failed: ${(e as Error)?.message ?? e}`);
    }
  }

  // ASSERTION 8a — Volcano formula lines for FC + p thresholds (PUB-06 / SC-2, B-2 / W-7)
  if (contract.expectVolcano && contract.expectFormulaLines && volcanoViewer != null) {
    const fcStr = String(expectedMeta.fcThreshold);
    const pStr = String(expectedMeta.pThreshold);
    const minusLogPStr = String(-Math.log10(expectedMeta.pThreshold));

    type FormulaLine = {formula?: string; title?: string};
    const collected: FormulaLine[] = [];

    try {
      const dfLines = ((reDf as any).meta?.formulaLines?.items ?? []) as FormulaLine[];
      if (Array.isArray(dfLines)) collected.push(...dfLines);
    } catch { /* fall through */ }
    try {
      const opts: any = (volcanoViewer as DG.Viewer).getOptions?.();
      const lookLines = opts?.look?.formulaLines;
      const parsed = typeof lookLines === 'string' ? JSON.parse(lookLines) : lookLines;
      if (Array.isArray(parsed)) collected.push(...(parsed as FormulaLine[]));
    } catch { /* fall through */ }
    try {
      const propLines = (volcanoViewer as any)?.props?.formulaLines;
      const parsed = typeof propLines === 'string' ? JSON.parse(propLines) : propLines;
      if (Array.isArray(parsed)) collected.push(...(parsed as FormulaLine[]));
    } catch { /* fall through */ }

    const formulas = collected
      .map((l) => (typeof l?.formula === 'string' ? l.formula : ''))
      .filter((f) => f.length > 0);

    const hasFcLine = formulas.some((f) => f.includes(fcStr) || f.includes(`-${fcStr}`));
    const hasPLine = formulas.some((f) => f.includes(minusLogPStr) || f.includes(pStr));

    if (!hasFcLine || !hasPLine) {
      throw new Error(
        `${FORMULA_LINE_ASSERTION_PREFIX} volcano reopened without formula lines for ` +
        `FC=${expectedMeta.fcThreshold} and/or p=${expectedMeta.pThreshold}. ` +
        `The look/filter config was stripped by the serializer (Phase 13 e527d07ba1 evidence). ` +
        `Plan 04 step 8 catch should invoke the post-open recovery hook from tags ` +
        `${PUBLISHED_TAGS.PUBLISHED_FC_THRESHOLD} / ${PUBLISHED_TAGS.PUBLISHED_P_THRESHOLD} ` +
        `and re-run this assertion once.`);
    }
  }

  // ASSERTION 9 — Enrichment DF presence (PUB-12, conditional)
  if (contract.expectEnrichment) {
    const tables = (grok.shell as any).tables as DG.DataFrame[] | undefined;
    if (!Array.isArray(tables) || tables.length < 2) {
      throw new Error(
        `assertPublishedShape: expectEnrichment=true but grok.shell.tables has ${tables?.length ?? 0} DataFrames (need ≥ 2)`);
    }
    const proteinDf = tables.find((d) => {
      try {
        const col = d.col(META_COLUMNS.PUBLISHED_INCLUDES_ENRICHMENT);
        return col != null && String(col.get(0)) === 'true';
      } catch { return false; }
    });
    if (!proteinDf)
      throw new Error('assertPublishedShape: protein DF (with _meta_published_includes_enrichment=true) not found among reopened tables');
    const enrichDf = tables.find((d) => d.getTag('proteomics.enrichment') === 'true');
    if (!enrichDf)
      throw new Error('assertPublishedShape: enrichment DF (proteomics.enrichment=true tag) not found among reopened tables');
  }

  // ASSERTION 10 — Project options bidirectional supersede pointers (PUB-10, conditional)
  const opts: any = (reopened as any).options ?? {};
  if (expectedMeta.supersedes != null) {
    const got = opts[PUBLISHED_TAGS.SUPERSEDES];
    if (got !== expectedMeta.supersedes) {
      throw new Error(
        `assertPublishedShape: project.options['${PUBLISHED_TAGS.SUPERSEDES}']: got "${got ?? 'null'}", ` +
        `expected "${expectedMeta.supersedes}"`);
    }
  }
  if (expectedMeta.supersededBy != null) {
    const got = opts[PUBLISHED_TAGS.SUPERSEDED_BY];
    if (got !== expectedMeta.supersededBy) {
      throw new Error(
        `assertPublishedShape: project.options['${PUBLISHED_TAGS.SUPERSEDED_BY}']: got "${got ?? 'null'}", ` +
        `expected "${expectedMeta.supersededBy}"`);
    }
  }
}
