import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, delay} from '@datagrok-libraries/test/src/test';

import {
  META_COLUMNS, PUBLISHED_TAGS, PublishOptions, PublishedMetadata,
  slugifyTarget,
} from './publish-state';
import {trimEnrichmentForPublish, trimForPublish} from './trim-dataframe';
import {
  FORMULA_LINE_ASSERTION_PREFIX, PublishedShapeContract, assertPublishedShape,
} from './assert-published-shape';
import {createVolcanoPlot} from '../viewers/volcano';
import {dockEnrichmentCharts, openEnrichmentVisualization} from '../viewers/enrichment-viewers';
import {reviewSpaceName, reviewNamePrefix, verifyPublishedDashboard} from './publish-settings';
import {DEFAULT_FC_THRESHOLD, DEFAULT_P_THRESHOLD} from '../utils/proteomics-types';

/**
 * Pure helper: applies (or replaces) the volcano threshold formula lines on
 * the DataFrame attached to {@link viewer}. Writes to `df.meta.formulaLines`
 * since that is where `applyThresholdLines` in `src/viewers/volcano.ts` stores
 * them (the platform's canonical location for scatter-plot formula lines).
 *
 * Exported so Plan 07's post-open recovery hook can re-apply the lines on a
 * reopened published Project's volcano when the serializer strips look config
 * (Phase 13 e527d07ba1 evidence). Idempotent — removes existing FC/y-axis
 * threshold lines before re-adding so repeat calls do not stack.
 */
export function applyVolcanoFormulaLines(
  viewer: DG.Viewer,
  fcThreshold: number,
  pThreshold: number,
): void {
  const df = (viewer as any).dataFrame as DG.DataFrame | undefined;
  if (!df) return;

  let yColName: string | null = null;
  try {
    const opts: any = (viewer as any).getOptions?.();
    yColName = opts?.look?.yColumnName ?? null;
  } catch { /* fall through */ }
  if (!yColName) {
    try { yColName = (viewer as any).props?.yColumnName ?? null; } catch { /* fall through */ }
  }
  if (!yColName) return;

  const yPrefix = `\${${yColName}}`;
  const fcPrefix = '${log2FC}';
  const hLine = -Math.log10(pThreshold);

  const meta: any = (df as any).meta;
  const formulaLines: any = meta?.formulaLines;
  if (!formulaLines) return;

  try {
    const items: any[] = Array.isArray(formulaLines.items) ? formulaLines.items : [];
    formulaLines.items = items.filter((line: any) => {
      const f = line?.formula ?? '';
      return !(typeof f === 'string' && (f.startsWith(yPrefix) || f.startsWith(fcPrefix)));
    });
  } catch { /* best effort */ }

  try {
    formulaLines.addLine({formula: `${yPrefix} = ${hLine}`, color: '#888888', width: 1, visible: true});
    formulaLines.addLine({formula: `${fcPrefix} = ${fcThreshold}`, color: '#888888', width: 1, visible: true});
    formulaLines.addLine({formula: `${fcPrefix} = ${-fcThreshold}`, color: '#888888', width: 1, visible: true});
  } catch { /* best effort */ }

  try { (viewer as any).props.showViewerFormulaLines = true; } catch { /* best effort */ }
}

function parsePriorVersion(priorName: string | null | undefined): number {
  if (!priorName) return 0;
  const m = /-v(\d+)-/.exec(priorName);
  if (!m) return 0;
  const n = parseInt(m[1], 10);
  return Number.isFinite(n) ? n : 0;
}

function generateUuid(): string {
  try {
    const c = (globalThis as any)?.crypto;
    if (c && typeof c.randomUUID === 'function') return c.randomUUID();
  } catch { /* fall through */ }
  // RFC4122 v4 fallback
  return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, (c) => {
    const r = (Math.random() * 16) | 0;
    const v = c === 'x' ? r : (r & 0x3) | 0x8;
    return v.toString(16);
  });
}

/**
 * Orchestrates the read-only publish flow. Source DataFrame is NEVER mutated
 * (Pitfall 1 / T-15-02 guarantee — Plan 02's `trimForPublish` owns the clone
 * boundary). Two non-negotiable gates: Step 7 verifies reviewer-group ACL via
 * `permissions.get` and rolls back on Edit-leak; Step 8 calls
 * {@link assertPublishedShape} against the reopened Project and rolls back on
 * contract violation — with a one-shot W-7 self-heal pass for the
 * formula-line-stripped failure mode.
 *
 * Publish-side belt-and-braces for PUB-06 / SC-2 formula lines:
 *   #1 (this orchestrator, step 6.5) — apply formula lines on the trimmed
 *      volcano BEFORE saving the view
 *   #2 (Plan 01 setPublishedTags via Plan 02 trimForPublish) — persist FC and
 *      p threshold values as `proteomics.published_fc_threshold` /
 *      `proteomics.published_p_threshold` tags on the trimmed DataFrame so
 *      the recovery hook can re-apply them later
 *   #3 (Plan 07 post-open recovery hook) — re-apply formula lines from the
 *      tags on reopen if the serializer stripped them; reuses the exported
 *      {@link applyVolcanoFormulaLines} helper
 *
 * Per spike 15-00 (`15-00-SUMMARY.md`):
 *   - A1: `permissions.get` returns `["view", "edit"]` keys only (no share/delete)
 *   - A2: Space-inheritance grants do NOT propagate to `permissions.get(project)`
 *     → Step 7 grants View directly on the Project, not via the parent Space
 *   - A4: `project.options[*]` survives round-trip → supersede pointer uses options
 */
export async function publishAnalysis(df: DG.DataFrame, opts: PublishOptions): Promise<DG.Project> {
  const pi = DG.TaskBarProgressIndicator.create('Publishing snapshot...');
  try {
    // ─── Step 1 — assemble PublishedMetadata (server-authoritative identity) ────
    pi.description = 'Preparing metadata...';
    const publishedBy = (grok.shell.user as any)?.friendlyName ?? '';
    const publishedByEmail = (grok.shell.user as any)?.email ?? null;
    const deMethod = (df.getTag('proteomics.de_method') ?? 't-test') as PublishedMetadata['deMethod'];

    let fcThreshold = DEFAULT_FC_THRESHOLD;
    let pThreshold = DEFAULT_P_THRESHOLD;
    try {
      const optsAny = opts as any;
      if (typeof optsAny.fcThreshold === 'number' && Number.isFinite(optsAny.fcThreshold))
        fcThreshold = optsAny.fcThreshold;
      if (typeof optsAny.pThreshold === 'number' && Number.isFinite(optsAny.pThreshold))
        pThreshold = optsAny.pThreshold;
    } catch { /* defaults */ }

    const slug = slugifyTarget(opts.target);
    const priorVersionN = parsePriorVersion((opts.priorVersion as any)?.name);
    const version = priorVersionN > 0 ? priorVersionN + 1 : 1;
    const publishId = generateUuid();
    const publishedAt = new Date();
    const dateStr = publishedAt.toISOString().slice(0, 10);

    const enrichSource: DG.DataFrame | null = (() => {
      const tables = (grok.shell as any).tables as DG.DataFrame[] | undefined;
      if (!Array.isArray(tables)) return null;
      return tables.find((t) => t.getTag('proteomics.enrichment') === 'true') ?? null;
    })();

    const meta: PublishedMetadata = {
      target: opts.target,
      publishedAt,
      publishedBy,
      publishedByEmail,
      deMethod,
      fcThreshold,
      pThreshold,
      version,
      publishId,
      includesEnrichment: !!enrichSource,
      supersedes: opts.priorVersion?.id ?? null,
      supersededBy: null,
    };

    // ─── Step 2 — trim protein DF (Pitfall 1: source not mutated) ────────────────
    pi.description = 'Trimming snapshot...';
    const frozen = trimForPublish(df, meta);
    const expectedName = frozen.name;

    // ─── Step 3 — opportunistic enrichment trim (D-05) ──────────────────────────
    let frozenEnrich: DG.DataFrame | null = null;
    if (enrichSource) {
      pi.description = 'Trimming enrichment...';
      frozenEnrich = trimEnrichmentForPublish(enrichSource, meta);
    }

    // ─── Step 4 — ensure umbrella Space (opts.umbrellaName override supported) ──
    pi.description = 'Ensuring umbrella Space...';
    const umbrellaName = opts.umbrellaName ?? reviewSpaceName();
    const dapiAny = grok.dapi as any;
    let umbrella: any = null;

    // Try create-first; if the Space already exists, fall back to enumeration.
    // The Datagrok smart-filter API does not reliably surface Space-flagged
    // Projects via `dapi.projects.filter(... and isSpace = true).first()`,
    // so we rely on `spaces.list()` enumeration as the lookup path.
    try {
      umbrella = await dapiAny.spaces.createRootSpace(umbrellaName);
    } catch (createErr: any) {
      const msg = String(createErr?.message ?? createErr ?? '');
      if (!/already exists|exists/i.test(msg)) throw createErr;
      try {
        const all: any[] = await dapiAny.spaces.list();
        umbrella = (all ?? []).find((p) =>
          p?.name === umbrellaName || (p as any)?.friendlyName === umbrellaName) ?? null;
      } catch { /* fall through */ }
      if (!umbrella) {
        throw new Error(
          `Umbrella Space '${umbrellaName}' exists on the server but could not be resolved ` +
          `by name from spaces.list() (msg: ${msg}). Ask an admin to inspect the conflicting Space.`);
      }
    }
    const umbrellaClient = dapiAny.spaces.id(umbrella.id);

    // ─── Step 5 — ensure per-target child Space ─────────────────────────────────
    pi.description = 'Ensuring per-target Space...';
    const childName = `${reviewNamePrefix()}-${slug}`;
    let childSpace: any = null;

    const findChildByName = async (): Promise<any> => {
      try {
        const kids: any[] = await umbrellaClient.children.filter('Project', false).list();
        return (kids ?? []).find((k) =>
          k?.friendlyName === childName || k?.name === childName) ?? null;
      } catch { return null; }
    };

    try {
      if (await umbrellaClient.subspaceExists(childName)) {
        childSpace = await findChildByName();
      }
    } catch { /* fall through to addSubspace */ }

    if (!childSpace) {
      try {
        childSpace = await umbrellaClient.addSubspace(childName);
      } catch (subErr: any) {
        const msg = String(subErr?.message ?? subErr ?? '');
        if (!/already exists|exists/i.test(msg)) throw subErr;
        childSpace = await findChildByName();
        if (!childSpace) {
          throw new Error(
            `Per-target child Space '${childName}' exists but could not be resolved via ` +
            `umbrella children.list() (msg: ${msg}).`);
        }
      }
    }
    const childClient = dapiAny.spaces.id(childSpace.id);

    // ─── Step 6 — create Project + addChild ─────────────────────────────────────
    pi.description = 'Saving Project...';
    const project = DG.Project.create();
    const projectName = `${reviewNamePrefix()}-${slug}-v${version}-${dateStr}`;
    (project as any).name = projectName;
    // Set friendlyName explicitly so the platform doesn't "humanize" the slug
    // (e.g. DMD → "DM D"). Use the original target verbatim for a clean label,
    // and derive the label's leading words from the (configurable) name prefix.
    (project as any).friendlyName =
      `${reviewNamePrefix().replace(/-/g, ' ')} ${opts.target} v${version} ${dateStr}`;

    try { (project as any).options[PUBLISHED_TAGS.PUBLISHED_ID] = publishId; } catch { /* swallow */ }
    // Step 9 writes the supersede pointer on the prior project; the NEW
    // project's `supersedes` pointer must be written BEFORE save so the
    // round-trip assertion in Step 8 sees it on reopen.
    if (opts.priorVersion != null) {
      try {
        (project as any).options[PUBLISHED_TAGS.SUPERSEDES] = (opts.priorVersion as any).id;
      } catch { /* swallow */ }
    }

    project.addChild(frozen.getTableInfo());
    if (frozenEnrich) project.addChild(frozenEnrich.getTableInfo());

    // ─── Step 6.5 — REVISION B-2: re-create volcano on trimmed DF + apply lines ─
    pi.description = 'Re-rendering volcano on trimmed snapshot...';
    const trimmedTv = grok.shell.addTableView(frozen);
    await delay(100);
    const volcano = createVolcanoPlot(frozen, {fcThreshold, pThreshold});
    trimmedTv.addViewer(volcano);
    applyVolcanoFormulaLines(volcano, fcThreshold, pThreshold);
    await delay(100);

    const viewInfo = trimmedTv.getInfo();
    project.addChild(viewInfo);

    // ─── Step 6.6 — re-create the enrichment Up/Down 2×2 on the published
    // enrichment view. The charts bind to frozenEnrich itself (dockEnrichmentCharts
    // adds the derived negLog10FDR / ~enrichChartTop columns onto it and splits
    // Up/Down via per-viewer formula filters), so there are no chart-backing
    // subset tables to bundle — the enrichment frame is fully self-contained and
    // the charts survive the project round-trip. ─────────────────────────────
    let enrichViewInfo: any = null;
    if (frozenEnrich) {
      pi.description = 'Re-rendering enrichment charts...';
      const enrichTv = grok.shell.addTableView(frozenEnrich);
      await delay(100);
      dockEnrichmentCharts(enrichTv, frozenEnrich);
      await delay(100);
      enrichViewInfo = enrichTv.getInfo();
      project.addChild(enrichViewInfo);
      // Land the user (and downstream readers of grok.shell.tv) back on the
      // protein/volcano view — the primary deliverable — not the enrichment tab.
      grok.shell.v = trimmedTv;
    }

    await grok.dapi.tables.uploadDataFrame(frozen);
    await grok.dapi.tables.save(frozen.getTableInfo());
    if (frozenEnrich) {
      // Upload AFTER dockEnrichmentCharts so the derived chart columns persist.
      await grok.dapi.tables.uploadDataFrame(frozenEnrich);
      await grok.dapi.tables.save(frozenEnrich.getTableInfo());
    }
    await grok.dapi.views.save(viewInfo);
    if (enrichViewInfo) await grok.dapi.views.save(enrichViewInfo);
    await grok.dapi.projects.save(project);
    // meta.publishId stays at the originally-generated UUID — it identifies the
    // *publish event* in the DataFrame's tag/column, distinct from project.id
    // (the Datagrok entity id). Step 9 supersede pointers use project.id directly.

    try { await childClient.addEntity(project.id); } catch (mvErr: any) {
      // Project move to child Space failed — surface but do not abort, since
      // the verify-and-rollback gate below grants directly on the Project
      // (not via Space inheritance per spike A2).
      grok.shell.warning(
        `Could not move Project into child Space "${childName}": ${mvErr?.message ?? mvErr}. ` +
        `Continuing — reviewer view grant is applied directly to the Project.`);
    }

    // ─── Step 7a — grant View directly on the Project (spike A2 — inheritance not visible) ──
    pi.description = 'Granting reviewer view-only access...';
    try {
      await grok.dapi.permissions.grant(project, opts.reviewerGroup, false);
    } catch (grantErr: any) {
      try { await grok.dapi.projects.delete(project); } catch (rb: any) {
        grok.shell.error(`Manual cleanup required: project id ${project.id} could not be deleted: ${rb?.message ?? rb}`);
      }
      throw new Error(`Reviewer group view grant failed: ${grantErr?.message ?? grantErr}`);
    }

    // ─── Step 7b — verify-and-rollback gate (T-15-01, NON-NEGOTIABLE) ───────────
    // Spike A2: `permissions.get(project)` does NOT surface Space-inherited grants.
    // So we check three rings: project, child Space, umbrella Space. Any Edit/Share/
    // Delete on the reviewer group at any level is a leak that fails the gate.
    pi.description = 'Verifying view-only access...';
    const reviewerId = (opts.reviewerGroup as any)?.id;
    const matchesGroup = (g: any): boolean => {
      if (g == null) return false;
      if (typeof g === 'string') return g === reviewerId;
      if (g?.id != null) return g.id === reviewerId;
      return false;
    };

    const projPerm: any = await grok.dapi.permissions.get(project);
    let childSpacePerm: any = null;
    let umbrellaPerm: any = null;
    try { if (childSpace) childSpacePerm = await grok.dapi.permissions.get(childSpace); } catch { /* swallow */ }
    try { if (umbrella) umbrellaPerm = await grok.dapi.permissions.get(umbrella); } catch { /* swallow */ }

    const inAny = (key: string): boolean => {
      for (const perm of [projPerm, childSpacePerm, umbrellaPerm]) {
        const arr = (perm as any)?.[key];
        if (Array.isArray(arr) && arr.some(matchesGroup)) return true;
      }
      return false;
    };

    const inView = inAny('view');
    const inEdit = inAny('edit');
    const inShare = inAny('share');
    const inDelete = inAny('delete');

    if (!inView || inEdit || inShare || inDelete) {
      try { await grok.dapi.projects.delete(project); } catch (rb: any) {
        grok.shell.error(`Manual cleanup required: project id ${project.id} could not be deleted: ${rb?.message ?? rb}`);
      }
      throw new Error(
        'Reviewer group already has elevated access via Space inheritance — ' +
        'publish aborted; ask an admin to scope the umbrella Space\'s permissions');
    }

    // ─── Step 8 — round-trip verification (optional; `verifyPublishedDashboard`,
    // default ON) with W-7 self-heal on formula lines. This is the heavy part:
    // closeAll + reopen the project + assert it survives a reload. Keep it on for
    // client deliverables; turn it off for faster demo shares. The reviewer-access
    // verify-and-rollback gate (Step 7b) is separate and always runs. ────────────
    // Per-share checkbox (opts.verify) wins; fall back to the package setting for
    // non-dialog callers.
    if (opts.verify ?? verifyPublishedDashboard()) {
      pi.description = 'Verifying round-trip survival...';
      // expectedAllowlist captures the actual post-trim post-volcano column state
      // (volcano factory in step 6.5 adds derived columns like '-log10(adj.p-value)'),
      // not just the trim allowlist. The trim allowlist is the floor; the orchestrator
      // may add derived columns; the contract reflects what actually ships.
      // Exclude _meta_* belt-and-braces columns AND ~-prefixed technical columns
      // (e.g. the volcano's '~Volcano label') — assertPublishedShape ignores both
      // on the reopened side, so the expected list must too or the count mismatches.
      const expectedAllowlist = frozen.columns.toList()
        .map((c) => c.name)
        .filter((n) => !n.startsWith('_meta_') && !n.startsWith('~'));
      const contract: PublishedShapeContract = {
        expectedName,
        expectedProjectName: projectName,
        expectedAllowlist,
        expectedMeta: {...meta},
        expectVolcano: true,
        expectEnrichment: !!frozenEnrich,
        expectFormulaLines: true,
      };

      let healedOnce = false;
      while (true) {
        try {
          await assertPublishedShape(project, contract);
          break;
        } catch (assertErr: any) {
          const msg = String(assertErr?.message ?? assertErr ?? '');
          if (!healedOnce && msg.startsWith(FORMULA_LINE_ASSERTION_PREFIX)) {
            healedOnce = true;
            const tv = grok.shell.tv;
            const reopenedVolcano = tv
              ? Array.from(tv.viewers).find((v) => v.type === DG.VIEWER.SCATTER_PLOT)
              : null;
            if (reopenedVolcano && tv) {
              const reDf = tv.dataFrame;
              const fcFromTag = parseFloat(reDf.getTag(PUBLISHED_TAGS.PUBLISHED_FC_THRESHOLD) ?? String(fcThreshold));
              const pFromTag = parseFloat(reDf.getTag(PUBLISHED_TAGS.PUBLISHED_P_THRESHOLD) ?? String(pThreshold));
              applyVolcanoFormulaLines(
                reopenedVolcano,
                Number.isFinite(fcFromTag) ? fcFromTag : fcThreshold,
                Number.isFinite(pFromTag) ? pFromTag : pThreshold,
              );
              try { await grok.dapi.views.save(tv.getInfo()); } catch { /* best-effort */ }
              await delay(100);
              continue;
            }
          }
          try { await grok.dapi.projects.delete(project); } catch (rb: any) {
            grok.shell.error(`Manual cleanup required: project id ${project.id} could not be deleted: ${rb?.message ?? rb}`);
          }
          throw new Error(`Round-trip shape verification failed: ${assertErr?.message ?? assertErr}`);
        }
      }
    }

    // ─── Step 9 — supersede chain (D-04 + W-8 dual-write, non-destructive) ──────
    if (opts.priorVersion != null) {
      pi.description = 'Updating supersede chain...';
      const priorProjectId = (opts.priorVersion as any).id;
      try {
        const priorReloaded = await grok.dapi.projects.find(priorProjectId);
        (priorReloaded as any).options[PUBLISHED_TAGS.SUPERSEDED_BY] = project.id;
        await grok.dapi.projects.save(priorReloaded);
      } catch (priorOptsErr: any) {
        grok.shell.warning(
          `Supersede pointer on prior version's options failed: ${priorOptsErr?.message ?? priorOptsErr} ` +
          `(DataFrame-tag dual-write path below provides fallback).`);
      }

      try {
        grok.shell.closeAll();
        await delay(100);
        const priorReopened = await grok.dapi.projects.find(priorProjectId);
        await priorReopened.open();
        await awaitCheck(
          () => !!grok.shell.tv && !!grok.shell.tv.dataFrame,
          'supersede stamp: prior project DF did not materialize',
          5000,
        );
        const priorDf = grok.shell.tv!.dataFrame;
        priorDf.setTag(PUBLISHED_TAGS.SUPERSEDED_BY, project.id);
        const supersededByCol = priorDf.col(META_COLUMNS.SUPERSEDED_BY);
        if (supersededByCol) {
          try { supersededByCol.set(0, project.id); } catch { /* best-effort */ }
        }
        await grok.dapi.tables.uploadDataFrame(priorDf);
        await grok.dapi.tables.save(priorDf.getTableInfo());
      } catch (stampErr: any) {
        grok.shell.warning(
          `Supersede tag could not be stamped on prior version DF: ${stampErr?.message ?? stampErr} ` +
          `(Project.options pointer was still written — Plan 06 panel falls back to that path)`);
      }

      try { (project as any).options[PUBLISHED_TAGS.SUPERSEDES] = priorProjectId; } catch { /* swallow */ }
      try { await grok.dapi.projects.save(project); } catch (e) {
        grok.shell.warning(`Could not write supersedes pointer on new project options: ${(e as Error)?.message ?? e}`);
      }
    }

    // ─── Step 10 — restore the user's workspace ─────────────────────────────────
    // Round-trip verification (assertPublishedShape) does closeAll() + reopens the
    // published project to prove it survives a reload, and a self-heal step re-saves
    // a viewer (leaving that reopened project dirty → a SAVE badge). The supersede
    // path likewise closes everything and opens the PRIOR version. Either way the
    // user's original analysis is gone and the workspace is cluttered with
    // verification artifacts. Tear all that down and put them back on their source
    // analysis (table + volcano, plus the enrichment view + charts they had open);
    // the shared copy is one click away in the success dialog's "Open shared analysis".
    pi.description = 'Restoring your analysis...';
    try {
      grok.shell.closeAll();
      await delay(100);
      // Restore the enrichment results view + charts (and its volcano cross-link)
      // FIRST if the user had it open. openEnrichmentVisualization focuses the
      // enrichment tab; adding the protein view next lands the user back on the
      // volcano — the primary deliverable — with enrichment one tab away, exactly
      // as it was pre-publish. Uses enrichSource (the original untrimmed table),
      // not the trimmed publish copy.
      if (enrichSource) {
        try { openEnrichmentVisualization(enrichSource, df); } catch { /* best-effort */ }
      }
      const restoredTv = grok.shell.addTableView(df);
      try {
        restoredTv.addViewer(createVolcanoPlot(df));
      } catch { /* volcano is best-effort — the table view is what matters */ }
    } catch { /* restore is best-effort — never fail a successful publish on cleanup */ }

    pi.description = 'Done.';
    // No toast here — the share-dialog caller emits the single user-facing
    // confirmation, so publishing stays quiet to avoid a duplicate.
    return project;
  } finally {
    pi.close();
  }
}
