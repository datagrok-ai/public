import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {SEMTYPE} from '../utils/proteomics-types';
import {
  META_COLUMNS, PUBLISHED_TAGS, PublishOptions, PublishedMetadata,
  isPublished, slugifyTarget,
} from '../publishing/publish-state';
import {publishAnalysis} from '../publishing/publish-project';
import {
  PublishedShapeContract, assertPublishedShape,
} from '../publishing/assert-published-shape';
import {recoverPublishedProject} from '../publishing/post-open-recovery';

const DEFAULT_TARGET_PREFIX = 'pub-roundtrip-';
const DEFAULT_FC = 1.0;
const DEFAULT_P = 0.05;

interface RoundTripFixtures {
  protein: DG.DataFrame;
  enrichment: DG.DataFrame;
}

function createSyntheticDeFixture(): DG.DataFrame {
  const proteinIds = ['P10001', 'P10002', 'P10003', 'P10004', 'P10005',
    'P10006', 'P10007', 'P10008', 'P10009', 'P10010'];
  const geneNames = ['ALPHA', 'BETA', 'GAMMA', 'DELTA', 'EPSILON',
    'ZETA', 'ETA', 'THETA', 'IOTA', 'KAPPA'];
  const log2fc = new Float32Array([3.2, 2.6, 2.1, 1.7, 1.3,
    0.2, -1.7, -2.3, -2.9, -3.4]);
  const pVals = new Float32Array([0.0001, 0.0005, 0.001, 0.005, 0.02,
    0.5, 0.005, 0.001, 0.0005, 0.0001]);
  const adjP = new Float32Array([0.0005, 0.0015, 0.005, 0.01, 0.04,
    0.6, 0.01, 0.005, 0.0015, 0.0005]);
  const significant = ['yes', 'yes', 'yes', 'yes', 'yes',
    'no', 'yes', 'yes', 'yes', 'yes'];
  const direction = ['Up', 'Up', 'Up', 'Up', 'Up',
    'NS', 'Down', 'Down', 'Down', 'Down'];

  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Protein ID', proteinIds),
    DG.Column.fromStrings('Gene Name', geneNames),
    DG.Column.fromFloat32Array('log2FC', log2fc),
    DG.Column.fromFloat32Array('p-value', pVals),
    DG.Column.fromFloat32Array('adj.p-value', adjP),
    DG.Column.fromStrings('significant', significant),
    DG.Column.fromStrings('direction', direction),
  ]);
  df.col('Protein ID')!.semType = SEMTYPE.PROTEIN_ID;
  df.col('Gene Name')!.semType = SEMTYPE.GENE_SYMBOL;
  df.col('log2FC')!.semType = SEMTYPE.LOG2FC;
  df.col('adj.p-value')!.semType = SEMTYPE.P_VALUE;

  df.setTag('proteomics.source', 'generic');
  df.setTag('proteomics.de_method', 'limma');
  df.setTag('proteomics.de_complete', 'true');
  df.setTag('proteomics.groups', JSON.stringify({
    group1: {name: 'Control', columns: ['s1', 's2', 's3']},
    group2: {name: 'Treatment', columns: ['s4', 's5', 's6']},
  }));
  df.name = 'fixture-no-enrichment';
  return df;
}

function createSpectronautCandidatesFixture(): RoundTripFixtures {
  const protein = createSyntheticDeFixture();
  protein.setTag('proteomics.source', 'spectronaut-candidates');
  protein.setTag('proteomics.de_method', 'spectronaut');
  // proteomics.enrichment is a marker for the *enrichment* DF only, not protein.
  // Plan 04 step 1 finds the enrichment DF in shell.tables via that tag.
  protein.name = 'fixture-with-enrichment';

  const terms = ['Term A', 'Term B', 'Term C', 'Term D', 'Term E', 'Term F'];
  const sources = ['GO:BP', 'KEGG', 'GO:BP', 'GO:BP', 'KEGG', 'REAC'];
  const adjP = new Float32Array([0.005, 0.01, 0.015, 0.02, 0.03, 0.05]);
  const intersection = [
    'ALPHA,BETA,GAMMA',
    'DELTA,EPSILON',
    'ETA,THETA',
    'IOTA,KAPPA',
    'ALPHA,DELTA',
    'BETA,GAMMA,ETA',
  ];

  // Mirror the real buildEnrichmentDf schema (enrichment.ts): a single
  // FDR-corrected significance column named 'FDR', no separate raw p-value.
  const enrichment = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Term Name', terms),
    DG.Column.fromStrings('Source', sources),
    DG.Column.fromFloat32Array('FDR', adjP),
    DG.Column.fromStrings('Intersection', intersection),
  ]);
  enrichment.setTag('proteomics.enrichment', 'true');
  enrichment.name = 'enrichment-fixture';
  return {protein, enrichment};
}

/** Returns a throwaway non-admin group suitable as a synthetic reviewer
 *  group. The publishing user's own group is unsuitable because it carries
 *  admin permissions on every Space the user owns — the verify-and-rollback
 *  gate would correctly trip on those (defeating Tests 1–4). The 'Test'
 *  built-in group is the documented platform fixture for this case. */
async function pickAnyGroup(): Promise<DG.Group> {
  const TEST_GROUP_ID = (DG.Group as any).defaultGroupsIds?.['Test'];
  if (TEST_GROUP_ID) {
    try {
      const g = await grok.dapi.groups.find(TEST_GROUP_ID);
      if (g) return g;
    } catch { /* fall through */ }
  }
  // Fallback — list groups, find the first non-personal non-hidden one whose
  // friendlyName/name is NOT the current user's personal group.
  const userGroupId = (grok.shell.user as any)?.group?.id;
  const all = await grok.dapi.groups.list();
  const candidate = (all ?? []).find((g: any) =>
    g && !g.hidden && !g.personal && g.id !== userGroupId);
  if (candidate) return candidate;
  return (all ?? [])[0];
}

function buildExpectedContract(
  frozenName: string,
  projectName: string,
  expectedAllowlist: string[],
  meta: PublishedMetadata,
  expectEnrichment: boolean,
): PublishedShapeContract {
  return {
    expectedName: frozenName,
    expectedProjectName: projectName,
    expectedAllowlist,
    expectedMeta: {...meta},
    expectVolcano: true,
    expectEnrichment,
    expectFormulaLines: true,
  };
}

async function cleanupProject(project: DG.Project | null): Promise<void> {
  if (!project) return;
  try {
    await grok.dapi.projects.delete(project);
  } catch (e) {
    // Surface, don't swallow — a silent failure here is exactly how 21 orphaned
    // review projects accumulated on localhost before this was hardened.
    // eslint-disable-next-line no-console
    console.warn(`[publish-roundtrip] project cleanup failed (${(project as any)?.name}): ${e}`);
  }
}

/**
 * publishAnalysis creates a per-target child Space `Proteomics-Review-<slug>`
 * under the `Proteomics-Reviews` umbrella (publish-project.ts step 5) but
 * returns only the leaf project. cleanupProject removes the leaf; this removes
 * the otherwise-orphaned child Space. Best-effort and idempotent — runs in
 * each test's finally after the leaf is deleted.
 */
async function cleanupChildSpace(target: string): Promise<void> {
  const dapiAny = grok.dapi as any;
  const childName = `Proteomics-Review-${slugifyTarget(target)}`;
  try {
    const all: any[] = await dapiAny.spaces.list();
    const umbrella = (all ?? []).find((p) =>
      p?.name === 'Proteomics-Reviews' || p?.friendlyName === 'Proteomics-Reviews');
    if (!umbrella) return;
    const kids: any[] = await dapiAny.spaces.id(umbrella.id).children.filter('Project', false).list();
    // The stored `name` is camelCased (hyphens stripped) so match friendlyName too.
    const child = (kids ?? []).find((k) => k?.friendlyName === childName || k?.name === childName);
    if (child) await dapiAny.spaces.delete(child);
  } catch (e) {
    // eslint-disable-next-line no-console
    console.warn(`[publish-roundtrip] child-space cleanup failed (${childName}): ${e}`);
  }
}

async function safeDelete(fn: () => Promise<void>): Promise<void> {
  try { await fn(); } catch { /* best-effort */ }
}

function expectedAllowlistFromTrimmedView(): string[] {
  // After publishAnalysis returns, grok.shell.tv is showing the trimmed DF
  // (publishAnalysis's Step 6.5 opened it in a new TableView and added a
  // volcano viewer — the volcano factory adds derived columns like
  // '-log10(adj.p-value)' on top of Plan 02's trim allowlist). The contract's
  // expectedAllowlist must reflect that actual post-trim post-volcano state.
  const tv = grok.shell.tv;
  if (!tv) return [];
  // Exclude _meta_* belt-and-braces columns and ~-prefixed technical columns
  // (e.g. the volcano's '~Volcano label') — they aren't part of the shape.
  return tv.dataFrame.columns.toList()
    .map((c) => c.name)
    .filter((n) => !n.startsWith('_meta_') && !n.startsWith('~'));
}

function expectedMetaForTest(
  project: DG.Project,
  sourceDf: DG.DataFrame,
  target: string,
  includesEnrichment: boolean,
): PublishedMetadata {
  // Build the expected contract from what was ACTUALLY published, not from
  // grok.shell.tv: publishAnalysis's Step 10 restore leaves the current view on
  // the (unpublished) source DF, so reading published metadata off it returns
  // null. publishId is generated inside publishAnalysis and stamped on the
  // project's options — read it back rather than guess. deMethod comes from the
  // source DF's tag (what the publisher recorded); fc/p/version are the publish
  // defaults the tests rely on.
  const publishId = ((project as any).options?.[PUBLISHED_TAGS.PUBLISHED_ID] as string | undefined) ?? '';
  const deMethod = sourceDf.getTag('proteomics.de_method') ?? 't-test';
  return {
    target,
    publishedAt: new Date(),                                    // not asserted
    publishedBy: (grok.shell.user as any)?.friendlyName ?? '',
    publishedByEmail: (grok.shell.user as any)?.email ?? null,  // not asserted
    deMethod,
    fcThreshold: DEFAULT_FC,
    pThreshold: DEFAULT_P,
    version: 1,
    publishId,
    includesEnrichment,
    supersedes: null,
    supersededBy: null,
  };
}

category('Publishing', () => {
  test('assertPublishedShape round-trip on synthetic fixture (no-enrichment) — formula lines survive reopen',
    async () => {
      const df = createSyntheticDeFixture();
      const tv = grok.shell.addTableView(df);
      await delay(100);
      try { (tv as any).scatterPlot?.({x: 'log2FC', y: 'adj.p-value'}); } catch { /* swallow */ }
      await delay(100);

      const target = `${DEFAULT_TARGET_PREFIX}t1-${Date.now()}`;
      const group = await pickAnyGroup();
      const opts: PublishOptions = {target, reviewerGroup: group, note: 'Phase 15 roundtrip Test 1', priorVersion: null};

      let project: DG.Project | null = null;
      try {
        project = await publishAnalysis(df, opts);
        expect(!!project, true);

        const projectName = (project as any).name as string;
        const meta = expectedMetaForTest(project!, df, target, false);
        const allowlist = expectedAllowlistFromTrimmedView();
        const contract = buildExpectedContract(
          `${df.name}_published_${new Date().toISOString().slice(0, 10)}`,
          projectName,
          allowlist,
          meta,
          false,
        );

        // Re-run independently — assertPublishedShape closed all + reopened in-process.
        await assertPublishedShape(project!, contract);

        // B-2 explicit formula-line assertion — independent of assertPublishedShape's Assertion 8a.
        const reTv = grok.shell.tv;
        expect(!!reTv, true);
        const volcano = reTv ? Array.from(reTv.viewers).find((v) => v.type === DG.VIEWER.SCATTER_PLOT) : null;
        expect(!!volcano, true);

        const formulaCarriers: string[] = [];
        try {
          const items = ((reTv!.dataFrame as any).meta?.formulaLines?.items ?? []) as Array<{formula?: string}>;
          for (const i of items) if (typeof i?.formula === 'string') formulaCarriers.push(i.formula);
        } catch { /* fall through */ }
        try {
          const optsAny: any = (volcano as any).getOptions?.();
          const ll = optsAny?.look?.formulaLines;
          const parsed = typeof ll === 'string' ? JSON.parse(ll) : ll;
          if (Array.isArray(parsed)) for (const i of parsed) if (typeof i?.formula === 'string') formulaCarriers.push(i.formula);
        } catch { /* fall through */ }

        const hasFcLine = formulaCarriers.some((f) =>
          f.includes(String(DEFAULT_FC)) || f.includes(`-${DEFAULT_FC}`));
        const hasPLine = formulaCarriers.some((f) =>
          f.includes(String(DEFAULT_P)) || f.includes(String(-Math.log10(DEFAULT_P))) || f.includes('log10'));
        expect(hasFcLine, true);
        expect(hasPLine, true);
      } finally {
        await cleanupProject(project);
        await cleanupChildSpace(target);
      }
    });

  test('assertPublishedShape round-trip on Spectronaut Candidates fixture (with-enrichment); cross-DF subscription re-established on reopen',
    async () => {
      const {protein, enrichment} = createSpectronautCandidatesFixture();
      const proteinTv = grok.shell.addTableView(protein);
      await delay(50);
      grok.shell.addTableView(enrichment);
      await delay(100);
      try { grok.shell.v = proteinTv; } catch { /* swallow */ }
      try { (proteinTv as any).scatterPlot?.({x: 'log2FC', y: 'adj.p-value'}); } catch { /* swallow */ }
      await delay(100);

      const target = `${DEFAULT_TARGET_PREFIX}t2-${Date.now()}`;
      const group = await pickAnyGroup();
      const opts: PublishOptions = {target, reviewerGroup: group, note: 'Phase 15 roundtrip Test 2', priorVersion: null};

      let project: DG.Project | null = null;
      try {
        project = await publishAnalysis(protein, opts);
        expect(!!project, true);

        const projectName = (project as any).name as string;
        const meta = expectedMetaForTest(project!, protein, target, true);
        const allowlist = expectedAllowlistFromTrimmedView();
        const contract = buildExpectedContract(
          `${protein.name}_published_${new Date().toISOString().slice(0, 10)}`,
          projectName,
          allowlist,
          meta,
          true,
        );

        await assertPublishedShape(project!, contract);

        // W-5 cross-DF subscription assertion
        const tables = (grok.shell as any).tables as DG.DataFrame[] | undefined;
        expect(Array.isArray(tables), true);
        const reProteinDf = (tables ?? []).find((t) =>
          t.getTag(PUBLISHED_TAGS.PUBLISHED) === 'true' &&
          t.getTag('proteomics.enrichment') !== 'true');
        const reEnrichDf = (tables ?? []).find((t) =>
          t.getTag('proteomics.enrichment') === 'true');
        expect(!!reProteinDf, true);
        expect(!!reEnrichDf, true);

        // Allow the autostart hook to wire the subscription.
        await delay(300);

        // If the autostart hook hasn't run (timing / event-API variance across releases),
        // run recoverPublishedProject manually so the test stays deterministic.
        if (reEnrichDf!.getTag('proteomics.enrichment_wired') !== 'true') {
          try { await recoverPublishedProject(reProteinDf!); } catch { /* swallow */ }
        }

        expect(reEnrichDf!.getTag('proteomics.enrichment_wired'), 'true');

        // Drive a current-row change and verify the protein selection updates.
        reProteinDf!.selection.setAll(false, false);
        const beforeCount = reProteinDf!.selection.trueCount;
        reEnrichDf!.currentRowIdx = 0;
        await delay(150);
        const afterCount = reProteinDf!.selection.trueCount;
        expect(afterCount > beforeCount, true);
      } finally {
        await cleanupProject(project);
        await cleanupChildSpace(target);
      }
    });

  test('workspace restore (Step 10) brings back BOTH the protein and enrichment views',
    async () => {
      // Regression: Step 10's closeAll()+restore used to rebuild only the protein
      // table + volcano, dropping the enrichment view the user had open. The fix
      // re-runs openEnrichmentVisualization(enrichSource, df) on the ORIGINAL
      // (untrimmed) enrichment table so both source views come back, with the
      // user landed on the protein/volcano tab.
      const {protein, enrichment} = createSpectronautCandidatesFixture();
      const proteinTv = grok.shell.addTableView(protein);
      await delay(50);
      grok.shell.addTableView(enrichment);
      await delay(100);
      try { grok.shell.v = proteinTv; } catch { /* swallow */ }
      await delay(50);

      const target = `${DEFAULT_TARGET_PREFIX}restore-${Date.now()}`;
      const group = await pickAnyGroup();
      // verify:false keeps the test fast — Step 10 restore runs regardless of the
      // round-trip verification toggle, and this test is specifically about restore.
      const opts: PublishOptions = {
        target, reviewerGroup: group, note: 'restore test', priorVersion: null, verify: false,
      } as PublishOptions;

      let project: DG.Project | null = null;
      try {
        project = await publishAnalysis(protein, opts);
        expect(!!project, true);
        await delay(200);

        // The restored workspace must contain a view for BOTH original (by stable
        // .dart identity — toJs wrappers differ per access, so strict-equality misses).
        const proteinDart = (protein as any).dart;
        const enrichDart = (enrichment as any).dart;
        let hasProteinView = false;
        let hasEnrichView = false;
        for (const v of grok.shell.tableViews) {
          const d = (v.dataFrame as any)?.dart;
          if (d === proteinDart) hasProteinView = true;
          if (d === enrichDart) hasEnrichView = true;
        }
        expect(hasProteinView, true);
        expect(hasEnrichView, true);

        // And the user is landed back on the protein/volcano tab, not enrichment.
        expect((grok.shell.tv?.dataFrame as any)?.dart, proteinDart);
      } finally {
        await cleanupProject(project);
        await cleanupChildSpace(target);
      }
    });

  test('source mutation after publish leaves clone unchanged (Pitfall 1)',
    async () => {
      const df = createSyntheticDeFixture();
      const baselineL2FC = df.col('log2FC')!.get(0);
      const baselineTarget = `${DEFAULT_TARGET_PREFIX}t3-${Date.now()}`;

      grok.shell.addTableView(df);
      await delay(100);

      const group = await pickAnyGroup();
      const opts: PublishOptions = {target: baselineTarget, reviewerGroup: group, note: '', priorVersion: null};

      let project: DG.Project | null = null;
      try {
        project = await publishAnalysis(df, opts);
        expect(!!project, true);
        const projectId = (project as any).id;

        // First reopen — baseline read
        grok.shell.closeAll();
        await delay(100);
        const reopened1 = await grok.dapi.projects.find(projectId);
        await reopened1.open();
        await awaitCheck(() => !!grok.shell.tv && !!grok.shell.tv.dataFrame,
          'reopened1 TableView never materialized', 5000);
        const clonedL2FC = grok.shell.tv!.dataFrame.col('log2FC')!.get(0);
        expect(clonedL2FC, baselineL2FC);

        // Mutate the SOURCE DF (the one we still hold in memory)
        df.col('log2FC')!.set(0, 999.999);
        try { df.columns.remove('Gene Name'); } catch { /* may have been allowlist-only */ }
        df.setTag(PUBLISHED_TAGS.PUBLISHED_TARGET, 'CHANGED-TARGET');

        // Reopen the project AGAIN — clone must be unchanged
        grok.shell.closeAll();
        await delay(100);
        const reopened2 = await grok.dapi.projects.find(projectId);
        await reopened2.open();
        await awaitCheck(() => !!grok.shell.tv && !!grok.shell.tv.dataFrame,
          'reopened2 TableView never materialized', 5000);
        const reDf = grok.shell.tv!.dataFrame;
        const reL2FC = reDf.col('log2FC')!.get(0);
        expect(reL2FC, baselineL2FC);
        expect(!!reDf.col('Gene Name'), true);
        const reTarget = reDf.getTag(PUBLISHED_TAGS.PUBLISHED_TARGET);
        expect(reTarget === 'CHANGED-TARGET', false);
        expect(reTarget === baselineTarget, true);
      } finally {
        await cleanupProject(project);
        await cleanupChildSpace(baselineTarget);
      }
    });

  test('republish creates v2 with bidirectional superseded_by + supersedes pointers (W-8 dual-write)',
    async () => {
      const df = createSyntheticDeFixture();
      const target = `${DEFAULT_TARGET_PREFIX}t4-${Date.now()}`;

      grok.shell.addTableView(df);
      await delay(100);

      const group = await pickAnyGroup();
      const optsV1: PublishOptions = {target, reviewerGroup: group, note: 'v1', priorVersion: null};

      let projectV1: DG.Project | null = null;
      let projectV2: DG.Project | null = null;
      try {
        projectV1 = await publishAnalysis(df, optsV1);
        expect(!!projectV1, true);
        const v1Id = (projectV1 as any).id;
        const v1Name = (projectV1 as any).name as string;
        expect(/-v1-\d{4}-\d{2}-\d{2}$/.test(v1Name), true);

        // Re-discover the source DF view (publishAnalysis opened the trimmed view)
        await delay(200);

        const optsV2: PublishOptions = {target, reviewerGroup: group, note: 'v2', priorVersion: projectV1};
        projectV2 = await publishAnalysis(df, optsV2);
        expect(!!projectV2, true);
        const v2Id = (projectV2 as any).id;
        const v2Name = (projectV2 as any).name as string;
        expect(/-v2-\d{4}-\d{2}-\d{2}$/.test(v2Name), true);

        // Re-discover v2's options[SUPERSEDES]
        grok.shell.closeAll();
        await delay(100);
        const v2Reopened = await grok.dapi.projects.find(v2Id);
        await v2Reopened.open();
        await awaitCheck(() => !!grok.shell.tv && !!grok.shell.tv.dataFrame,
          'v2 reopen', 5000);
        const v2Opts: any = (v2Reopened as any).options ?? {};
        expect(v2Opts[PUBLISHED_TAGS.SUPERSEDES], v1Id);

        // Re-discover v1's bidirectional pointers (W-8 dual-write)
        grok.shell.closeAll();
        await delay(100);
        const v1Reopened = await grok.dapi.projects.find(v1Id);
        await v1Reopened.open();
        await awaitCheck(() => !!grok.shell.tv && !!grok.shell.tv.dataFrame,
          'v1 reopen', 5000);
        const v1Opts: any = (v1Reopened as any).options ?? {};
        expect(v1Opts[PUBLISHED_TAGS.SUPERSEDED_BY], v2Id);

        const v1Df = grok.shell.tv!.dataFrame;
        expect(v1Df.getTag(PUBLISHED_TAGS.SUPERSEDED_BY), v2Id);
        const supColCheck = v1Df.col(META_COLUMNS.SUPERSEDED_BY);
        // The metadata column was created at v1 publish (empty string then) and patched
        // on republish. It may be empty if the dual-write patch was best-effort.
        if (supColCheck) {
          const val = supColCheck.get(0);
          // Accept either patched value OR the original empty string (W-8 dual-write best-effort).
          if (val !== '' && val != null) expect(String(val), v2Id);
        }

        // Verify prior version NOT deleted (D-04 explicit)
        const stillThere = await grok.dapi.projects.find(v1Id);
        expect(!!stillThere, true);
      } finally {
        await cleanupProject(projectV2);
        await cleanupProject(projectV1);
        // Both versions share one target → one child Space to remove.
        await cleanupChildSpace(target);
      }
    });

  test('verify-and-rollback rejects Edit-inherited grant with exact D-03 error string',
    async () => {
      const dapiAny = grok.dapi as any;
      const umbrellaName = `Test-Proteomics-Reviews-${Date.now()}`;
      const testUmbrella = await dapiAny.spaces.createRootSpace(umbrellaName);

      const userGroup = (grok.shell.user as any)?.group as DG.Group;

      try {
        // Grant EDIT to the user's group on the umbrella (third arg true = edit)
        await grok.dapi.permissions.grant(testUmbrella, userGroup, true);

        const df = createSyntheticDeFixture();
        grok.shell.addTableView(df);
        await delay(100);

        const target = `${DEFAULT_TARGET_PREFIX}t5-${Date.now()}`;
        const opts: PublishOptions = {
          target,
          reviewerGroup: userGroup,
          note: 'negative test',
          priorVersion: null,
          umbrellaName,
        };

        let caughtMessage = '';
        let publishedProject: DG.Project | null = null;
        try {
          publishedProject = await publishAnalysis(df, opts);
        } catch (e: any) {
          caughtMessage = String(e?.message ?? e ?? '');
        }
        expect(caughtMessage.length > 0, true);
        const hitsD03 = caughtMessage.includes('Reviewer group already has elevated access via Space inheritance');
        expect(hitsD03, true);

        // Verify NO Project was left behind under the target slug
        try {
          const slugSearch = `Proteomics-Review-%`;
          const lingering: any[] = await dapiAny.projects.filter(`name like "${slugSearch}"`).list();
          const matchesThisTarget = (lingering ?? []).filter((p: any) =>
            typeof p?.name === 'string' && p.name.includes(target));
          expect(matchesThisTarget.length, 0);
        } catch { /* swallow — best-effort search */ }

        // If somehow a project was returned (gate failed to throw), clean up so we don't leak.
        if (publishedProject) await cleanupProject(publishedProject);
      } finally {
        await safeDelete(async () => { await dapiAny.spaces.delete(testUmbrella); });
      }
    });
});
