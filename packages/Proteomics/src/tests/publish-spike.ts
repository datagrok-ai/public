import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, delay, test} from '@datagrok-libraries/test/src/test';
import {SEMTYPE} from '../utils/proteomics-types';

const SPIKE_TAG_KEYS = [
  'proteomics.source',
  'proteomics.de_method',
  'proteomics.de_complete',
  'proteomics.groups',
  'proteomics.published',
  'proteomics.published_at',
  'proteomics.published_by',
  'proteomics.published_target',
  'proteomics.published_de_method',
  'proteomics.published_fc_threshold',
  'proteomics.published_p_threshold',
  'proteomics.published_version',
  'proteomics.published_id',
  'proteomics.published_includes_enrichment',
] as const;

const SPIKE_SEMTYPE_COLS: Array<[string, string]> = [
  ['Protein ID', SEMTYPE.PROTEIN_ID],
  ['Gene Name', SEMTYPE.GENE_SYMBOL],
  ['log2FC', SEMTYPE.LOG2FC],
  ['adj.p-value', SEMTYPE.P_VALUE],
];

function emit(label: string, value: unknown): void {
  let serialized: string;
  try {
    serialized = typeof value === 'string' ? value : JSON.stringify(value);
  } catch {
    serialized = String(value);
  }
  const line = `[publish-spike] ${label}: ${serialized}`;
  console.log(line);
  grok.shell.info(line);
}

function buildFixtureDataFrame(publishedId: string, nowIso: string): DG.DataFrame {
  const proteinIds = Array.from({length: 10}, (_, i) => `P${String(i + 1).padStart(5, '0')}`);
  const geneNames = Array.from({length: 10}, (_, i) => `GENE${i + 1}`);
  const log2fc = new Float32Array([3.0, 2.5, 2.0, 1.5, 0.5, -0.4, -1.6, -2.2, -2.8, 0.1]);
  const pVals = new Float32Array([0.0001, 0.001, 0.005, 0.02, 0.04, 0.1, 0.001, 0.0005, 0.0001, 0.3]);
  const adjP = new Float32Array([0.001, 0.005, 0.01, 0.04, 0.08, 0.2, 0.005, 0.001, 0.0005, 0.5]);
  const significant = ['up', 'up', 'up', 'up', 'ns', 'ns', 'down', 'down', 'down', 'ns'];
  const direction = ['Enriched in Treatment', 'Enriched in Treatment', 'Enriched in Treatment',
    'Enriched in Treatment', 'Not significant', 'Not significant',
    'Enriched in Control', 'Enriched in Control', 'Enriched in Control', 'Not significant'];

  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Protein ID', proteinIds),
    DG.Column.fromStrings('Gene Name', geneNames),
    DG.Column.fromFloat32Array('log2FC', log2fc),
    DG.Column.fromFloat32Array('p-value', pVals),
    DG.Column.fromFloat32Array('adj.p-value', adjP),
    DG.Column.fromStrings('significant', significant),
    DG.Column.fromStrings('direction', direction),
  ]);

  for (const [colName, sem] of SPIKE_SEMTYPE_COLS)
    df.col(colName)!.semType = sem;

  df.name = 'spike-source';

  const groupsPayload = {
    group1: {name: 'Control', columns: []},
    group2: {name: 'Treatment', columns: []},
  };

  const tagValues: Record<string, string> = {
    'proteomics.source': 'generic',
    'proteomics.de_method': 'limma',
    'proteomics.de_complete': 'true',
    'proteomics.groups': JSON.stringify(groupsPayload),
    'proteomics.published': 'true',
    'proteomics.published_at': nowIso,
    'proteomics.published_by': 'spike-operator',
    'proteomics.published_target': 'SPIKE-TARGET-A',
    'proteomics.published_de_method': 'limma',
    'proteomics.published_fc_threshold': '1.0',
    'proteomics.published_p_threshold': '0.05',
    'proteomics.published_version': '1',
    'proteomics.published_id': publishedId,
    'proteomics.published_includes_enrichment': 'false',
  };

  for (const [k, v] of Object.entries(tagValues))
    df.setTag(k, v);

  return df;
}

function summarizeViewers(): Array<Record<string, unknown>> {
  const out: Array<Record<string, unknown>> = [];
  const tv = grok.shell.tv;
  if (!tv) return out;
  for (const v of tv.viewers) {
    let props: any = null;
    try {
      props = (v as any).getOptions ? (v as any).getOptions()?.look ?? (v as any).getOptions() : null;
    } catch (e) {
      props = `getOptions-threw:${(e as Error)?.message ?? e}`;
    }
    let xCol: unknown;
    let yCol: unknown;
    let colorCol: unknown;
    try { xCol = (v as any).props?.xColumnName; } catch { /* swallow */ }
    try { yCol = (v as any).props?.yColumnName; } catch { /* swallow */ }
    try { colorCol = (v as any).props?.colorColumnName; } catch { /* swallow */ }
    out.push({
      type: v.type,
      xColumnName: xCol,
      yColumnName: yCol,
      colorColumnName: colorCol,
      props,
    });
  }
  return out;
}

async function safeDelete(label: string, fn: () => Promise<void>): Promise<void> {
  try {
    await fn();
    emit(`cleanup:${label}`, 'ok');
  } catch (e) {
    emit(`cleanup:${label}:error`, `${(e as Error)?.message ?? e}`);
  }
}

category('Publishing-Spike', () => {
  test('save+open roundtrip — enumerate survivors', async () => {
    const ts = Date.now();
    const publishedId = `spike-pub-${ts}`;
    const nowIso = new Date(ts).toISOString();
    const projectName = `spike-publish-roundtrip-${ts}`;

    emit('boot', {projectName, publishedId, nowIso});

    // (1) Build the fixture DataFrame and prime tags + semTypes.
    const df = buildFixtureDataFrame(publishedId, nowIso);

    // (2) Open in a TableView and add a scatter plot mimicking the volcano shape.
    const tv = grok.shell.addTableView(df);
    await delay(100);
    let sp: DG.ScatterPlotViewer | null = null;
    try {
      const tvAny = tv as any;
      if (typeof tvAny.scatterPlot === 'function') {
        sp = tvAny.scatterPlot({x: 'log2FC', y: 'adj.p-value'});
      } else {
        sp = DG.Viewer.scatterPlot(df, {x: 'log2FC', y: 'adj.p-value'});
        let already = false;
        try { already = Array.from(tv.viewers).includes(sp as any); } catch { /* ignore */ }
        if (sp && !already) tv.addViewer(sp as DG.Viewer);
      }
    } catch (e) {
      emit('scatter-add-error', `${(e as Error)?.message ?? e}`);
    }
    await delay(100);

    // (3) Create Project, set TWO options entries (A4 probe), attach children.
    const project = DG.Project.create();
    project.name = projectName;
    try {
      (project as any).options['proteomics.superseded_by'] = 'dummy-id-A';
      (project as any).options['custom_option'] = 'hello';
    } catch (e) {
      emit('project-options-set-error', `${(e as Error)?.message ?? e}`);
    }

    const tableInfo = tv.dataFrame.getTableInfo();
    const layoutInfo = tv.getInfo();
    project.addChild(tableInfo);
    project.addChild(layoutInfo);

    let tableInfoId: string | null = null;
    let projectId: string | null = null;
    let savedSuccessfully = false;
    try {
      await grok.dapi.tables.uploadDataFrame(tv.dataFrame);
      await grok.dapi.tables.save(tableInfo);
      tableInfoId = tableInfo.id;
      await grok.dapi.views.save(layoutInfo);
      await grok.dapi.projects.save(project);
      projectId = project.id;
      savedSuccessfully = true;
      emit('save:ok', {projectId, tableInfoId});
    } catch (e) {
      emit('save:error', `${(e as Error)?.message ?? e}`);
    }

    if (!savedSuccessfully || !projectId) {
      emit('halt', 'save sequence failed — skipping reopen and downstream probes');
      try { grok.shell.closeAll(); } catch { /* noop */ }
      return;
    }

    // (4) Capture id, close shell, find + open.
    try { grok.shell.closeAll(); } catch (e) { emit('closeAll-error', `${(e as Error)?.message ?? e}`); }
    await delay(100);

    let reopenedProject: DG.Project | null = null;
    try {
      reopenedProject = await grok.dapi.projects.find(projectId);
      await reopenedProject.open();
    } catch (e) {
      emit('reopen:error', `${(e as Error)?.message ?? e}`);
    }
    if (!reopenedProject) {
      emit('halt', 'reopen failed — skipping enumeration');
      return;
    }

    try {
      await awaitCheck(
        () => !!grok.shell.tv && !!grok.shell.tv.dataFrame,
        'TableView did not materialize after Project.open()',
        5000,
      );
    } catch (e) {
      emit('await-tv-error', `${(e as Error)?.message ?? e}`);
    }

    // (5) Enumerate survivors.
    const reTv = grok.shell.tv;
    const reDf = reTv?.dataFrame ?? null;

    emit('df.name survived', reDf?.name ?? '(no dataFrame)');

    if (reDf) {
      const tagSurvival: Record<string, string | null> = {};
      for (const k of SPIKE_TAG_KEYS)
        tagSurvival[k] = reDf.getTag(k) ?? null;
      emit('proteomics.* tags survived', tagSurvival);

      const semTypeSurvival: Record<string, string | null> = {};
      for (const [colName] of SPIKE_SEMTYPE_COLS) {
        const c = reDf.col(colName);
        semTypeSurvival[colName] = c ? (c.semType ?? null) : '(column missing)';
      }
      emit('Proteomics-* semTypes survived', semTypeSurvival);

      const allCols = reDf.columns.toList().map((c) => ({name: c.name, type: c.type, semType: c.semType ?? null}));
      emit('reopened columns shape', allCols);
    }

    // project.options enumeration (A4)
    try {
      const opts = (reopenedProject as any).options ?? {};
      const o1 = opts['proteomics.superseded_by'];
      const o2 = opts['custom_option'];
      let optsKeys: string[] = [];
      try { optsKeys = Object.keys(opts); } catch { /* swallow */ }
      emit('project.options', {
        'proteomics.superseded_by': o1 ?? null,
        'custom_option': o2 ?? null,
        keys: optsKeys,
        rawJson: (() => { try { return JSON.stringify(opts); } catch { return '(unstringifiable)'; } })(),
      });
    } catch (e) {
      emit('project.options:error', `${(e as Error)?.message ?? e}`);
    }

    // Viewer survival + bindings (A6)
    emit('viewers after reopen', summarizeViewers());

    // permissions.get shape (A1)
    try {
      const perm: any = await grok.dapi.permissions.get(reopenedProject);
      let keys: string[] = [];
      try { keys = Object.keys(perm); } catch { /* swallow */ }
      let permJson = '(unstringifiable)';
      try { permJson = JSON.stringify(perm); } catch { /* swallow */ }
      emit('permissions.get keys', keys);
      emit('permissions.get json', permJson);
    } catch (e) {
      emit('permissions.get:error', `${(e as Error)?.message ?? e}`);
    }

    // (6) Space-inheritance probe (A2).
    let umbrellaSpace: any = null;
    let inheritedChildSpace: any = null;
    try {
      umbrellaSpace = await (grok.dapi as any).spaces.createRootSpace(`Spike-Reviews-${ts}`);
      emit('umbrella-space', {id: umbrellaSpace.id, friendlyName: (umbrellaSpace as any).friendlyName});
      const umbrellaClient = (grok.dapi as any).spaces.id(umbrellaSpace.id);
      inheritedChildSpace = await umbrellaClient.addSubspace(`Spike-Review-Inherit-${ts}`);
      emit('inherited-child-space', {id: inheritedChildSpace.id, friendlyName: (inheritedChildSpace as any).friendlyName});
      const inheritedClient = (grok.dapi as any).spaces.id(inheritedChildSpace.id);
      try {
        await inheritedClient.addEntity(reopenedProject.id);
        emit('inherited-child:addEntity', 'ok');
      } catch (e) {
        emit('inherited-child:addEntity:error', `${(e as Error)?.message ?? e}`);
      }
      try {
        const me = await grok.dapi.users.current();
        const myGroup: any = (me as any).group ?? null;
        if (myGroup) {
          await grok.dapi.permissions.grant(umbrellaSpace, myGroup, false);
          emit('umbrella-space:grant-view', {group: (myGroup as any).friendlyName ?? (myGroup as any).name ?? '(unnamed)'});
        } else {
          emit('umbrella-space:grant-skipped', 'no user.group found');
        }
      } catch (e) {
        emit('umbrella-space:grant:error', `${(e as Error)?.message ?? e}`);
      }
      try {
        const childPerm = await grok.dapi.permissions.get(reopenedProject);
        let cpKeys: string[] = [];
        try { cpKeys = Object.keys(childPerm as any); } catch { /* swallow */ }
        emit('permissions.get(childProject) after umbrella-grant', {
          keys: cpKeys,
          json: (() => { try { return JSON.stringify(childPerm); } catch { return '(unstringifiable)'; } })(),
        });
      } catch (e) {
        emit('childProject:permissions.get:error', `${(e as Error)?.message ?? e}`);
      }
    } catch (e) {
      emit('space-probe:error', `${(e as Error)?.message ?? e}`);
    }

    // (7) Smart-filter syntax probe (A8).
    try {
      const likeList: any = await grok.dapi.projects.filter('name like "spike-publish-roundtrip-%"').list();
      emit('filter "like" result', {count: Array.isArray(likeList) ? likeList.length : '(non-array)',
        firstName: Array.isArray(likeList) && likeList.length > 0 ? (likeList[0] as any).name ?? null : null});
    } catch (e) {
      emit('filter:like:error', `${(e as Error)?.message ?? e}`);
    }
    try {
      const containsList: any = await grok.dapi.projects.filter('name contains "spike-publish-roundtrip"').list();
      emit('filter "contains" result', {count: Array.isArray(containsList) ? containsList.length : '(non-array)',
        firstName: Array.isArray(containsList) && containsList.length > 0 ? (containsList[0] as any).name ?? null : null});
    } catch (e) {
      emit('filter:contains:error', `${(e as Error)?.message ?? e}`);
    }

    // (8) Delete-cascade probe (A3).
    const savedTableInfoId = tableInfoId;
    try {
      await grok.dapi.projects.delete(reopenedProject);
      emit('projects.delete', 'ok');
    } catch (e) {
      emit('projects.delete:error', `${(e as Error)?.message ?? e}`);
    }
    if (savedTableInfoId) {
      try {
        const t = await grok.dapi.tables.find(savedTableInfoId);
        emit('tables.find(tableInfoId) after project.delete', {
          resolved: !!t,
          name: t ? (t as any).name ?? null : null,
        });
      } catch (e) {
        emit('tables.find:after-delete:error', `${(e as Error)?.message ?? e}`);
      }
    }

    // (9) Cleanup — idempotent.
    if (inheritedChildSpace)
      await safeDelete('inherited-child-space', () => (grok.dapi as any).spaces.delete(inheritedChildSpace));
    if (umbrellaSpace)
      await safeDelete('umbrella-space', () => (grok.dapi as any).spaces.delete(umbrellaSpace));
    // Best-effort: delete the project again in case the A3 probe earlier already ran but didn't cascade.
    await safeDelete('project (idempotent)', async () => {
      try {
        const stillThere = await grok.dapi.projects.find(projectId!);
        if (stillThere) await grok.dapi.projects.delete(stillThere);
      } catch { /* not found is fine */ }
    });

    emit('done', 'spike enumeration complete — read output above');
  });
});
