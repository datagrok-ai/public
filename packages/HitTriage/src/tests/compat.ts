/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {_package} from '../package';
import {loadCampaigns, toFormatedDateString} from '../app/utils';
import {CampaignJsonName, CampaignTableName, HitSelectionColName} from '../app/consts';
import {HitDesignApp} from '../app/hit-design-app';
import {HitTriageApp} from '../app/hit-triage-app';
import {allUsers, designCampaign, designColumns, designCsv, expectSame, hasBalloon, noComputeTemplate, oldDate,
  Sandbox, triageTemplate} from './fixtures';

/** The keys a Hit Design / PeptiHit save writes. */
const designKeys = ['name', 'friendlyName', 'templateName', 'status', 'createDate', 'campaignFields', 'columnSemTypes',
  'rowCount', 'filteredRowCount', 'savePath', 'columnTypes', 'columnTags', 'version', 'authorUserId',
  'authorUserFriendlyName', 'lastModifiedUserName', 'permissions', 'layout', 'template'];

const sketch = JSON.stringify({'#type': 'SketchState', 'elementStates': [{left: 10, top: 20, width: 100, height: 20,
  type: 'field', viewerSettings: {table: 'Table', column: 'Stage'}}]});
const mergeConfig = {mode: 'smiles', molColName: 'Molecule', addNewRows: true, runComputeOnNewRows: false, clashStrategy: 'skip'};

function expectKeys(json: {[key: string]: unknown}, keys: string[]): void {
  const missing = keys.filter((k) => !(k in json));
  expect(missing.join(', '), '', 'keys missing after the save');
}

async function listed(sb: Sandbox, id: string) {
  const all = await loadCampaigns(sb.appName, []);
  if (!all[id])
    throw new Error(`${id} is not listed; the list has ${Object.keys(all).join(', ')}`);
  return all[id];
}

category('Compat: campaign.json', () => {
  test('Hit Design: the current shape round-trips', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = noComputeTemplate('ZZB3 Design', 'ZZB3');
      await sb.writeTemplate(t);
      const before = designCampaign(sb, 'ZZB3-1', t, {
        columnTags: {Note: {'.choices': '["first","second"]'}}, tilesViewerFormSketch: sketch, mergeConfig});
      await sb.writeCampaign(before, designCsv);
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      await app.saveCampaign(false);

      const after = await sb.readCampaign('ZZB3-1');
      expectKeys(after, [...designKeys, 'tilesViewerFormSketch', 'mergeConfig']);
      for (const k of ['name', 'friendlyName', 'templateName', 'status', 'createDate', 'savePath', 'tilesViewerFormSketch',
        'authorUserId', 'authorUserFriendlyName'])
        expect(after[k], (before as {[key: string]: unknown})[k], k);
      expectSame(after.campaignFields, before.campaignFields, 'campaignFields');
      expectSame(after.permissions, before.permissions, 'permissions');
      expectSame(after.mergeConfig, before.mergeConfig, 'mergeConfig');
      expectSame(after.template, t, 'template');
      expect(after.version, 4, 'version');
      expect(after.lastModifiedUserName, DG.User.current().friendlyName, 'lastModifiedUserName');
      expect(after.rowCount, 2, 'rowCount');
      expect(after.columnSemTypes.Molecule, 'Molecule', 'columnSemTypes.Molecule');
      expect(after.columnTags?.Note?.['.choices'], '["first","second"]', 'columnTags.Note');
      expect(typeof after.layout === 'string' && DG.ViewLayout.fromViewState(after.layout) != null, true, 'layout is a view state');
      const df = await grok.dapi.files.readCsv(before.savePath);
      expect(df.rowCount, 2, 'rows at savePath');
      expect(['Molecule', 'Stage', 'V-iD', 'Note'].every((c) => df.columns.contains(c)), true, 'columns at savePath');
      expect((await listed(sb, 'ZZB3-1')).friendlyName, before.friendlyName, 'listed after the save');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: a legacy campaign without author, friendly name, column tags or template (DMT1-1)', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = noComputeTemplate('ZZB3 Legacy', 'ZZB3');
      await sb.writeTemplate(t);
      const before = {name: 'ZZB3-1', templateName: t.name, status: 'In Progress', createDate: '2023/11/11',
        campaignFields: {Chemist: 'mnmnmn', Deadline: oldDate}, ...designColumns, rowCount: 2, filteredRowCount: 2,
        savePath: sb.savePath('ZZB3-1'), version: 1, permissions: {view: [allUsers], edit: [allUsers]}};
      await sb.writeCampaign(before, designCsv);
      expect((await listed(sb, 'ZZB3-1')).name, 'ZZB3-1', 'listed before any save');
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      expectSame(app.template, t, 'the template read from its file');
      await app.saveCampaign(false);

      const after = await sb.readCampaign('ZZB3-1');
      expectKeys(after, designKeys);
      expect(after.friendlyName, 'ZZB3-1', 'friendlyName falls back to the code');
      expect(after.createDate, '2023/11/11', 'createDate');
      expectSame(after.campaignFields, before.campaignFields, 'campaignFields, the old Date shape included');
      expectSame(after.template, t, 'the template file becomes the snapshot');
      expect(after.version, 2, 'version');
      expect(after.authorUserId, DG.User.current().id, 'the first save records its user as the author');
      expect(after.authorUserFriendlyName, DG.User.current().friendlyName, 'authorUserFriendlyName');
      expectSame(Object.keys(after.columnTags).sort(), ['Molecule', 'Note', 'Stage', 'V-iD'], 'columnTags of every column');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: missing optionals get their defaults', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = noComputeTemplate('ZZB3 Minimal', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign({name: 'ZZB3-1', templateName: t.name, campaignFields: {Chemist: 'x'}} as {name: string}, designCsv);
      await listed(sb, 'ZZB3-1');
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      expect(app.dataFrame!.rowCount, 2, 'rows read from the default savePath');
      await app.saveCampaign(false);

      const after = await sb.readCampaign('ZZB3-1');
      expectKeys(after, designKeys);
      expect(after.status, 'In Progress', 'status');
      expect(after.createDate, toFormatedDateString(new Date()), 'createDate is today in YYYY/M/D');
      expect(after.version, 1, 'version');
      expect(after.savePath, sb.savePath('ZZB3-1'), 'savePath is the default, in the dot form');
      expect(after.friendlyName, 'ZZB3-1', 'friendlyName');
      expect(after.permissions?.view?.length > 0 && after.permissions?.edit?.length > 0, true, 'permissions get the defaults');
      expect(after.authorUserId, DG.User.current().id, 'authorUserId');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: a campaign whose template file is gone opens from its snapshot', async () => {
    const sb = await Sandbox.create('compat');
    try {
      await sb.writeTemplate(noComputeTemplate('ZZB3 Present', 'ZZB3P'));
      const gone = noComputeTemplate('ZZB3 Gone', 'ZZB3');
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', gone), designCsv);
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      await app.saveCampaign(false);
      const after = await sb.readCampaign('ZZB3-1');
      expect(after.templateName, 'ZZB3 Gone', 'templateName');
      expectSame(after.template, gone, 'template snapshot');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a legacy campaign whose template is gone says so and logs no error; an unreadable campaign.json is logged', async () => {
    const sb = await Sandbox.create('compat');
    const logged: unknown[] = [];
    const logError = _package.logger.error;
    _package.logger.error = (e: unknown) => logged.push(e);
    try {
      const {template: _, ...legacy} = designCampaign(sb, 'ZZB3-1', noComputeTemplate('ZZB3 Gone', 'ZZB3'));
      await sb.writeCampaign(legacy as {name: string}, designCsv);
      await sb.writeCampaign({name: 'ZZB3-2'}, designCsv);
      await _package.files.writeAsText(sb.path('ZZB3-2'), '{not json');
      const app = await sb.hitDesign();
      await app.infoView.setCampaign('ZZB3-1');
      expect(hasBalloon('Campaign ZZB3-1 could not be opened: its template is missing'), true, 'the user is told');
      expect(logged.length, 0, `nothing logged as an error for an expected state; logged ${logged.map(String).join('; ')}`);
      await app.infoView.setCampaign('ZZB3-2');
      expect(hasBalloon('Campaign ZZB3-2 could not be opened'), true, 'the user is told of the failure too');
      expect(logged.length, 1, 'a real failure (campaign.json does not parse) is logged as an error');
    } finally {
      _package.logger.error = logError;
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('the landing and the app\'s view carry the app\'s README as their help, as the platform shows a view\'s help', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const design = await sb.hitDesign();
      const triage = await sb.hitTriage();
      for (const [app, readme] of [[design, 'README_HD.md'], [triage, 'README_HT.md']] as const) {
        expect(app.infoView.helpUrl, _package.webRoot + readme, `${readme}: the landing's help`);
        expect(app.multiView.helpUrl, _package.webRoot + readme, `${readme}: the help of the view the shell shows`);
      }
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: a savePath in another folder is read and written there', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = noComputeTemplate('ZZB3 Elsewhere', 'ZZB3');
      await sb.writeTemplate(t);
      const savePath = sb.savePath('ZZB3-1', 'elsewhere');
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {savePath}), designCsv, savePath);
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      expect(app.dataFrame!.rowCount, 2, 'rows read from savePath');
      app.dataFrame!.set('Note', 0, 'edited');
      await app.saveCampaign(false);
      const after = await sb.readCampaign('ZZB3-1');
      expect(after.savePath, savePath, 'savePath');
      expect((await grok.dapi.files.readCsv(savePath)).get('Note', 0), 'edited', 'the table written at savePath');
      expect(await grok.dapi.files.exists(sb.ingestPath('ZZB3-1')), false, 'no table in the campaign folder');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: a merge file kept with Auto-merge is merged on open without a save; one kept without it is not', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = noComputeTemplate('ZZB3 Merge', 'ZZB3');
      await sb.writeTemplate(t);
      const filePath = `System:AppData/HitTriage/${sb.appName}/merge-source.csv`;
      await grok.dapi.files.writeAsText(filePath, 'Molecule,Note\nOCC,merged\n');
      const merge = {mode: 'smiles', filePath, molColName: 'Molecule', addNewRows: false, runComputeOnNewRows: false,
        clashStrategy: 'overwrite'};
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {mergeConfig: {...merge, autoMergeOnOpen: true}}), designCsv);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-2', t, {mergeConfig: {...merge, autoMergeOnOpen: false}}), designCsv);
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      expect(app.dataFrame!.get('Note', 0), 'merged', 'the file\'s value in the matched row (OCC is CCO)');
      expect(app.dataFrame!.get('Note', 1), 'second', 'the row the file does not hold');
      expect(await sb.lastModified('ZZB3-1'), null, 'no save logged by the open');
      await sb.openDesign(app, 'ZZB3-2');
      expect(app.dataFrame!.get('Note', 0), 'first', 'without Auto-merge the open reads the campaign as saved');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: column tags are restored on open and saved back', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = noComputeTemplate('ZZB3 Tags', 'ZZB3');
      await sb.writeTemplate(t);
      const order = ['Tested', 'Made', 'Idea'];
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {columnTags: {Stage: {'.category-order': JSON.stringify(order)},
        Note: {'.choices': '["first","second"]'}}}), designCsv);
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      expectSame(app.dataFrame!.col('Stage')!.getCategoryOrder(), order, 'Stage category order');
      expect(app.dataFrame!.col('Note')!.getTag('.choices'), '["first","second"]', 'Note .choices');
      await app.saveCampaign(false);
      const after = await sb.readCampaign('ZZB3-1');
      expectSame(JSON.parse(after.columnTags.Stage['.category-order']), order, 'saved category order');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('PeptiHit: the PHD shape round-trips', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = {...noComputeTemplate('ZZB3 Peptides', 'ZZB3'), stages: ['Design', 'Synthesis']};
      await sb.writeTemplate(t);
      const csv = 'Helm,Molecule,Stage,V-iD\nPEPTIDE1{A.G.L.V.K}$$$$V2.0,CC(N)C(=O)O,Design,\nPEPTIDE1{G.A.V}$$$$V2.0,NCC(=O)O,Synthesis,\n';
      const before = {name: 'ZZB3-1', templateName: t.name, status: 'In Progress', createDate: '2024/9/4',
        campaignFields: {Chemist: 'Davit'},
        columnSemTypes: {'Helm': 'Macromolecule', 'Molecule': 'Molecule', 'Stage': null, 'V-iD': null},
        columnTypes: {'Helm': 'string', 'Molecule': 'string', 'Stage': 'string', 'V-iD': 'string'},
        rowCount: 2, filteredRowCount: 2, savePath: sb.savePath('ZZB3-1'), tilesViewerFormSketch: sketch, version: 26,
        permissions: {view: [allUsers], edit: [allUsers]}, template: t};
      await sb.writeCampaign(before, csv);
      const app = await sb.peptiHit();
      await sb.openDesign(app, 'ZZB3-1');
      expect(app.helmColName, 'Helm', 'helm column');
      await app.saveCampaign(false);

      const after = await sb.readCampaign('ZZB3-1');
      expectKeys(after, [...designKeys, 'tilesViewerFormSketch']);
      expect(after.tilesViewerFormSketch, sketch, 'tilesViewerFormSketch');
      expect(after.version, 27, 'version');
      expect(after.columnSemTypes.Helm, 'Macromolecule', 'Helm semType');
      expect(after.columnSemTypes.Molecule, 'Molecule', 'Molecule semType');
      expectSame(after.template, t, 'template');
      expect(after.friendlyName, 'ZZB3-1', 'friendlyName');
      const df = await grok.dapi.files.readCsv(before.savePath);
      expect(df.get('Helm', 0), 'PEPTIDE1{A.G.L.V.K}$$$$V2.0', 'Helm written back');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('Hit Triage: the HTD shape round-trips: empty-object filters, colon-form ingest, old dates', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = triageTemplate('ZZB3 Triage', 'ZZB3');
      await sb.writeTemplate(t);
      const before = {name: 'ZZB3-1', templateName: t.name, filters: {},
        ingest: {type: 'File', query: sb.ingestPath('ZZB3-1'), molColName: 'smiles'}, status: 'In Progress',
        createDate: '2024/10/17', campaignFields: {Chemist: 'Davit', Deadline: oldDate, Scaffold: 'c1cc2ccccc2cc1'},
        columnSemTypes: {smiles: 'Molecule', MolWt: null}, rowCount: 2, filteredRowCount: 2, template: t,
        columnTypes: {smiles: 'string', MolWt: 'double'}, authorUserId: DG.User.current().id,
        permissions: {view: [allUsers], edit: [allUsers]}};
      await sb.writeCampaign(before, 'smiles,MolWt\nCCO,46.07\nc1ccccc1,78.11\n');
      const app = await sb.hitTriage();
      await sb.openTriage(app, 'ZZB3-1');
      await app.saveCampaign(undefined, false);

      const after = await sb.readCampaign('ZZB3-1');
      expect(Array.isArray(after.filters), true, 'empty-object filters are saved as the filter panel list');
      expectSame(after.ingest, {type: 'File', query: sb.ingestPath('ZZB3-1'), molColName: 'smiles'}, 'ingest, colon form');
      expectSame(after.campaignFields, before.campaignFields, 'campaignFields');
      for (const k of ['templateName', 'status', 'createDate', 'authorUserId'])
        expect(after[k], (before as {[key: string]: unknown})[k], k);
      expect('friendlyName' in after, false, 'friendlyName stays absent');
      expectSame(after.permissions, before.permissions, 'permissions');
      expectSame(after.template, t, 'template = the file plus its data source');
      expect(after.lastModifiedUserName, DG.User.current().friendlyName, 'lastModifiedUserName');
      const df = await grok.dapi.files.readCsv(sb.ingestPath('ZZB3-1'));
      expect(df.rowCount, 2, 'rows');
      expect(df.columns.contains(HitSelectionColName), true, `"${HitSelectionColName}" written`);
      await listed(sb, 'ZZB3-1');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});
});

type DemoApp = 'Hit Design' | 'PeptiHit' | 'Hit Triage' | 'PepTriage';

/** Copies a shipped demo campaign (and its template file) into [sb], pointing its table paths at the copy. */
async function copyDemo(sb: Sandbox, app: DemoApp, id: string): Promise<{[key: string]: any}> {
  const json = JSON.parse(await _package.files.readAsText(`${app}/campaigns/${id}/${CampaignJsonName}`));
  const csv = await _package.files.readAsText(`${app}/campaigns/${id}/${CampaignTableName}`);
  if ('savePath' in json)
    json.savePath = sb.savePath(id);
  if (json.ingest)
    json.ingest.query = sb.ingestPath(id);
  const templatePath = `${app}/templates/${json.templateName}.json`;
  await sb.writeTemplate(await _package.files.exists(templatePath) ? JSON.parse(await _package.files.readAsText(templatePath)) : json.template);
  await sb.writeCampaign(json, csv);
  return json;
}

async function openDemo(sb: Sandbox, app: DemoApp, id: string): Promise<HitDesignApp | HitTriageApp> {
  if (app === 'Hit Design' || app === 'PeptiHit') {
    const a = app === 'Hit Design' ? await sb.hitDesign() : await sb.peptiHit();
    await sb.openDesign(a, id);
    await a.saveCampaign(false);
    return a;
  }
  const a = app === 'Hit Triage' ? await sb.hitTriage() : await sb.pepTriage();
  await sb.openTriage(a, id);
  await a.saveCampaign(undefined, false);
  return a;
}

category('Compat: campaign.json', () => {
  test('Hit Design keeps campaign.json keys it does not know', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = noComputeTemplate('ZZB3 Keys', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {zzFuture: {kept: true}}), designCsv);
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      await app.saveCampaign(false);
      const after = await sb.readCampaign('ZZB3-1');
      expect(after.version, 4, 'the save happened');
      expect(after.zzFuture?.kept, true, 'the unknown key zzFuture kept after a save');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Triage keeps campaign.json keys it does not know', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = triageTemplate('ZZB3 Triage Keys', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign({name: 'ZZB3-1', templateName: t.name, filters: [],
        ingest: {type: 'File', query: sb.ingestPath('ZZB3-1'), molColName: 'smiles'}, status: 'In Progress',
        createDate: '2024/10/17', campaignFields: {Chemist: 'x'}, columnSemTypes: {smiles: 'Molecule'},
        columnTypes: {smiles: 'string'}, template: t, zzFuture: {kept: true}} as {name: string}, 'smiles\nCCO\n');
      const app = await sb.hitTriage();
      await sb.openTriage(app, 'ZZB3-1');
      await app.saveCampaign(undefined, false);
      const after = await sb.readCampaign('ZZB3-1');
      expect(after.lastModifiedUserName, DG.User.current().friendlyName, 'the save happened');
      expect(after.zzFuture?.kept, true, 'the unknown key zzFuture kept after a save');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Triage: campaign.json names no hidden column; one an older save named still opens', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = triageTemplate('ZZB3 Triage Hidden', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign({name: 'ZZB3-1', templateName: t.name, filters: [],
        ingest: {type: 'File', query: sb.ingestPath('ZZB3-1'), molColName: 'smiles'}, status: 'In Progress',
        createDate: '2024/10/17', campaignFields: {Chemist: 'x'},
        columnSemTypes: {'smiles': 'Molecule', '~smiles.canonicalSmiles': 'Molecule'},
        columnTypes: {'smiles': 'string', '~smiles.canonicalSmiles': 'string'}, template: t} as {name: string},
      'smiles\nCCO\nc1ccccc1\n');
      const app = await sb.hitTriage();
      await sb.openTriage(app, 'ZZB3-1');
      const df = app.dataFrame!;
      const columns = df.columns.names();
      expectSame(columns.filter((n) => n.startsWith('~')), [], 'hidden columns of the table an older campaign.json opens');
      // what Chem keeps once a filter or a search read the molecules
      const twin = df.columns.addNewString('~smiles.canonicalSmiles');
      twin.init((i) => df.get('smiles', i));
      twin.semType = DG.SEMTYPE.MOLECULE;
      await app.saveCampaign(undefined, false);
      const json = await sb.readCampaign('ZZB3-1');
      expectSame([...Object.keys(json.columnSemTypes), ...Object.keys(json.columnTypes)].filter((n) => n.startsWith('~')),
        [], 'hidden columns campaign.json names');
      const reopened = await sb.hitTriage();
      await sb.openTriage(reopened, 'ZZB3-1');
      expectSame(reopened.dataFrame!.columns.names(), columns, 'the columns the campaign reopens with');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a closed Hit Triage app stops reacting to view changes', async () => {
    const sb = await Sandbox.create('compat');
    try {
      const t = triageTemplate('ZZB3 Leak', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign({name: 'ZZB3-1', templateName: t.name, filters: [],
        ingest: {type: 'File', query: sb.ingestPath('ZZB3-1'), molColName: 'smiles'}, status: 'In Progress',
        createDate: '2024/10/17', campaignFields: {Chemist: 'x'}, columnSemTypes: {smiles: 'Molecule'},
        columnTypes: {smiles: 'string'}} as {name: string}, 'smiles\nCCO\n');
      const app = await sb.hitTriage();
      const pick = await sb.openTriage(app, 'ZZB3-1');
      const pickName = pick.name;
      pick.close();
      app.multiView.close();
      history.replaceState(history.state, '', `${location.pathname}?zzb3=probe`);
      const probe = DG.View.create();
      probe.name = pickName;
      let shown = false;
      const sub = grok.events.onCurrentViewChanged.subscribe(() => shown ||= grok.shell.v?.name === pickName);
      grok.shell.addView(probe);
      try {
        await waitUntil(() => shown, 'the probe view never became current', 10000);
        expect(location.search.includes('campaignId'), false, `the URL left alone once the app is closed: ${location.href}`);
      } finally {
        sub.unsubscribe();
        probe.close();
      }
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});
});

category('Compat: demo campaigns', () => {
  const demos: [DemoApp, string][] = [['Hit Design', 'DMT1-1'], ['Hit Design', 'DMT2-1'], ['PeptiHit', 'PHD-1'],
    ['PeptiHit', 'PHD-2'], ['Hit Triage', 'HTD-1'], ['Hit Triage', 'HTD-2'], ['PepTriage', 'DMT-1']];
  for (const [app, id] of demos) {
    test(`${app} ${id}: the shipped campaign opens and saves with every key kept`, async () => {
      const sb = await Sandbox.create('demo');
      try {
        const before = await copyDemo(sb, app, id);
        await listed(sb, id);
        await openDemo(sb, app, id);
        const after = await sb.readCampaign(id);
        expectKeys(after, Object.keys(before));
        for (const k of ['name', 'templateName', 'status', 'createDate'])
          expect(after[k], before[k], k);
        expect(after.friendlyName ?? after.name, before.friendlyName ?? before.name, 'the shown name');
        expectSame(after.campaignFields, before.campaignFields, 'campaignFields');
        expectSame(after.permissions, before.permissions, 'permissions');
        if (before.ingest) {
          expect(Array.isArray(after.filters), true, 'filters are a list');
          expectSame(after.ingest, before.ingest, 'ingest');
          expect(after.template?.name === before.templateName && after.template?.key === before.template?.key, true, 'the template');
        } else {
          expect(after.version, (before.version ?? 0) + 1, 'version');
          expect(after.savePath, before.savePath, 'savePath');
          expectSame(after.template, before.template, 'template snapshot');
          expect(after.tilesViewerFormSketch ?? null, before.tilesViewerFormSketch ?? null, 'tilesViewerFormSketch');
        }
        const df = await grok.dapi.files.readCsv(before.savePath ?? before.ingest.query);
        expect(df.rowCount, before.rowCount, 'the table rows');
        await listed(sb, id);
      } finally {
        await sb.dispose();
      }
    }, {timeout: 150000});
  }
});
