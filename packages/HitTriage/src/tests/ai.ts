/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, expectFloat, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {_package} from '../package';
import * as api from '../package-api';
import {HitDesignApp} from '../app/hit-design-app';
import {HitDesignTemplate} from '../app/types';
import {allUsers, designCampaign, designCsv, expectSame, noComputeTemplate, Sandbox, triageTemplate} from './fixtures';

const infoNames = ['getHitAppInfo', 'listHitCampaigns', 'getHitCampaignDetails', 'openHitCampaign', 'listHitTemplates',
  'getHitTemplateDetails', 'createHitCampaign', 'deleteHitCampaign'];
const sharedNames = ['getHitCampaignState', 'saveHitCampaign', 'setHitCampaignStatus', 'submitHitCampaign',
  'searchHitTableFunctions', 'callHitTableFunction'];
const designNames = [...sharedNames, 'setHitMolecule', 'setHitStage', 'setHitStages', 'recalculateHitRow'];

const names = (v: DG.ViewBase) => v.getFunctions().map((f) => f.name).sort();

/** Calls a view function as the assistant does: found in the current view's getFunctions(), that view injected. */
async function ai(name: string, params: {[key: string]: unknown} = {}): Promise<any> {
  const v = grok.shell.v;
  const f = v.getFunctions().find((x) => x.name === name);
  if (!f)
    throw new Error(`the view "${v.name}" does not report ${name}; it reports ${names(v).join(', ')}`);
  return await f.apply({view: v, ...params});
}

function expectFailure(res: any, contains: string, what: string): void {
  expect(res?.success === false && String(res?.error).includes(contains), true, `${what}: ${JSON.stringify(res)}`);
}

let sb: Sandbox;
let app: HitDesignApp;
const template: HitDesignTemplate = {...noComputeTemplate('ZZB3 AI Template', 'ZZB3'), compute: {
  descriptors: {enabled: false, args: []}, scripts: [], queries: [],
  functions: [{name: 'addChemPropertiesColumns', package: 'Chem', args: {MW: true, HBA: false, HBD: false, logP: false,
    logS: false, PSA: false, rotatableBonds: false, stereoCenters: false, moleculeCharge: false, molecularFormula: false}}],
}};

const design = () => grok.shell.v as DG.TableView;
function closeDialogs(): void {
  for (const d of DG.Dialog.getOpenDialogs())
    d.close();
}

// One landing and one campaign session for the whole category, as the assistant meets them: the tests run in order.
category('AI', () => {
  before(async () => {
    sb = await Sandbox.create('ai');
    await sb.writeTemplate(template);
    await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', template), designCsv);
    await sb.writeCampaign(designCampaign(sb, 'ZZB3-2', template, {campaignFields: {Chemist: 'Bo'}, status: 'Approved'}), designCsv);
    app = await sb.hitDesign();
  });

  after(async () => {
    await sb.dispose();
  });

  test('the landing reports the info functions to the assistant', async () => {
    expect(grok.shell.v === app.multiView, true, 'the landing is the current view');
    expectSame(names(grok.shell.v), [...infoNames].sort(), 'getFunctions() of the landing');
    expect(String((grok.shell.v as unknown as {aiDescription?: string}).aiDescription).includes('getHitAppInfo'), true, 'aiDescription');
  });

  test('getHitAppInfo', async () => {
    expectSame(await ai('getHitAppInfo'), {app: sb.appName, campaigns: 2, templates: [template.name],
      currentCampaign: null, canCreateCampaignsHere: true}, 'getHitAppInfo');
  });

  test('listHitCampaigns', async () => {
    const res = await ai('listHitCampaigns');
    expect(res.total, 2, 'total');
    const one = res.campaigns.find((c: any) => c.code === 'ZZB3-1');
    expectSame(one, {code: 'ZZB3-1', name: 'Friendly ZZB3-1', template: template.name, status: 'In Progress',
      author: DG.User.current().friendlyName, created: '2024/9/4', rows: 2}, 'the brief of ZZB3-1');
  });

  test('getHitCampaignDetails: by code, by friendly name in any case, and an unknown code', async () => {
    const res = await ai('getHitCampaignDetails', {campaignCode: 'ZZB3-1'});
    expect(res.code, 'ZZB3-1', 'code');
    expectSame(res.campaignFields, designCampaign(sb, 'ZZB3-1', template).campaignFields, 'campaignFields');
    expectSame(res.stages, template.stages, 'stages');
    expectSame(res.compute, {descriptors: [], functions: ['Chem:addChemPropertiesColumns'], scripts: [], queries: []}, 'compute');
    expect(res.canEdit, true, 'canEdit');
    expect(res.lastModifiedBy, 'Someone Else', 'lastModifiedBy');
    expect((await ai('getHitCampaignDetails', {campaignCode: 'friendly zzb3-2'})).code, 'ZZB3-2', 'by friendly name');
    expectFailure(await ai('getHitCampaignDetails', {campaignCode: 'ZZB3-99'}), 'listHitCampaigns', 'an unknown code');
  });

  test('listHitTemplates and getHitTemplateDetails', async () => {
    expectSame(await ai('listHitTemplates'), {templates: [template.name]}, 'listHitTemplates');
    const res = await ai('getHitTemplateDetails', {templateName: template.name});
    expectSame(res, {name: template.name, key: 'ZZB3',
      campaignFields: template.campaignFields.map((f) => ({name: f.name, type: f.type, required: !!f.required})),
      stages: template.stages, compute: {descriptors: [], functions: ['Chem:addChemPropertiesColumns'], scripts: [], queries: []},
      submitFunction: 'HitTriage:demoFileSubmit'}, 'getHitTemplateDetails');
    expectFailure(await ai('getHitTemplateDetails', {templateName: 'No such template'}), 'listHitTemplates', 'an unknown template');
  });

  test('createHitCampaign: refuses a missing required field, then creates and opens the campaign', async () => {
    expectFailure(await ai('createHitCampaign', {templateName: template.name, campaignName: 'AI Three', campaignFields: {}}),
      'Missing required campaign fields: Chemist', 'without Chemist');
    const res = await ai('createHitCampaign', {templateName: template.name, campaignName: 'AI Three', campaignFields: {Chemist: 'Eve'}});
    expectSame(res, {success: true, campaign: 'ZZB3-3', note: 'The campaign is created and its design view is open'}, 'createHitCampaign');
    const saved = await sb.readCampaign('ZZB3-3');
    expect(saved.friendlyName, 'AI Three', 'friendlyName');
    expectSame(saved.campaignFields, {Chemist: 'Eve'}, 'campaignFields');
    expect(design() === app.designView && design().name === 'AI Three', true, 'its design view is current');
    await waitUntil(() => DG.Dialog.getOpenDialogs().length > 0, 'the sketcher did not open for the empty first row', 15000);
    closeDialogs();
    grok.shell.v = app.multiView;
  });

  test('openHitCampaign', async () => {
    expectFailure(await ai('openHitCampaign', {campaignCode: 'ZZB3-99'}), 'listHitCampaigns', 'an unknown code');
    expectSame(await ai('openHitCampaign', {campaignCode: 'ZZB3-1'}), {success: true, campaign: 'ZZB3-1',
      note: 'The campaign view is now open'}, 'openHitCampaign');
    expect(app.campaignId, 'ZZB3-1', 'the app\'s campaign');
    expect(design() === app.designView && design().name === 'Friendly ZZB3-1', true, 'its design view is current');
  });

  test('Q4: grok.shell.v of the campaign view reports the campaign functions', async () => {
    expectSame(names(grok.shell.v), [...designNames].sort(), 'getFunctions() of grok.shell.v');
    expect(String((grok.shell.v as unknown as {aiDescription?: string}).aiDescription).includes('getHitCampaignState'), true, 'aiDescription');
  });

  test('getHitCampaignState', async () => {
    const res = await ai('getHitCampaignState');
    expect(res.app, sb.appName, 'app');
    expect(res.code, 'ZZB3-1', 'code');
    expect(res.campaign?.code, 'ZZB3-1', 'campaign brief');
    expect(res.moleculeColumn, 'Molecule', 'moleculeColumn');
    expectSame(res.stages, template.stages, 'stages');
    expect(res.rows, 2, 'rows');
    expect(res.filteredRows, 2, 'filteredRows');
    expect(['Molecule', 'Stage', 'V-iD', 'Note'].every((c) => res.columns.some((x: any) => x.name === c)), true,
      `columns: ${JSON.stringify(res.columns)}`);
    expect(res.canEdit, true, 'canEdit');
    expect(res.submitFunction, 'HitTriage:demoFileSubmit', 'submitFunction');
    expectSame(res.compute.functions, ['Chem:addChemPropertiesColumns'], 'compute');
  });

  test('saveHitCampaign', async () => {
    const version = (await sb.readCampaign('ZZB3-1')).version;
    expectSame(await ai('saveHitCampaign'), {success: true, campaign: 'ZZB3-1'}, 'saveHitCampaign');
    expect((await sb.readCampaign('ZZB3-1')).version, version + 1, 'saved');
  });

  test('setHitCampaignStatus', async () => {
    expectFailure(await ai('setHitCampaignStatus', {status: ' '}), 'status is required', 'an empty status');
    expectSame(await ai('setHitCampaignStatus', {status: 'Approved'}), {success: true, status: 'Approved'}, 'setHitCampaignStatus');
    expect((await sb.readCampaign('ZZB3-1')).status, 'Approved', 'saved status');
  });

  test('setHitMolecule: a new row and an existing row get the molecule, a V-iD and the computed MW', async () => {
    const df = app.dataFrame!;
    const added = await ai('setHitMolecule', {molecule: 'CCN'});
    await sb.savesDone();
    closeDialogs();
    expect(added.success === true && added.row === 2 && /^V\d{6,}$/.test(added.vId), true, `a new row: ${JSON.stringify(added)}`);
    expect(df.get('Molecule', 2), 'CCN', 'the new row\'s molecule');
    expectFloat(df.get('MW', 2), 45.085, 0.01, 'the new row\'s MW');
    const campaigns = await api.queries.getCampaignsByVid(added.vId, sb.appName);
    expect(campaigns.rowCount === 1 && campaigns.get('campaign_id', 0) === 'ZZB3-1', true, 'the V-iD is registered to the campaign');
    const existing = await ai('setHitMolecule', {molecule: 'Oc1ccccc1', row: 0});
    await sb.savesDone();
    expect(existing.success === true && existing.row === 0 && /^V\d{6,}$/.test(existing.vId), true, `row 0: ${JSON.stringify(existing)}`);
    expectFloat(df.get('MW', 0), 94.113, 0.01, 'row 0 MW');
    expect((await grok.dapi.files.readCsv(sb.savePath('ZZB3-1'))).rowCount, 3, 'the rows saved');
    expectFailure(await ai('setHitMolecule', {molecule: 'CCO', row: 99}), 'row must be between 0 and 2', 'a row out of range');
  });

  test('setHitStage', async () => {
    expectSame(await ai('setHitStage', {row: 1, stage: 'Tested'}), {success: true, row: 1, stage: 'Tested'}, 'setHitStage');
    expect((await grok.dapi.files.readCsv(sb.savePath('ZZB3-1'))).get('Stage', 1), 'Tested', 'saved stage');
    expectFailure(await ai('setHitStage', {row: 1, stage: 'Nope'}), 'the stages are: Idea, Made, Tested', 'an unknown stage');
  });

  test('setHitStages', async () => {
    const stages = ['Idea', 'Made', 'Tested', 'Shipped'];
    expectSame(await ai('setHitStages', {stages}), {success: true, stages}, 'setHitStages');
    expectSame((await sb.readCampaign('ZZB3-1')).template.stages, stages, 'saved stages');
    expectFailure(await ai('setHitStages', {stages: []}), 'non-empty list', 'no stages');
  });

  test('recalculateHitRow', async () => {
    app.dataFrame!.set('MW', 1, null);
    expectSame(await ai('recalculateHitRow', {row: 1}), {success: true, row: 1}, 'recalculateHitRow');
    await sb.savesDone();
    expectFloat(app.dataFrame!.get('MW', 1), 78.114, 0.01, 'benzene MW recalculated');
    expectFailure(await ai('recalculateHitRow', {row: 5}), 'row must be between 0 and 2', 'a row out of range');
  });

  test('searchHitTableFunctions and callHitTableFunction reach the table view\'s own commands', async () => {
    const found = await ai('searchHitTableFunctions', {query: 'select all', limit: 50});
    expect(found.functions.some((f: any) => f.name === 'Select all'), true, `found: ${found.functions.map((f: any) => f.name).join(', ')}`);
    const df = app.dataFrame!;
    df.selection.setAll(false);
    expect((await ai('callHitTableFunction', {name: 'Select all', parameters: {}})).success, true, 'Select all ran');
    expect(df.selection.trueCount, df.rowCount, 'every row selected');
    df.selection.setAll(false);
    expectFailure(await ai('callHitTableFunction', {name: 'No such command', parameters: {}}), 'searchHitTableFunctions', 'an unknown command');
  });

  test('submitHitCampaign', async () => {
    expectSame(await ai('submitHitCampaign', {status: 'Submitted by the assistant'}),
      {success: true, submittedTo: 'HitTriage:demoFileSubmit'}, 'submitHitCampaign');
    expect((await sb.readCampaign('ZZB3-1')).status, 'Submitted by the assistant', 'saved status');
  });

  test('deleteHitCampaign', async () => {
    grok.shell.v = app.multiView;
    expectFailure(await ai('deleteHitCampaign', {campaignCode: 'ZZB3-99'}), 'listHitCampaigns', 'an unknown code');
    expectSame(await ai('deleteHitCampaign', {campaignCode: 'ZZB3-2'}), {success: true, deleted: 'ZZB3-2'}, 'deleteHitCampaign');
    expect(await _package.files.exists(`${sb.appName}/campaigns/ZZB3-2`), false, 'the folder is gone');
    expectSame((await ai('listHitCampaigns')).campaigns.map((c: any) => c.code).sort(), ['ZZB3-1', 'ZZB3-3'], 'the list after the delete');
  });

  test('Hit Triage: the pick view reports the shared campaign functions', async () => {
    const ht = await Sandbox.create('ai');
    try {
      const t = triageTemplate('ZZB3 AI Triage', 'ZZB3');
      await ht.writeTemplate(t);
      await ht.writeCampaign({name: 'ZZB3-1', templateName: t.name, filters: [],
        ingest: {type: 'File', query: ht.ingestPath('ZZB3-1'), molColName: 'smiles'}, status: 'In Progress',
        createDate: '2024/10/17', campaignFields: {Chemist: 'Ann'}, columnSemTypes: {smiles: 'Molecule'},
        columnTypes: {smiles: 'string'}, rowCount: 1, filteredRowCount: 1, authorUserId: DG.User.current().id,
        permissions: {view: [allUsers], edit: [allUsers]}} as {name: string}, 'smiles\nCCO\n');
      const triage = await ht.hitTriage();
      expectSame(names(grok.shell.v), [...infoNames].sort(), 'getFunctions() of the Hit Triage landing');
      const pick = await ht.openTriage(triage, 'ZZB3-1');
      grok.shell.v = pick;
      expectSame(names(grok.shell.v), [...sharedNames].sort(), 'getFunctions() of the pick view');
      const state = await ai('getHitCampaignState');
      expect(state.code === 'ZZB3-1' && state.moleculeColumn === 'smiles' && state.canEdit === true, true, JSON.stringify(state));
      expectSame(await ai('saveHitCampaign'), {success: true, campaign: 'ZZB3-1'}, 'saveHitCampaign');
      expect((await ht.readCampaign('ZZB3-1')).lastModifiedUserName, DG.User.current().friendlyName, 'saved');
    } finally {
      await ht.dispose();
      grok.shell.v = app.multiView;
    }
  }, {timeout: 120000});
}, {clear: false});
