/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {_package} from '../package';
import {loadTemplate} from '../app/utils';
import {HitDesignApp} from '../app/hit-design-app';
import {HitTriageApp} from '../app/hit-triage-app';
import {hitInfoViewFunctions} from '../app/ai/view-functions';
import {allUsers, expectSame, noComputeTemplate, pepTriageTemplate, Sandbox, startDesignCampaign,
  triageTemplate} from './fixtures';

/** The template file name is the template name, spaces included, and every reader takes it back unchanged. */
async function expectTemplateFile(sb: Sandbox, t: {name: string}): Promise<void> {
  const files = (await _package.files.list(`${sb.appName}/templates`)).map((f) => f.name);
  expect(files.includes(`${t.name}.json`), true, `the templates folder lists "${t.name}.json"; it has ${files.join(', ')}`);
  expectSame(await loadTemplate(`${sb.appName}/templates/${t.name}.json`), t, 'the template read back');
}

function expectPicked(app: HitDesignApp | HitTriageApp, name: string): void {
  const options = Array.from(app.infoView.root.querySelectorAll('select option')).map((o) => o.textContent);
  expect(options.includes(name), true, `the template picker offers "${name}"; it offers ${options.join(', ')}`);
}

async function aiTemplateDetails(app: HitDesignApp | HitTriageApp, name: string): Promise<{[key: string]: any}> {
  const f = hitInfoViewFunctions().find((x) => x.name === 'getHitTemplateDetails')!;
  return await f.apply({view: app.multiView, templateName: name});
}

async function startTriage(sb: Sandbox, app: HitTriageApp, id: string, templateName: string, csv: string,
  columns: {columnSemTypes: {[c: string]: string | null}, columnTypes: {[c: string]: string}}, molColName: string) {
  await sb.writeCampaign({name: id, templateName, filters: [],
    ingest: {type: 'File', query: sb.ingestPath(id), molColName}, status: 'In Progress', createDate: '2026/6/17',
    campaignFields: {Chemist: 'Ann'}, ...columns, rowCount: 1, filteredRowCount: 1,
    permissions: {view: [allUsers], edit: [allUsers]}} as {name: string}, csv);
  await sb.openTriage(app, id);
  await app.saveCampaign(undefined, false);
  return await sb.readCampaign(id);
}

category('Compat: templates', () => {
  test('Hit Design: with no template, the landing offers to create the first one', async () => {
    const sb = await Sandbox.create('tmpl');
    try {
      const app = await sb.hitDesign();
      const empty = app.infoView.root.querySelector<HTMLElement>('[data-u2-name="template-picker"] [data-u2="empty-state"]');
      expect(empty?.textContent?.includes('No templates yet') === true, true, `the picker's empty state; it shows "${empty?.textContent}"`);
      const create = Array.from(empty!.querySelectorAll('button')).find((b) => b.textContent === 'Create template');
      expect(create !== undefined, true, 'the Create template button');
      create!.click();
      await waitUntil(() => grok.shell.v?.name === 'New Template', 'the template editor did not open', 15000);
      grok.shell.v.close();
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Design: the template file is listed, offered, and snapshotted by a new campaign', async () => {
    const sb = await Sandbox.create('tmpl');
    try {
      const t = noComputeTemplate('ZZB3 Design Template', 'ZZB3');
      await sb.writeTemplate(t);
      await expectTemplateFile(sb, t);
      const app = await sb.hitDesign();
      expectPicked(app, t.name);
      const details = await aiTemplateDetails(app, t.name);
      expect(details.key, 'ZZB3', 'key');
      expectSame(details.stages, t.stages, 'stages');
      const id = await startDesignCampaign(app, t, {Chemist: 'Ann'}, 'Created from the file');
      expect(id, 'ZZB3-1', 'the first code of the key');
      const after = await sb.readCampaign(id);
      expect(after.templateName, t.name, 'templateName');
      expectSame(after.template, t, 'template snapshot');
      expect(after.friendlyName, 'Created from the file', 'friendlyName');
      expectSame(after.campaignFields, {Chemist: 'Ann'}, 'campaignFields');
      expect(after.version, 1, 'version');
      const df = await grok.dapi.files.readCsv(after.savePath);
      expectSame(df.columns.names(), ['Molecule', 'Stage', 'V-iD'], 'columns of a new campaign');
      expect(df.get('Stage', 0), 'Idea', 'the first stage is the default');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('PeptiHit: the template file is listed, offered, and snapshotted by a new campaign', async () => {
    const sb = await Sandbox.create('tmpl');
    try {
      const t = {...noComputeTemplate('ZZB3 Peptide Template', 'ZZB3'), stages: ['Design', 'Consideration']};
      await sb.writeTemplate(t);
      await expectTemplateFile(sb, t);
      const app = await sb.peptiHit();
      expectPicked(app, t.name);
      const id = await startDesignCampaign(app, t, {Chemist: 'Ann'}, undefined, 'PEPTIDE1{A.G.L}$$$$V2.0');
      const after = await sb.readCampaign(id);
      expect(after.templateName, t.name, 'templateName');
      expectSame(after.template, t, 'template snapshot');
      const df = await grok.dapi.files.readCsv(after.savePath);
      expectSame(df.columns.names(), ['Helm', 'Molecule', 'Stage', 'V-iD'], 'columns of a new campaign');
      expect(df.get('Stage', 0), 'Design', 'the first stage is the default');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('Hit Triage: the template file opens a campaign without a snapshot and becomes its snapshot', async () => {
    const sb = await Sandbox.create('tmpl');
    try {
      const t = triageTemplate('ZZB3 Triage Template', 'ZZB3');
      await sb.writeTemplate(t);
      await expectTemplateFile(sb, t);
      const app = await sb.hitTriage();
      expectPicked(app, t.name);
      const details = await aiTemplateDetails(app, t.name);
      expect(details.dataSourceType, 'Query', 'dataSourceType');
      expect(details.dataSourceFunction, 'Demo Molecules variable', 'the data source by its friendly name');
      const after = await startTriage(sb, app, 'ZZB3-1', t.name, 'smiles\nCCO\n',
        {columnSemTypes: {smiles: 'Molecule'}, columnTypes: {smiles: 'string'}}, 'smiles');
      expect(after.templateName, t.name, 'templateName');
      expectSame(after.template, t, 'template snapshot');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('PepTriage: the template file opens a campaign and becomes its snapshot', async () => {
    const sb = await Sandbox.create('tmpl');
    try {
      const t = pepTriageTemplate('ZZB3 PepTriage Template', 'ZZB3');
      await sb.writeTemplate(t);
      await expectTemplateFile(sb, t);
      const app = await sb.pepTriage();
      expectPicked(app, t.name);
      const details = await aiTemplateDetails(app, t.name);
      expect(details.dataSourceFunction, 'Demo Peptide Sequences', 'the data source by its friendly name');
      const after = await startTriage(sb, app, 'ZZB3-1', t.name, 'Helm,molfile(Helm)\nPEPTIDE1{A.G.L.V.K}$$$$V2.0,CC(N)C(=O)O\n',
        {columnSemTypes: {'Helm': 'Macromolecule', 'molfile(Helm)': 'Molecule'},
          columnTypes: {'Helm': 'string', 'molfile(Helm)': 'string'}}, 'molfile(Helm)');
      expectSame(after.template, t, 'template snapshot, sequence and molecule columns included');
      expect(after.ingest.molColName, 'molfile(Helm)', 'ingest.molColName');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});
});
