import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {byKind, byName, expectOperable, state, u2test, user} from '@datagrok-libraries/u2/testing';
import {openMergeTableDialog} from '../app/dialogs/merge-table-dialog';
import {HitDesignApp} from '../app/hit-design-app';
import {designCampaign, designCsv, expectSame, noComputeTemplate, Sandbox} from './fixtures';

type Host = Parameters<typeof openMergeTableDialog>[1];

const host = (app: HitDesignApp) => (app as unknown as {_session: {host: Host}})._session.host;

// a table saved after a filter or a search read its molecules carries Chem's twin
const TWIN_HEADER = '~Molecule.canonicalSmiles,Molecule,V-iD,Note';

const OPTIONS = ['merge-by', 'molecule-column', 'vid-column', 'add-rows', 'run-compute', 'on-clash', 'auto-merge'];

async function mergeDialog(sb: Sandbox, mergeConfig?: object): Promise<{app: HitDesignApp, dialog: HTMLElement,
  filePath: string}> {
  const t = noComputeTemplate('ZZB3 Merge', 'ZZB3');
  await sb.writeTemplate(t);
  const filePath = `System:AppData/HitTriage/${sb.appName}/merge-source.csv`;
  await grok.dapi.files.writeAsText(filePath, 'Molecule,V-iD,Note\nOCC,,merged\n');
  await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, mergeConfig ? {mergeConfig} : {}), designCsv);
  const app = await sb.hitDesign();
  await sb.openDesign(app, 'ZZB3-1');
  await openMergeTableDialog(app, host(app));
  return {app, dialog: byKind(document, 'dialog').find((d) => state(d).name === 'merge-table')!, filePath};
}

const checkbox = (dialog: HTMLElement, name: string) => byName(dialog, name).querySelector('input')!;

async function typePath(dialog: HTMLElement, filePath: string): Promise<void> {
  user.type(byName(dialog, 'source-file'), filePath);
  await waitUntil(() => state(byName(dialog, 'molecule-column')).value === 'Molecule',
    'the typed share file was never read', 15000);
}

category('Campaign view: merge', () => {
  test('Merge with table neither offers nor picks a hidden column of the incoming table', u2test(async () => {
    const sb = await Sandbox.create('merge');
    try {
      const t = noComputeTemplate('ZZB3 Merge', 'ZZB3');
      await sb.writeTemplate(t);
      const filePath = `System:AppData/HitTriage/${sb.appName}/merge-source.csv`;
      await grok.dapi.files.writeAsText(filePath, `${TWIN_HEADER}\nOCC,OCC,,merged\n`);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {mergeConfig: {mode: 'smiles', filePath,
        addNewRows: false, runComputeOnNewRows: false, clashStrategy: 'skip'}}), designCsv);
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      await openMergeTableDialog(app, host(app));
      const dialog = byKind(document, 'dialog').find((d) => state(d).name === 'merge-table')!;
      const offered = (name: string) => state(byName(dialog, name)).items?.list.map((i) => i.label) ?? [];
      await waitUntil(() => offered('molecule-column').length > 1, 'the dialog never listed the file\'s columns',
        15000);
      const columns = ['(none)', 'Molecule', 'V-iD', 'Note'];
      expectSame(offered('molecule-column'), columns, 'Molecule column (incoming) choices');
      expectSame(offered('vid-column'), columns, 'V-iD column (incoming) choices');
      expectSame([state(byName(dialog, 'molecule-column')).value, state(byName(dialog, 'vid-column')).value],
        ['Molecule', 'V-iD'], 'the columns picked by default');
      user.click(byName(dialog, 'Cancel'));
      await waitUntil(() => !dialog.isConnected, 'CANCEL left the dialog open');
    } finally {
      await sb.dispose();
    }
  }), {timeout: 120000});

  test('a merge carries no hidden column of the incoming table into the campaign: Auto-merge on open, a dropped table',
    async () => {
      const sb = await Sandbox.create('merge');
      try {
        const t = noComputeTemplate('ZZB3 Merge', 'ZZB3');
        await sb.writeTemplate(t);
        const filePath = `System:AppData/HitTriage/${sb.appName}/merge-source.csv`;
        await grok.dapi.files.writeAsText(filePath, `${TWIN_HEADER}\nOCC,OCC,,merged\n`);
        await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {mergeConfig: {mode: 'smiles', filePath,
          molColName: 'Molecule', addNewRows: false, runComputeOnNewRows: false, clashStrategy: 'overwrite',
          autoMergeOnOpen: true}}), designCsv);
        const app = await sb.hitDesign();
        await sb.openDesign(app, 'ZZB3-1');
        const hidden = () => app.dataFrame!.columns.names().filter((n) => n.startsWith('~'));
        expect(app.dataFrame!.get('Note', 0), 'merged', 'the file\'s value in the matched row (OCC is CCO)');
        expectSame(hidden(), [], 'hidden columns of the campaign\'s table after Auto-merge on open');
        const dropped = DG.DataFrame.fromCsv(`${TWIN_HEADER}\nCCN,CCN,,dropped\n`);
        await grok.data.detectSemanticTypes(dropped);
        await app.handleJoiningDataframe(dropped);
        expect(app.dataFrame!.get('Note', app.dataFrame!.rowCount - 1), 'dropped', 'the dropped table\'s row added');
        expectSame(hidden(), [], 'hidden columns of the campaign\'s table after a dropped table');
      } finally {
        await sb.dispose();
      }
    }, {timeout: 120000});

  test('Merge with table: an option that is off says why, as the file, the key and Add unmatched rows decide',
    u2test(async () => {
      const sb = await Sandbox.create('merge');
      try {
        const {dialog, filePath} = await mergeDialog(sb);
        const why = () => Object.fromEntries(OPTIONS.map((n) => [n, state(byName(dialog, n)).disabledReason ?? null]));
        const ok = () => byName(dialog, 'OK').getAttribute('aria-description');
        const first = 'Pick a source file first';
        expectOperable(dialog);
        expectSame(why(), Object.fromEntries(OPTIONS.map((n) => [n, first])), 'why each option is off before a file');
        expect(ok(), 'Pick a source file to merge from.', 'why OK is off before a file');
        await typePath(dialog, filePath);
        expectSame(why(), Object.fromEntries(OPTIONS.map((n) => [n, null])), 'what is off with a share file, by V-iD');
        user.pick(byName(dialog, 'merge-by'), 'By canonical SMILES');
        user.click(checkbox(dialog, 'add-rows'));
        expectSame(why(), {...Object.fromEntries(OPTIONS.map((n) => [n, null])),
          'vid-column': 'Used only when merging by V-iD', 'run-compute': 'Turn on Add unmatched rows first'},
        'what is off by canonical SMILES, with no unmatched rows added');
        user.pick(byName(dialog, 'molecule-column'), '(none)');
        expectSame(why(), {...Object.fromEntries(OPTIONS.map((n) => [n, null])),
          'vid-column': 'Used only when merging by V-iD', 'add-rows': 'Pick the molecule column (incoming) first',
          'run-compute': 'Pick the molecule column (incoming) first'}, 'what is off with no molecule column');
        expect(ok(), 'Pick the molecule column from the incoming table.', 'why OK is off with no molecule column');
        expectOperable(dialog);
        user.click(byName(dialog, 'Cancel'));
        await waitUntil(() => !dialog.isConnected, 'CANCEL left the dialog open');
      } finally {
        await sb.dispose();
      }
    }), {timeout: 120000});

  test('Merge with table: Source file offers the platform\'s Files; a share file is kept to merge again on open',
    u2test(async () => {
      const sb = await Sandbox.create('merge');
      try {
        const {app, dialog, filePath} = await mergeDialog(sb);
        const filesDialog = () => DG.Dialog.getOpenDialogs().find((d) => d.title === 'Select a file');
        user.click(byName(dialog, 'Add file from Files'));
        await waitUntil(() => filesDialog() !== undefined, 'Add file from Files opened no "Select a file"', 15000);
        filesDialog()!.close();
        await waitUntil(() => filesDialog() === undefined, '"Select a file" stayed open', 10000);
        expect(state(byName(dialog, 'source-file')).value ?? null, null, 'nothing picked: the box stays empty');
        await typePath(dialog, filePath);
        user.pick(byName(dialog, 'merge-by'), 'By canonical SMILES');
        user.click(checkbox(dialog, 'auto-merge'));
        user.click(byName(dialog, 'OK'));
        await waitUntil(() => !dialog.isConnected, 'OK left the dialog open', 30000);
        await sb.savesDone();
        expect(app.dataFrame!.get('Note', 0), 'merged', 'the file\'s value in the matched row (OCC is CCO)');
        const kept = (await sb.readCampaign('ZZB3-1')).mergeConfig;
        expectSame([kept?.filePath, kept?.autoMergeOnOpen], [filePath, true], 'the file the campaign merges on open');
      } finally {
        await sb.dispose();
      }
    }), {timeout: 120000});
});
