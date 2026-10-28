/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {_package} from '../package';
import * as api from '../package-api';
import {HitDesignApp} from '../app/hit-design-app';
import {HitDesignTemplate} from '../app/types';
import {obfuscateSmiles} from '../app/utils/molreg';
import {mergeIntoCampaign} from '../app/utils/merge-table';
import {addRows, allUsers, countQueryCalls, designCampaign, designColumns, designCsv, expectSame, hasBalloon, hitdesignSql,
  noComputeTemplate, otherUserId, Sandbox, sqlText, startDesignCampaign, startTriageCampaign, triageTemplate,
  unknownGroupId} from './fixtures';

const acquire = (sb: Sandbox, id: string, lockedBy: string) => api.queries.acquireCampaignLock(sb.appName, id, lockedBy);
const release = (sb: Sandbox, id: string) => api.queries.releaseCampaignLock(sb.appName, id);

async function lockRow(sb: Sandbox, id: string): Promise<{ttl: number, lockedBy: string} | null> {
  const df = await hitdesignSql('SELECT EXTRACT(EPOCH FROM expires_at - NOW()) AS ttl, locked_by FROM hitdesign.campaign_locks ' +
    `WHERE app_name = ${sqlText(sb.appName)} AND campaign_id = ${sqlText(id)}`);
  return df.rowCount ? {ttl: df.get('ttl', 0), lockedBy: df.get('locked_by', 0)} : null;
}

/** A campaign on {@link designCsv} at version 1, open in [sessions] separate app instances. */
async function openShared(sb: Sandbox, sessions: number, extra: {[key: string]: unknown} = {}): Promise<HitDesignApp[]> {
  const t = noComputeTemplate('ZZB3 Shared', 'ZZB3');
  await sb.writeTemplate(t);
  await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {version: 1, ...extra}), designCsv);
  const apps: HitDesignApp[] = [];
  for (let i = 0; i < sessions; i++) {
    const app = await sb.hitDesign();
    await sb.openDesign(app, 'ZZB3-1');
    apps.push(app);
  }
  return apps;
}

/** The server table of the shared campaign, keyed by canonical SMILES; a molecule held twice fails. */
async function serverRows(sb: Sandbox): Promise<Map<string, {[col: string]: unknown}>> {
  const df = await grok.dapi.files.readCsv(sb.savePath('ZZB3-1'));
  const rows = new Map<string, {[col: string]: unknown}>();
  for (let i = 0; i < df.rowCount; i++) {
    const smiles = _package.convertToSmiles(df.get('Molecule', i));
    if (rows.has(smiles))
      throw new Error(`the server table holds ${smiles} twice`);
    rows.set(smiles, Object.fromEntries(df.columns.names().map((c) => [c, df.get(c, i)])));
  }
  return rows;
}

const otherCsv = 'Molecule,Stage,V-iD,Note\nCCCC,Idea,,other\n';

/** A second campaign of the shared template, as files: what a test checks no save of the first one touched. */
async function writeOther(sb: Sandbox): Promise<() => Promise<void>> {
  await sb.writeCampaign(designCampaign(sb, 'ZZB3-2', noComputeTemplate('ZZB3 Shared', 'ZZB3'),
    {friendlyName: 'Other campaign', version: 1}), otherCsv);
  const json = await sb.readCampaignText('ZZB3-2');
  const csv = await grok.dapi.files.readAsText(sb.savePath('ZZB3-2'));
  return async () => {
    expect(await sb.readCampaignText('ZZB3-2'), json, 'ZZB3-2\'s campaign.json untouched');
    expect(await grok.dapi.files.readAsText(sb.savePath('ZZB3-2')), csv, 'ZZB3-2\'s table untouched');
  };
}

/** The view state of a table view over [csv] with a scatter plot of [column], after [edit] (say, closing its grid). */
function layoutOf(csv: string, column: string, edit: (tv: DG.TableView) => void = () => {}): string {
  const tv = grok.shell.addTableView(DG.DataFrame.fromCsv(csv));
  try {
    tv.addViewer(DG.VIEWER.SCATTER_PLOT, {xColumnName: column, yColumnName: column});
    edit(tv);
    return tv.saveLayout().viewState;
  } finally {
    const df = tv.dataFrame;
    tv.close();
    grok.shell.closeTable(df);
  }
}

type SessionInternals = {_session?: {openTracker(): {getButton(text: string): HTMLButtonElement | undefined}}};

function u2Button(name: string, text: string): HTMLButtonElement | undefined {
  return Array.from(document.querySelectorAll<HTMLButtonElement>(`[data-u2-name="${name}"] button`))
    .find((b) => b.textContent?.trim() === text);
}

/** A molfile as a sketcher draws it: a chain at angles no layout picks, so a save that redrew it would show. */
function drawn(atoms: [string, number, number][]): string {
  const n = (v: number, width: number) => String(v).padStart(width);
  return ['', '  ZZB3 drawn', '', `${n(atoms.length, 3)}${n(atoms.length - 1, 3)}  0  0  0  0  0  0  0  0999 V2000`,
    ...atoms.map(([el, x, y]) => `${x.toFixed(4).padStart(10)}${y.toFixed(4).padStart(10)}    0.0000 ${el.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`),
    ...atoms.slice(1).map((_, i) => `${n(i + 1, 3)}${n(i + 2, 3)}  1  0  0  0  0`), 'M  END', ''].join('\n');
}

const drawnEthanol = drawn([['C', 0.5, 1.25], ['C', 1.75, 0.4], ['O', 1.9, -1.1]]);
const drawnPropylamine = drawn([['C', 0, 0], ['C', 1.2, 0.9], ['C', 1.5, -0.6], ['N', 2.9, -0.3]]);
/** What OpenChemLib's sketcher gives on OK with nothing drawn. */
const emptySketch = '\nActelion Java MolfileCreator 1.0\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n';

const units = (app: HitDesignApp) => app.dataFrame!.col('Molecule')!.meta.units;

/** The template with the JS-only Chemical Properties compute, MW alone. */
const withMolWeight = (t: HitDesignTemplate): HitDesignTemplate => ({...t, compute: {...t.compute, functions: [
  {name: 'addChemPropertiesColumns', package: 'Chem', args: {MW: true, HBA: false, HBD: false, logP: false, logS: false,
    PSA: false, rotatableBonds: false, stereoCenters: false, moleculeCharge: false, molecularFormula: false}}]}});

/** The sketcher a campaign opens on its empty first row, closed once it is up. */
async function closeSketcher(): Promise<void> {
  await waitUntil(() => DG.Dialog.getOpenDialogs().length > 0, 'the sketcher never opened on the empty row', 30000);
  for (const d of DG.Dialog.getOpenDialogs())
    d.close();
}

category('Save and lock', () => {
  test('the campaign lock is exclusive and lasts 30 s', async () => {
    const sb = await Sandbox.create('lock');
    try {
      const first = await acquire(sb, 'ZZB3-1', 'session A');
      expect(first.rowCount, 1, 'the first acquire takes the lock');
      const second = await acquire(sb, 'ZZB3-1', 'session B');
      expect(second.rowCount, 0, 'a second acquire is refused while the lock is held');
      const row = await lockRow(sb, 'ZZB3-1');
      expect(row?.lockedBy, 'session A', 'locked_by');
      expect(row!.ttl > 25 && row!.ttl <= 30.5, true, `the lock expires in 30 s, it has ${row!.ttl} s left`);
    } finally {
      await sb.dispose();
    }
  });

  // The query's DELETE of the expired row and its INSERT share one snapshot, so the first call only clears the row and
  // the save's retry 200 ms later takes the lock; the outcome pinned here holds with or without that extra round.
  test('an expired lock is taken over by the next attempt at the latest', async () => {
    const sb = await Sandbox.create('lock');
    try {
      await acquire(sb, 'ZZB3-1', 'session A');
      await hitdesignSql(`WITH u AS (UPDATE hitdesign.campaign_locks SET expires_at = NOW() - INTERVAL '1 second' ` +
        `WHERE app_name = ${sqlText(sb.appName)} RETURNING 1) SELECT COUNT(*) AS n FROM u`);
      let taken = await acquire(sb, 'ZZB3-1', 'session B');
      if (taken.rowCount === 0)
        taken = await acquire(sb, 'ZZB3-1', 'session B');
      expect(taken.rowCount, 1, 'the expired lock is replaced');
      expect((await lockRow(sb, 'ZZB3-1'))?.lockedBy, 'session B', 'the new holder');
    } finally {
      await sb.dispose();
    }
  });

  test('releasing the lock logs the save time that GetLastModified returns', async () => {
    const sb = await Sandbox.create('lock');
    try {
      expect(await sb.lastModified('ZZB3-1'), null, 'no save logged yet');
      await acquire(sb, 'ZZB3-1', 'session A');
      const released = (await release(sb, 'ZZB3-1')).col('updated_at')!.getString(0);
      expect(await sb.lastModified('ZZB3-1'), released, 'GetLastModified');
      expect(await lockRow(sb, 'ZZB3-1'), null, 'the lock row is gone');
      const again = await release(sb, 'ZZB3-1');
      expect(again.rowCount, 0, 'releasing a lock nobody holds logs nothing');
      expect(await sb.lastModified('ZZB3-1'), released, 'GetLastModified after the empty release');
    } finally {
      await sb.dispose();
    }
  });

  test('two saves in sequence raise the version by one each and log each save', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [app] = await openShared(sb, 1);
      await app.saveCampaign(false);
      expect((await sb.readCampaign('ZZB3-1')).version, 2, 'version after the first save');
      const t1 = await sb.lastModified('ZZB3-1');
      await app.saveCampaign(false);
      expect((await sb.readCampaign('ZZB3-1')).version, 3, 'version after the second save');
      const t2 = await sb.lastModified('ZZB3-1');
      expect(t1 != null && t2 != null && Date.parse(t2) > Date.parse(t1), true, `the save time moves on: ${t1} → ${t2}`);
      expect(await lockRow(sb, 'ZZB3-1'), null, 'no lock left behind');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a save waits while another session holds the lock', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [app] = await openShared(sb, 1);
      await acquire(sb, 'ZZB3-1', 'another session');
      await countQueryCalls('AcquireCampaignLock', null, async (attempts) => {
        const save = app.saveCampaign(false);
        await waitUntil(() => attempts.length >= 3, `the save asked for the lock ${attempts.length} times`, 15000);
        expect((await sb.readCampaign('ZZB3-1')).version, 1, 'nothing is written while the lock is held');
        await release(sb, 'ZZB3-1');
        await save;
      });
      expect((await sb.readCampaign('ZZB3-1')).version, 2, 'the save went through after the release');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a stale save adds the other session\'s new rows by canonical SMILES', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [a, b] = await openShared(sb, 2);
      const score = a.dataFrame!.columns.addNewInt('Score');
      score.set(0, 1);
      score.set(1, 2);
      addRows(a.dataFrame!, [{'Molecule': 'CCN', 'Stage': 'Idea', 'Score': 3}, {'Molecule': 'OCC', 'Stage': 'Idea', 'Score': 9}]);
      await a.saveCampaign(false);
      addRows(b.dataFrame!, [{Molecule: 'CCCl', Stage: 'Idea'}]);
      await b.saveCampaign(false);

      expect((await sb.readCampaign('ZZB3-1')).version, 3, 'the stale save goes past the other one');
      const rows = await serverRows(sb);
      expect(Array.from(rows.keys()).sort().join(' '), ['CCCl', 'CCN', 'CCO', 'c1ccccc1'].sort().join(' '),
        'rows: ours, plus the other session\'s new CCN; its OCC is our CCO');
      expect(rows.get('CCO')!['Score'], 9, 'a new column is copied into the rows we share; of two rows of one molecule the later wins');
      expect(rows.get('c1ccccc1')!['Score'], 2, 'Score of benzene');
      expect(rows.get('CCN')!['Score'], 3, 'Score of the added row');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('a view-only save is refused and writes nothing', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [app] = await openShared(sb, 1, {authorUserId: otherUserId, authorUserFriendlyName: 'Someone Else',
        permissions: {view: [allUsers], edit: [unknownGroupId]}});
      expect(app.canEdit, false, 'the campaign is view-only');
      const json = await sb.readCampaignText('ZZB3-1');
      const csv = await grok.dapi.files.readAsText(sb.savePath('ZZB3-1'));
      app.dataFrame!.set('Note', 0, 'a local edit');
      await app.saveCampaign(false);
      expect(await sb.readCampaignText('ZZB3-1'), json, 'campaign.json unchanged');
      expect(await grok.dapi.files.readAsText(sb.savePath('ZZB3-1')), csv, 'the table unchanged');
      expect(hasBalloon('You do not have permission to modify this campaign'), true, 'the refusal is shown');
      expect(await lockRow(sb, 'ZZB3-1'), null, 'no lock was taken');
      expect(app.dataFrame!.get('Note', 0), 'a local edit', 'the local edit stays');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a new campaign takes the next free number of its template key', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = noComputeTemplate('ZZB3 Codes', 'ZZB3');
      await sb.writeTemplate(t);
      const app = await sb.hitDesign();
      const folder = `${sb.appName}/campaigns`;
      expect(await app.getNewCampaignName(folder, 'ZZB3'), 'ZZB3-1', 'with no campaigns of the key');
      for (const id of ['ZZB3-1', 'ZZB3-7', 'ZZB3-x', 'OTHER-9', 'ZZB3X-12'])
        await sb.writeCampaign(designCampaign(sb, id, t), designCsv);
      expect(await app.getNewCampaignName(folder, 'ZZB3'), 'ZZB3-8', 'the highest number + 1');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('two sessions creating from one template get distinct codes', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = noComputeTemplate('ZZB3 Race', 'ZZB3');
      await sb.writeTemplate(t);
      const a = await sb.hitDesign();
      const b = await sb.hitDesign();
      const ids = await Promise.all([startDesignCampaign(a, t, {Chemist: 'A'}, 'Race A', 'CCO'),
        startDesignCampaign(b, t, {Chemist: 'B'}, 'Race B', 'CCN')]);
      expectSame([...ids].sort(), ['ZZB3-1', 'ZZB3-2'], 'two codes of the key');
      expect(await lockRow(sb, ids[0]) === null && await lockRow(sb, ids[1]) === null, true, 'no lock is left');
      for (const [i, id] of ids.entries()) {
        const saved = await sb.readCampaign(id);
        expect(saved.friendlyName, i === 0 ? 'Race A' : 'Race B', `${id} name`);
        const df = await grok.dapi.files.readCsv(saved.savePath);
        expect(df.rowCount, 1, `${id} holds only its own row`);
        expect(_package.convertToSmiles(df.get('Molecule', 0)), i === 0 ? 'CCO' : 'CCN', `${id} molecule`);
      }
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('two Hit Triage sessions starting from one template get distinct codes', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = {...triageTemplate('ZZB3 Triage Race', 'ZZB3'), dataSourceType: 'File' as const};
      await sb.writeTemplate(t);
      const [a, b] = [await sb.hitTriage(), await sb.hitTriage()];
      const ids = await Promise.all([startTriageCampaign(a, t, 'CCO', 'Race A'), startTriageCampaign(b, t, 'CCN', 'Race B')]);
      expectSame([...ids].sort(), ['ZZB3-1', 'ZZB3-2'], 'two codes of the key');
      expect(await lockRow(sb, ids[0]) === null && await lockRow(sb, ids[1]) === null, true, 'no lock is left');
      for (const [i, id] of ids.entries()) {
        const saved = await sb.readCampaign(id);
        expect(saved.friendlyName, i === 0 ? 'Race A' : 'Race B', `${id} name`);
        const df = await grok.dapi.files.readCsv(saved.ingest.query);
        expect(df.rowCount, 1, `${id} holds only its own row`);
        expect(_package.convertToSmiles(df.get('smiles', 0)), i === 0 ? 'CCO' : 'CCN', `${id} molecule`);
      }
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('a campaign that fails to open leaves the open one in place: its next save writes itself', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [app] = await openShared(sb, 1);
      const t = noComputeTemplate('ZZB3 Shared', 'ZZB3');
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-2', t, {status: 'Broken', savePath: sb.savePath('ZZB3-2', 'elsewhere')}),
        designCsv);
      await app.infoView.setCampaign('ZZB3-2');
      expect(hasBalloon('Campaign ZZB3-2 could not be opened'), true, 'the failed open says so');
      expect(app.campaign?.name, 'ZZB3-1', 'the app still holds the open campaign');
      app.dataFrame!.set('Note', 0, 'after the failed open');
      await app.saveCampaign(false);
      expect((await serverRows(sb)).get('CCO')!['Note'], 'after the failed open', 'the open campaign\'s table took the edit');
      const own = await sb.readCampaign('ZZB3-1');
      expect(own.savePath, sb.savePath('ZZB3-1'), 'its campaign.json keeps its own table');
      expect(own.status, 'In Progress', 'and its own status');
      expect(await grok.dapi.files.exists(sb.savePath('ZZB3-2', 'elsewhere')), false, 'nothing written where the broken one points');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Hit Triage: a campaign that fails to open leaves the open one saving as itself', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = triageTemplate('ZZB3 Triage Open', 'ZZB3');
      await sb.writeTemplate(t);
      const campaign = (id: string, status: string, templateName: string) => ({name: id, templateName, filters: [],
        ingest: {type: 'File', query: sb.ingestPath(id), molColName: 'smiles'}, status, createDate: '2024/10/17',
        campaignFields: {Chemist: id}, columnSemTypes: {smiles: 'Molecule'}, columnTypes: {smiles: 'string'},
        permissions: {view: [allUsers], edit: [allUsers]}});
      await sb.writeCampaign(campaign('ZZB3-1', 'In Progress', t.name), 'smiles\nCCO\n');
      await sb.writeCampaign(campaign('ZZB3-2', 'Broken', 'ZZB3 Gone Template'), 'smiles\nCCN\n');
      const app = await sb.hitTriage();
      await sb.openTriage(app, 'ZZB3-1');
      await app.infoView.setCampaign('ZZB3-2');
      expect(hasBalloon('Campaign ZZB3-2 could not be opened'), true, 'the failed open says so');
      expect(app.campaign?.name, 'ZZB3-1', 'the app still holds the open campaign');
      await app.saveCampaign(undefined, false);
      const own = await sb.readCampaign('ZZB3-1');
      expect(own.status, 'In Progress', 'its own status, not the broken one\'s');
      expectSame(own.campaignFields, {Chemist: 'ZZB3-1'}, 'its own fields');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('deleting a campaign removes its last save, its lock and its V-iD links', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = noComputeTemplate('ZZB3 Doomed', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t), designCsv);
      await acquire(sb, 'ZZB3-1', 'an earlier save');
      await release(sb, 'ZZB3-1');
      await acquire(sb, 'ZZB3-1', 'a session that went away');
      await api.queries.addMolecule(obfuscateSmiles(_package.convertToSmiles(sb.ownMolecule('CCCCCCCCCO'))), sb.appName,
        'ZZB3-1', 'test');
      const where = `WHERE app_name = ${sqlText(sb.appName)} AND campaign_id = 'ZZB3-1'`;
      const rows = async () => (await hitdesignSql(`SELECT (SELECT COUNT(*) FROM hitdesign.update_logs ${where}) + ` +
        `(SELECT COUNT(*) FROM hitdesign.campaign_locks ${where}) + (SELECT COUNT(*) FROM hitdesign.campaign_vids ${where}) AS n`))
        .get('n', 0);
      expect(await rows(), 3, 'a save logged, a lock held, a V-iD linked');
      const app = await sb.hitDesign();
      await app.infoView.deleteCampaignAndRefresh(sb.appName, 'ZZB3-1');
      expect(await _package.files.exists(sb.path('ZZB3-1')), false, 'the campaign\'s files are gone');
      expect(await rows(), 0, 'and its hitdesign rows');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('dropping a lock without a save lets go only of the caller\'s own lock', async () => {
    const sb = await Sandbox.create('save');
    try {
      await acquire(sb, 'ZZB3-1', 'session A');
      await api.queries.dropCampaignLock(sb.appName, 'ZZB3-1', 'session B');
      expect((await lockRow(sb, 'ZZB3-1'))?.lockedBy, 'session A', 'another session\'s lock stays');
      await api.queries.dropCampaignLock(sb.appName, 'ZZB3-1', 'session A');
      expect(await lockRow(sb, 'ZZB3-1'), null, 'the holder\'s drop takes it');
      expect(await sb.lastModified('ZZB3-1'), null, 'and logs no save');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 60000});

  test('a stale save keeps the other session\'s edits to the rows both hold', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [a, b] = await openShared(sb, 2);
      a.dataFrame!.set('Stage', 0, 'Tested');
      await a.saveCampaign(false);
      b.dataFrame!.set('Note', 1, 'edited by B');
      b.dataFrame!.set('Note', 0, 'B wins a cell both edited');
      a.dataFrame!.set('Note', 0, 'A edited it first');
      await a.saveCampaign(false);
      await b.saveCampaign(false);
      const rows = await serverRows(sb);
      expect(rows.get('c1ccccc1')!['Note'], 'edited by B', 'the stale session\'s own edit is saved');
      expect(rows.get('CCO')!['Stage'], 'Tested', 'the other session\'s Stage edit of CCO is kept');
      expect(rows.get('CCO')!['Note'], 'B wins a cell both edited', 'a cell both edited takes the later save');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('a refused view-only save does not log a change', async () => {
    const sb = await Sandbox.create('save');
    try {
      await acquire(sb, 'ZZB3-1', 'an earlier save');
      await release(sb, 'ZZB3-1');
      const before = await sb.lastModified('ZZB3-1');
      const [app] = await openShared(sb, 1, {authorUserId: otherUserId, permissions: {view: [allUsers], edit: [unknownGroupId]}});
      app.dataFrame!.set('Note', 0, 'a local edit');
      await app.saveCampaign(false);
      expect(hasBalloon('You do not have permission to modify this campaign'), true, 'the save was refused');
      expect(await sb.lastModified('ZZB3-1'), before, 'update_logs stays as it was');
      expect(await lockRow(sb, 'ZZB3-1'), null, 'no lock was taken');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a save whose table cannot be written lets the lock go and says why', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [app] = await openShared(sb, 1);
      const before = await sb.lastModified('ZZB3-1');
      const toCsv = DG.DataFrame.prototype.toCsv;
      DG.DataFrame.prototype.toCsv = () => {
        throw new Error('ZZB3 unwritable table');
      };
      let saved: unknown;
      try {
        saved = await app.saveCampaign(false);
      } finally {
        DG.DataFrame.prototype.toCsv = toCsv;
      }
      expect(saved === undefined, true, 'the save answers that nothing was saved');
      expect(hasBalloon('Failed to save the campaign: ZZB3 unwritable table'), true, 'the failure is shown');
      expect(app.isSaving, false, 'the save has ended');
      expect(await lockRow(sb, 'ZZB3-1'), null, 'the lock is let go');
      expect(await sb.lastModified('ZZB3-1'), before, 'nothing written, no save logged');
      expect((await sb.readCampaign('ZZB3-1')).version, 1, 'campaign.json untouched');
      expect(await app.saveCampaign(false) !== undefined, true, 'the next save goes through');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a save whose table the share refuses lets the lock go and logs no save', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [app] = await openShared(sb, 1);
      const before = await sb.lastModified('ZZB3-1');
      const files = Object.getPrototypeOf(grok.dapi.files);
      const writeAsText = files.writeAsText;
      files.writeAsText = function(file: DG.FileInfo | string, data: string) {
        return file === sb.savePath('ZZB3-1') ? Promise.reject(new Error('ZZB3 share refused the table')) :
          writeAsText.call(this, file, data);
      };
      let saved: unknown;
      try {
        saved = await app.saveCampaign(false);
      } finally {
        files.writeAsText = writeAsText;
      }
      expect(saved === undefined, true, 'the save answers that nothing was saved');
      expect(hasBalloon('Failed to save the campaign: ZZB3 share refused the table'), true, 'the failure is shown');
      expect(await lockRow(sb, 'ZZB3-1'), null, 'the lock is let go');
      expect(await sb.lastModified('ZZB3-1'), before, 'nothing written, no save logged');
      expect((await sb.readCampaign('ZZB3-1')).version, 1, 'campaign.json untouched');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('saves asked for while one runs fold into one more write', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [app] = await openShared(sb, 1);
      const {calls} = await countQueryCalls('ReleaseCampaignLock', null, async () => {
        const first = app.saveCampaign(false);
        app.dataFrame!.set('Note', 0, 'second save');
        const second = app.saveCampaign(false);
        app.dataFrame!.set('Note', 1, 'third save');
        const third = app.saveCampaign(false);
        expect(second === third, true, 'the third save is folded into the second');
        expect(app.isSaving, true, 'saving while the first runs');
        await first;
        expect(app.isSaving, true, 'still saving: the folded save runs after the first');
        await second;
        expect(app.isSaving, false, 'done after the last');
      });
      expect(calls.length, 2, 'two writes for three saves');
      expect((await sb.readCampaign('ZZB3-1')).version, 3, 'the version went up by two');
      const rows = await serverRows(sb);
      expect(rows.get('CCO')!['Note'], 'second save', 'the folded save wrote the second edit');
      expect(rows.get('c1ccccc1')!['Note'], 'third save', 'and the third');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a save asked for while another campaign opens writes the campaign it was asked for', async () => {
    const sb = await Sandbox.create('save');
    const files = Object.getPrototypeOf(grok.dapi.files) as {readCsv: (...args: unknown[]) => Promise<DG.DataFrame>};
    const readCsv = files.readCsv;
    try {
      const [app] = await openShared(sb, 1);
      const untouched = await writeOther(sb);
      let release!: () => void;
      const gate = new Promise<void>((resolve) => release = resolve);
      let reading = false;
      files.readCsv = async function(this: unknown, ...args: unknown[]) {
        if (args[0] === sb.savePath('ZZB3-2')) {
          reading = true;
          await gate;
        }
        return readCsv.apply(this, args);
      };
      app.dataFrame!.set('Note', 0, 'asked for mid-open');
      const opening = app.infoView.setCampaign('ZZB3-2');
      await waitUntil(() => reading, 'the open never read ZZB3-2\'s table', 30000);
      const saved = await app.saveCampaign(false);
      release();
      await opening;
      expect(saved?.name, 'ZZB3-1', 'the save answers the campaign it was asked for');
      const own = await sb.readCampaign('ZZB3-1');
      expect(own.friendlyName, 'Friendly ZZB3-1', 'ZZB3-1\'s campaign.json is its own');
      expect(own.savePath, sb.savePath('ZZB3-1'), 'and points at its own table');
      expect((await serverRows(sb)).get('CCO')!['Note'], 'asked for mid-open', 'ZZB3-1\'s table took the edit');
      await untouched();
      expect(app.campaign?.name, 'ZZB3-2', 'the open went on');
    } finally {
      files.readCsv = readCsv;
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('a save queued behind one waiting for the lock writes its campaign though another opened meanwhile', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [app] = await openShared(sb, 1);
      const untouched = await writeOther(sb);
      await acquire(sb, 'ZZB3-1', 'another session');
      app.dataFrame!.set('Note', 0, 'first save');
      const first = app.saveCampaign(false);
      app.dataFrame!.set('Note', 1, 'queued save');
      const queued = app.saveCampaign(false);
      await app.infoView.setCampaign('ZZB3-2');
      expect(app.campaign?.name, 'ZZB3-2', 'ZZB3-2 opened while the saves waited');
      await release(sb, 'ZZB3-1');
      expect((await first)?.name, 'ZZB3-1', 'the first save');
      expect((await queued)?.name, 'ZZB3-1', 'the queued save wrote the campaign it was asked for');
      expect((await sb.readCampaign('ZZB3-1')).version, 3, 'two saves of ZZB3-1');
      const rows = await serverRows(sb);
      expect(rows.get('CCO')!['Note'], 'first save', 'the first edit');
      expect(rows.get('c1ccccc1')!['Note'], 'queued save', 'the queued edit');
      await untouched();
      expect(app.campaign?.version, 1, 'the open campaign is ZZB3-2 as it was read');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('a lock query that fails says why; only a lock another session keeps blames someone else', async () => {
    const sb = await Sandbox.create('save');
    const someoneElse = 'Someone else might be editing the campaign';
    try {
      const [app] = await openShared(sb, 1);
      const lock = app as unknown as {aquireAndWaitCampaignLock: (id: string) => Promise<{aquired: boolean}>};
      try {
        lock.aquireAndWaitCampaignLock = async () => {
          throw new Error('ZZB3 query refused');
        };
        expect(await app.saveCampaign(false) === undefined, true, 'nothing saved');
        expect(hasBalloon('Failed to save the campaign: ZZB3 query refused'), true, 'the failure is shown with its cause');
        expect(hasBalloon(someoneElse), false, 'nobody else is blamed');
        lock.aquireAndWaitCampaignLock = async () => ({aquired: false});
        expect(await app.saveCampaign(false) === undefined, true, 'nothing saved while the lock is kept');
        expect(hasBalloon(someoneElse), true, 'a lock another session kept blames someone else');
      } finally {
        delete (lock as {aquireAndWaitCampaignLock?: unknown}).aquireAndWaitCampaignLock;
      }
      expect((await sb.readCampaign('ZZB3-1')).version, 1, 'campaign.json untouched');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a stale save keeps its own cells of a molecule two rows hold', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = noComputeTemplate('ZZB3 Shared', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {version: 1, rowCount: 3}),
        'Molecule,Stage,V-iD,Note\nCCO,Idea,,cco-1\nOCC,Idea,,cco-2\nc1ccccc1,Made,,benzene\n');
      const [a, b] = [await sb.hitDesign(), await sb.hitDesign()];
      await sb.openDesign(a, 'ZZB3-1');
      await sb.openDesign(b, 'ZZB3-1');
      a.dataFrame!.set('Note', 1, 'A edited the second CCO');
      a.dataFrame!.set('Note', 2, 'A edited benzene');
      await a.saveCampaign(false);
      b.dataFrame!.set('Stage', 0, 'Made');
      await b.saveCampaign(false);
      const df = await grok.dapi.files.readCsv(sb.savePath('ZZB3-1'));
      expect(df.rowCount, 3, 'rows');
      expect(df.get('Note', 0), 'cco-1', 'the first CCO row');
      expect(df.get('Note', 1), 'cco-2', 'the second CCO row keeps ours: which of the other session\'s rows is it is a guess');
      expect(df.get('Stage', 0), 'Made', 'our own edit');
      expect(df.get('Note', 2), 'A edited benzene', 'a molecule one row holds takes the other session\'s edit');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});
});

category('Campaign view: saves', () => {
  test('Modify Stages saves a layout without the progress tracker\'s Tile Viewer', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [app] = await openShared(sb, 1);
      const session = (app as unknown as SessionInternals)._session!;
      session.openTracker().getButton('Modify Stages')!.click();
      await waitUntil(() => u2Button('modify-stages', 'OK') !== undefined, 'Modify Stages did not open', 10000);
      u2Button('modify-stages', 'OK')!.click();
      await waitUntil(() => app.campaign?.version === 2 && !app.isSaving, 'Modify Stages did not save', 30000);
      const layout = (await sb.readCampaign('ZZB3-1')).layout;
      expect(typeof layout, 'string', 'a layout is saved');
      expect(layout.includes(`"${DG.VIEWER.TILE_VIEWER}"`), false, 'the saved layout has no Tile Viewer');
      expect(Array.from(app.designView!.viewers).some((v) => v.type === DG.VIEWER.TILE_VIEWER), true,
        'the tracker shows its lanes again');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a campaign whose layout has no grid opens with one, and its edits still save', async () => {
    const sb = await Sandbox.create('save');
    try {
      const layout = layoutOf(designCsv, 'Note', (tv) => tv.grid.close());
      expect(layout.includes('"Grid"'), false, 'the layout has no grid');
      const [app] = await openShared(sb, 1, {layout});
      const view = app.designView!;
      expect(view.grid != null, true, 'the view has a grid');
      view.grid.cell('Stage', 0).setValue('Tested', true);
      await waitUntil(() => app.campaign?.version === 2 && !app.isSaving, 'the Stage edit was not saved', 30000);
      expect((await serverRows(sb)).get('CCO')!['Stage'], 'Tested', 'the edit is in the saved table');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('Calculate additional properties: OK keeps the campaign\'s compute as the editor reads it, and saves', async () => {
    const sb = await Sandbox.create('save');
    try {
      const compute = {descriptors: {enabled: true, args: ['MolWt']}, functions: [], scripts: [], queries: []};
      const t = {...noComputeTemplate('ZZB3 Shared', 'ZZB3'), compute};
      await sb.writeTemplate(t);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {version: 1}), designCsv);
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      app.editCompute();
      await waitUntil(() => u2Button('compute', 'OK')?.disabled === false, 'the compute dialog never let OK go', 30000);
      u2Button('compute', 'OK')!.click();
      await waitUntil(() => app.campaign?.version === 2 && !app.isSaving, 'the campaign was not saved', 30000);
      expectSame((await sb.readCampaign('ZZB3-1')).template.compute, compute, 'the stored compute');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('PeptiHit: a column only the layout binds is not saved while it is empty', async () => {
    const sb = await Sandbox.create('save');
    try {
      const csv = 'Helm,Molecule,Stage,V-iD\nPEPTIDE1{A.G.L.V.K}$$$$V2.0,CC(N)C(=O)O,Design,\n';
      const t = {...noComputeTemplate('ZZB3 Peptides', 'ZZB3'), stages: ['Design', 'Synthesis'], layoutViewState:
        layoutOf('Helm,Molecule,Stage,V-iD,MolWt\nPEPTIDE1{A.G.L.V.K}$$$$V2.0,CC(N)C(=O)O,Design,,71.08\n', 'MolWt')};
      await sb.writeTemplate(t);
      const campaign = {name: 'ZZB3-1', templateName: t.name, status: 'In Progress', createDate: '2024/9/4',
        campaignFields: {Chemist: 'Davit'}, columnSemTypes: {'Helm': 'Macromolecule', 'Molecule': 'Molecule'},
        columnTypes: {'Helm': 'string', 'Molecule': 'string', 'Stage': 'string', 'V-iD': 'string'}, rowCount: 1,
        filteredRowCount: 1, savePath: sb.savePath('ZZB3-1'), version: 1, permissions: {view: [allUsers], edit: [allUsers]},
        template: t};
      await sb.writeCampaign(campaign, csv);
      const app = await sb.peptiHit();
      await sb.openDesign(app, 'ZZB3-1');
      expect(app.dataFrame!.columns.contains('MolWt'), true, 'the layout\'s scatter plot has its column, empty');
      await app.saveCampaign(false);
      const df = await grok.dapi.files.readCsv(sb.savePath('ZZB3-1'));
      expectSame(df.columns.names(), ['Helm', 'Molecule', 'Stage', 'V-iD'], 'the saved table\'s columns');
      const after = await sb.readCampaign('ZZB3-1');
      expect('MolWt' in after.columnTypes || 'MolWt' in after.columnSemTypes || 'MolWt' in after.columnTags, false,
        'campaign.json does not name it');
      const again = await sb.peptiHit();
      await sb.openDesign(again, 'ZZB3-1');
      expect(again.dataFrame!.columns.contains('MolWt'), true, 'reopened, the saved layout\'s column is back, empty');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 180000});

  test('the Stage cells offer the template\'s stages: set on open over a stored copy, changed by Modify Stages, never saved', async () => {
    const sb = await Sandbox.create('save');
    try {
      const [app] = await openShared(sb, 1, {columnTags: {Stage: {[DG.TAGS.CHOICES]: '["Old"]'}}});
      const choices = (a: HitDesignApp) => a.dataFrame!.col('Stage')!.getTag(DG.TAGS.CHOICES);
      expect(choices(app), JSON.stringify(['Idea', 'Made', 'Tested']), 'the template\'s stages, not the stored copy');
      await app.setStages(['Idea', 'Made', 'Tested', 'Shipped']);
      expect(choices(app), JSON.stringify(['Idea', 'Made', 'Tested', 'Shipped']), 'Modify Stages changes them');
      const saved = await sb.readCampaign('ZZB3-1');
      expect(saved.version, 2, 'Modify Stages saved');
      expect(DG.TAGS.CHOICES in (saved.columnTags.Stage ?? {}), false, 'campaign.json keeps no copy');
      expect(saved.layout.includes(DG.TAGS.CHOICES), false, 'nor does its layout');
      const again = await sb.hitDesign();
      await sb.openDesign(again, 'ZZB3-1');
      expect(choices(again), JSON.stringify(['Idea', 'Made', 'Tested', 'Shipped']), 'reopened, from the saved template');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a reopened number column shows as the grid shows numbers, unless the campaign chose its format', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = noComputeTemplate('ZZB3 Shared', 'ZZB3');
      await sb.writeTemplate(t);
      const readerFormat = `0.${'0'.repeat(13)}`;
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {version: 1, rowCount: 1,
        columnTypes: {...designColumns.columnTypes, MW: 'double', Score: 'double', Old: 'double'},
        columnTags: {Score: {format: '0.000'}, Old: {'format': readerFormat, '.source-precision': '13'}}}),
      'Molecule,Stage,V-iD,Note,MW,Score,Old\nCCO,Idea,,first,446.9100036621094,1.5,4.2755999565124\n');
      const app = await sb.hitDesign();
      await sb.openDesign(app, 'ZZB3-1');
      const df = app.dataFrame!;
      expect(df.col('MW')!.getTag(DG.TAGS.FORMAT), null, 'MW carries no format the CSV reader made up');
      expect(df.col('MW')!.getString(0), '446.91', 'MW shows as it did when it was calculated');
      expect(df.col('Old')!.getTag(DG.TAGS.FORMAT), null, 'nor one an earlier save kept');
      expect(df.col('Score')!.getTag(DG.TAGS.FORMAT), '0.000', 'the campaign\'s own format stays');
      await app.saveCampaign(false);
      const tags = (await sb.readCampaign('ZZB3-1')).columnTags;
      expect(DG.TAGS.FORMAT in tags.MW || DG.TAGS.SOURCE_PRECISION in tags.MW, false, 'the save keeps no reader format');
      expect(DG.TAGS.FORMAT in tags.Old, false, 'nor the one an earlier save kept');
      expect(tags.Score.format, '0.000', 'and keeps the campaign\'s own');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a new campaign\'s Molecule column is molblock: a drawn molecule is registered, calculated, saved and reopened as drawn', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = withMolWeight(noComputeTemplate('ZZB3 Molblock', 'ZZB3'));
      await sb.writeTemplate(t);
      const app = await sb.hitDesign();
      const id = await startDesignCampaign(app, t, {Chemist: 'Ann'}, undefined, drawnEthanol);
      expect(units(app), DG.UNITS.Molecule.MOLBLOCK, 'the new campaign\'s units');
      expect((await sb.readCampaign(id)).columnTags.Molecule.units, DG.UNITS.Molecule.MOLBLOCK, 'the units campaign.json keeps');
      expect((await grok.dapi.files.readCsv(sb.savePath(id))).get('Molecule', 0), drawnEthanol, 'the molblock it started with, saved as it is');
      app.designView!.grid.cell('Molecule', 0).setValue(drawnPropylamine, true);
      await waitUntil(() => app.campaign?.version === 2 && !app.isSaving, 'the drawn molecule was not saved', 60000);
      const vid = app.dataFrame!.get('V-iD', 0);
      expect(/^V\d{6}$/.test(vid), true, `the V-iD registered: ${vid}`);
      expect(Math.round(app.dataFrame!.get('MW', 0)), 59, 'the MW calculated');
      const saved = await grok.dapi.files.readCsv(sb.savePath(id));
      expect(saved.get('Molecule', 0), drawnPropylamine, 'the molblock saved as drawn');
      expect(saved.get('V-iD', 0), vid, 'the V-iD saved');
      const again = await sb.hitDesign();
      await sb.openDesign(again, id);
      expect(units(again), DG.UNITS.Molecule.MOLBLOCK, 'the units reopened');
      expect(again.dataFrame!.get('Molecule', 0), drawnPropylamine, 'the molecule reopened as drawn');
      expect(again.dataFrame!.get('V-iD', 0), vid, 'the V-iD reopened');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('a new campaign saved before anything is drawn reopens molblock; an empty sketch is neither registered nor calculated', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = withMolWeight(noComputeTemplate('ZZB3 Molblock', 'ZZB3'));
      await sb.writeTemplate(t);
      const app = await sb.hitDesign();
      const id = await startDesignCampaign(app, t, {Chemist: 'Ann'}, undefined, '');
      await closeSketcher();
      expect((await sb.readCampaign(id)).columnTags.Molecule.units, DG.UNITS.Molecule.MOLBLOCK, 'the units campaign.json keeps');
      const again = await sb.hitDesign();
      await sb.openDesign(again, id);
      await closeSketcher();
      expect(units(again), DG.UNITS.Molecule.MOLBLOCK, 'the empty column reopened');
      again.dataFrame!.set('Molecule', 0, emptySketch);
      await again.registerRow(0);
      await again.performSingleCellCalculations(0, emptySketch);
      expect(again.dataFrame!.col('V-iD')!.isNone(0), true, 'no V-iD for an empty sketch');
      expect(again.dataFrame!.columns.contains('MW'), false, 'no MW calculated for an empty sketch');
      expect(again.isSaving || again.campaign?.version !== 1, false, 'nothing saved');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('an existing SMILES campaign keeps its SMILES and its units through an open and a save', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = noComputeTemplate('ZZB3 Shared', 'ZZB3');
      await sb.writeTemplate(t);
      const smilesTags = {'quality': 'Molecule', 'units': DG.UNITS.Molecule.SMILES, 'cell.renderer': 'Molecule'};
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {version: 1, columnTags: {Molecule: smilesTags}}), designCsv);
      // the shape of a campaign saved before campaign.json kept column tags (DMT1-1): its units stay unset
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-2', t, {version: 1}), designCsv);
      for (const [id, expected] of [['ZZB3-1', DG.UNITS.Molecule.SMILES], ['ZZB3-2', null]]) {
        const app = await sb.hitDesign();
        await sb.openDesign(app, id!);
        expect(units(app) ?? null, expected, `${id} units on open`);
        await app.saveCampaign(false);
        expect((await sb.readCampaign(id!)).columnTags.Molecule.units ?? null, expected, `${id} units saved`);
        expectSame((await grok.dapi.files.readCsv(sb.savePath(id!))).col('Molecule')!.toList(), ['CCO', 'c1ccccc1'], `${id} molecules saved`);
      }
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('a SMILES file merged into a molblock campaign matches its molecules by canonical SMILES', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = noComputeTemplate('ZZB3 Molblock', 'ZZB3');
      await sb.writeTemplate(t);
      const app = await sb.hitDesign();
      await startDesignCampaign(app, t, {Chemist: 'Ann'}, undefined, drawnEthanol);
      await mergeIntoCampaign(app, DG.DataFrame.fromCsv('Molecule,Note\nOCC,merged\nCCN,new\n'), {mode: 'smiles',
        molColName: 'Molecule', addNewRows: true, runComputeOnNewRows: false, clashStrategy: 'overwrite'});
      const df = app.dataFrame!;
      expect(df.rowCount, 2, 'the matched row kept, the new molecule added');
      expect(df.get('Molecule', 0), drawnEthanol, 'the drawn molecule kept as drawn');
      expect(df.get('Note', 0), 'merged', 'OCC matched the drawn ethanol');
      expect(df.get('Molecule', 1), 'CCN', 'the new molecule');
      expect(/^V\d{6}$/.test(df.get('V-iD', 1)), true, `the new molecule's V-iD: ${df.get('V-iD', 1)}`);
      expect(units(app), DG.UNITS.Molecule.MOLBLOCK, 'the units after the merge');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 150000});

  test('PeptiHit: the molecule a HELM converts to is a molblock, saved and reopened as it is', async () => {
    const sb = await Sandbox.create('save');
    try {
      const t = {...noComputeTemplate('ZZB3 Peptides', 'ZZB3'), stages: ['Design', 'Synthesis']};
      await sb.writeTemplate(t);
      const helm = 'PEPTIDE1{A.G.L}$$$$V2.0';
      const app = await sb.peptiHit();
      const id = await startDesignCampaign(app, t, {Chemist: 'Ann'}, undefined, helm);
      expect(units(app), DG.UNITS.Molecule.MOLBLOCK, 'the new campaign\'s units');
      await app.performSingleCellCalculations(0, helm);
      await waitUntil(() => app.campaign?.version === 2 && !app.isSaving, 'the converted molecule was not saved', 60000);
      const mol = app.dataFrame!.get('Molecule', 0);
      expect(_package.convertToSmiles(mol) !== '' && mol.includes('M  END'), true, `a molblock of the HELM: ${mol}`);
      expect((await grok.dapi.files.readCsv(sb.savePath(id))).get('Molecule', 0), mol, 'the molblock saved as it is');
      const again = await sb.peptiHit();
      await sb.openDesign(again, id);
      expect(units(again), DG.UNITS.Molecule.MOLBLOCK, 'the units reopened');
      expect(again.dataFrame!.get('Molecule', 0), mol, 'the molecule reopened as it is');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 180000});
});

/** Each card of the open progress tracker: the V-iD it shows, its structure's picture, its element. */
function trackerCards(): {vid: string, picture: string, card: Element}[] {
  return Array.from(document.querySelectorAll('[data-u2-name="progress-tracker"] .d4-tile-viewer-form')).map((card) => ({
    vid: Array.from(card.querySelectorAll('input')).map((i) => i.value).find((v) => /^ZV\d$/.test(v)) ?? '',
    picture: card.querySelector('canvas')?.toDataURL() ?? '',
    card,
  }));
}

const picturesByVid = (cards: {vid: string, picture: string}[]) =>
  cards.map((c) => `${c.vid}:${c.picture.length}:${c.picture.slice(-48)}`).sort().join(' | ');

category('Campaign view: tracker', () => {
  test('a card the tracker refills for another row shows that row\'s structure', async () => {
    const sb = await Sandbox.create('tracker');
    try {
      const t = noComputeTemplate('ZZB3 Shared', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {version: 1, rowCount: 3}),
        'Molecule,Stage,V-iD,Note\nCCO,Idea,ZV1,a\nc1ccccc1,Idea,ZV2,b\nCCN,Made,ZV3,c\n');
      const app = await sb.hitDesign();
      const view = await sb.openDesign(app, 'ZZB3-1');
      (app as unknown as SessionInternals)._session!.openTracker();
      let last = '';
      await waitUntil(() => {
        const cards = trackerCards();
        const now = picturesByVid(cards);
        const settled = cards.length === 3 && cards.every((c) => c.vid !== '' && c.picture !== '') && now === last;
        last = now;
        return settled;
      }, 'the tracker\'s three cards were not drawn', 20000);
      const before = last;
      const hosts = trackerCards().map((c) => c.card);
      // the lanes laid out again: the cards go to the spare ones and come back refilled for other rows
      view.dataFrame.rows.requestFilter();
      await waitUntil(() => trackerCards().every((c) => !hosts.includes(c.card)), 'the lanes were not laid out again', 10000);
      await waitUntil(() => picturesByVid(trackerCards()) === before,
        `a refilled card shows another row's structure: ${picturesByVid(trackerCards())} vs ${before}`, 10000);
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('closing the campaign view with the tracker open keeps its moves, and saves them', async () => {
    const sb = await Sandbox.create('tracker');
    try {
      const [app] = await openShared(sb, 1);
      (app as unknown as SessionInternals)._session!.openTracker();
      await waitUntil(() => u2Button('progress-tracker', 'CLOSE') !== undefined, 'the tracker did not open', 10000);
      expect(u2Button('progress-tracker', 'CANCEL') === undefined && u2Button('progress-tracker', 'OK') === undefined, true,
        'CLOSE is the tracker\'s one action');
      const df = app.dataFrame!;
      df.set('Stage', 0, 'Tested');
      app.designView!.close();
      await waitUntil(() => app.campaign?.version === 2 && !app.isSaving, 'the move was not saved', 30000);
      expect(df.get('Stage', 0), 'Tested', 'the move stays in the table');
      expect((await serverRows(sb)).get('CCO')!['Stage'], 'Tested', 'the move is in the saved table');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});
});
