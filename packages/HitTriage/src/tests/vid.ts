/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package';
import * as api from '../package-api';
import {deobfuscateSmiles, obfuscateSmiles, registerMol, registerMolsBatch} from '../app/utils/molreg';
import {countQueryCalls, expectSame, hitdesignSql, Sandbox, sqlText} from './fixtures';

/** Long primary amines: valid, distinct, and found in no real campaign. */
const chainAmines = (count: number) => Array.from({length: count}, (_, i) => `${'C'.repeat(20 + i)}N`);

category('V-iD', () => {
  test('the SMILES obfuscation round-trips and hides the SMILES', async () => {
    for (const s of ['CCO', 'c1ccccc1', 'C[C@H](N)C(=O)O', '[Na+].[Cl-]', 'O=C(O)c1ccccc1OC(C)=O', 'C'.repeat(300)]) {
      const o = obfuscateSmiles(s);
      expect(deobfuscateSmiles(o), s, `round trip of ${s}`);
      expect(o !== s && /^[A-Za-z0-9+/]+=*$/.test(o), true, `${s} is stored as base64, not as itself`);
      expect(obfuscateSmiles(s), o, `the obfuscation of ${s} is stable`);
    }
  });

  test('registerMolsBatch registers in groups of 50 and keeps the order', async () => {
    const sb = await Sandbox.create('vid');
    try {
      const generated = chainAmines(120).map((m) => sb.ownMolecule(m));
      const mols = [...generated, 'CCO', 'OCC'];
      const {result: vids, calls} = await countQueryCalls('AddMolecules', 'smiles',
        () => registerMolsBatch(mols, 'ZZB3-1', sb.appName));
      expectSame(calls, [50, 50, 21], 'molecules per AddMolecules call (121 distinct canonical SMILES)');
      expect(mols.every((m) => /^V\d{6,}$/.test(vids.get(m) ?? '')), true, 'every molecule has a V-iD');
      expect(vids.get('OCC'), vids.get('CCO'), 'two notations of one molecule share a V-iD');
      expect(new Set(vids.values()).size, 121, 'distinct V-iDs');

      const list = Array.from(new Set(vids.values()));
      const stored = await hitdesignSql(`SELECT vid, mh_string FROM hitdesign.vid_dictionary WHERE vid IN (${list.map(sqlText).join(',')})`);
      const byVid = new Map<string, string>();
      for (let i = 0; i < stored.rowCount; i++)
        byVid.set(stored.get('vid', i), deobfuscateSmiles(stored.get('mh_string', i)));
      const wrong = mols.filter((m) => byVid.get(vids.get(m)!) !== _package.convertToSmiles(m));
      expect(wrong.length, 0, `molecules whose V-iD holds another molecule: ${wrong.slice(0, 3).join(', ')}`);

      const first = vids.get(mols[0])!;
      const byQuery = await api.queries.getMoleculeByVid(first);
      expect(deobfuscateSmiles(byQuery.get('mh_string', 0)), _package.convertToSmiles(mols[0]), 'GetMoleculeByVid');
      const single = await countQueryCalls('AddMolecule', null, () => registerMol(mols[0], 'ZZB3-2', sb.appName));
      expect(single.result, first, 'registering a known molecule returns its V-iD');
      expect(single.calls.length, 1, 'AddMolecule calls');

      const campaigns = await api.queries.getCampaignsByVid(first, sb.appName);
      const pairs = Array.from({length: campaigns.rowCount}, (_, i) =>
        `${campaigns.get('app_name', i)}|${campaigns.get('campaign_id', i)}|${campaigns.get('created_by', i)}`).sort();
      const login = DG.User.current().login;
      expectSame(pairs, [`${sb.appName}|ZZB3-1|${login}`, `${sb.appName}|ZZB3-2|${login}`], 'GetCampaignsByVid');
      expect((await api.queries.getCampaignsByVid(first, 'ZZB3 nowhere')).rowCount, 0, 'GetCampaignsByVid of another app');
      expect((await api.queries.getCampaignsByVid(first, null)).rowCount >= 2, true, 'GetCampaignsByVid of every app');
      const count = await hitdesignSql(`SELECT COUNT(*) AS n FROM hitdesign.campaign_vids WHERE app_name = ${sqlText(sb.appName)}`);
      expect(count.get('n', 0), 122, 'campaign_vids rows: 121 in ZZB3-1, 1 in ZZB3-2');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('empty and malformed molecules get no V-iD and send no query', async () => {
    const sb = await Sandbox.create('vid');
    try {
      const {result: vids, calls} = await countQueryCalls('AddMolecules', 'smiles',
        () => registerMolsBatch(['', 'not a molecule'], 'ZZB3-1', sb.appName));
      expect(vids.size, 0, 'V-iDs');
      expect(calls.length, 0, 'AddMolecules calls');
      const single = await countQueryCalls('AddMolecule', null, () => registerMol('', 'ZZB3-1', sb.appName));
      expect(single.result, '', 'registerMol of an empty cell');
      expect(single.calls.length, 0, 'AddMolecule calls');
    } finally {
      await sb.dispose();
    }
  });
});
