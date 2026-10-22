import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, category, delay, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {domains, DomainApp, DomainPick, DomainSearch} from '@datagrok-libraries/ems/src/dg';
import {Control, Filters, Rows} from '@datagrok-libraries/u2';
import {moltrackEmsDb} from '../generated/db';
import {Chemistry, MatchRule} from '../ems/chemistry';
import {MolTrackEms, VIEWS} from '../ems/app';
import {CompoundsApp} from '../ems/compounds-app';
import {MolTrackTables} from '../ems/tables';
import {Registration} from '../ems/registration';
import {RegisterDialog} from '../ems/register-dialog';
import {DemoData} from '../ems/demo-data';
import {StructureSearch} from '../ems/structure-search';

// Expected values are RDKit 2026.03.6's (rdkit.Chem.RegistrationHash, MolToInchiKey, CalcMolFormula)
// for the same inputs — what Crux reproduces.
category('MolTrack EMS: chemistry', () => {
  test('aspirin: standardized structure, hash, InChIKey and formula', async () => {
    const id = await Chemistry.identify('CC(=O)Oc1ccccc1C(=O)O', 'ALL_LAYERS');
    expect(id.structure, 'CC(=O)Oc1ccccc1C(=O)O');
    expect(id.hash_mol, '60ec019cc6977893a0e2fe68e05fb6f163dc98a3');
    expect(id.inchikey, 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N');
    expect(id.formula, 'C9H8O4');
  });

  test('a salt is the same compound as its parent; InChI and formula are of the structure as given', async () => {
    const salt = await Chemistry.identify('CC(=O)Oc1ccccc1C(=O)[O-].[Na+]', 'ALL_LAYERS');
    expect(salt.hash_mol, '60ec019cc6977893a0e2fe68e05fb6f163dc98a3');
    expect(salt.inchikey, 'JZLOKWGVGHYBKD-UHFFFAOYSA-M');
    expect(salt.formula, 'C9H7NaO4');
  });

  test('the matching rule decides what is the same: stereo and tautomers', async () => {
    const hash = async (smiles: string, rule: MatchRule) =>
      (await Chemistry.identify(smiles, rule)).hash_mol;
    expect(await hash('C[C@H](N)C(=O)O', 'ALL_LAYERS') === await hash('C[C@@H](N)C(=O)O', 'ALL_LAYERS'), false);
    expect(await hash('C[C@H](N)C(=O)O', 'STEREO_INSENSITIVE_LAYERS'),
      await hash('C[C@@H](N)C(=O)O', 'STEREO_INSENSITIVE_LAYERS'));
    expect(await hash('O=c1cccc[nH]1', 'ALL_LAYERS') === await hash('Oc1ccccn1', 'ALL_LAYERS'), false);
    expect(await hash('O=c1cccc[nH]1', 'TAUTOMER_INSENSITIVE_LAYERS'),
      await hash('Oc1ccccn1', 'TAUTOMER_INSENSITIVE_LAYERS'));
  });

  test('substructure, superstructure and similarity', async () => {
    const set = ['CC(=O)Oc1ccccc1C(=O)O', 'c1ccccc1', 'CCO', 'not a molecule'];
    expectArray(await Chemistry.substructure(set, 'c1ccccc1'), [0, 1]);
    // ethanol is the C-C-O of aspirin's acetyl ester
    expectArray(await Chemistry.superstructure(set, 'CC(=O)Oc1ccccc1C(=O)O'), [0, 1, 2]);
    const scores = await Chemistry.similarity(set, 'CC(=O)Oc1ccccc1C(=O)O');
    expect(scores[0], 1);
    expect(scores[2] < 0.2, true);
    expect(Number.isNaN(scores[3]), true);
  });

  test('a structure Crux cannot read is refused by name', async () => {
    let message = '';
    try {
      await Chemistry.identify('C1CC', 'ALL_LAYERS');
    } catch (e) {
      message = (e as Error).message;
    }
    expect(message, 'Cannot read the structure');
  });
});

// Real rows in moltrack_ems, soft-deleted in after(); each test registers its own unusual structures.
category('MolTrack EMS: registration', () => {
  const tag = `T${Date.now()}`;
  const chain = (n: number, end = 'C(=O)O') => `${'C'.repeat(n)}${end}`;
  const created = {compounds: new Set<string>(), batches: new Set<string>()};

  /** The compounds a run registered — not the ones it found registered already, which may be the demo's. */
  const remember = async (report: DG.DataFrame) => {
    const ids = report.col(Registration.COMPOUND_ID)!;
    for (let i = 0; i < report.rowCount; i++) {
      if (report.get(Registration.STATUS, i) !== 'success')
        continue;
      const row = await moltrackEmsDb.compounds.first({filter: `molregno = ${Number(DomainPick.digits(ids.get(i)))}`});
      if (row !== null)
        created.compounds.add(row.id);
    }
  };
  const register = async (rows: Record<string, string>[], rejectAll = false) => {
    const report = await MolTrackEms.register(DG.DataFrame.fromObjects(rows)!, rejectAll);
    await remember(report);
    return report;
  };
  const statuses = (report: DG.DataFrame) => report.col(Registration.STATUS)!.toList();
  const messages = (report: DG.DataFrame) => report.col(Registration.MESSAGE)!.toList();
  // an empty cell is a hole in the list
  const grokIds = (report: DG.DataFrame) =>
    Array.from(report.col(Registration.COMPOUND_ID)!.toList(), (id) => id ?? '');
  const compoundOf = async (grokId: string) =>
    (await moltrackEmsDb.compounds.first({filter: `molregno = ${Number(DomainPick.digits(grokId))}`}))!;

  after(async () => {
    const failed: string[] = [];
    for (const id of created.batches)
      await moltrackEmsDb.compoundBatches.delete(id).catch((e) => failed.push(`batch ${id}: ${e}`));
    for (const id of created.compounds)
      await moltrackEmsDb.compounds.delete(id).catch((e) => failed.push(`compound ${id}: ${e}`));
    if (failed.length > 0)
      throw new Error(`Cleanup failed: ${failed.join('; ')}`);
  });

  test('bulk: a new structure registers once, its salt and a registered one are reported with their Grok IDs',
    async () => {
      const acid = await register([{'smiles': chain(25), 'Common Name': `${tag} acid`}]);
      expectArray(statuses(acid), ['success']);
      const acidId = grokIds(acid)[0];
      expect(/^DG-\d{6}$/.test(acidId), true, `corporate compound id ${acidId}`);
      const report = await register([
        {'smiles': chain(26), 'Common Name': `${tag} new`},
        {'smiles': chain(25, 'C(=O)[O-].[Na+]'), 'Common Name': `${tag} registered salt`},
        {'smiles': 'C1CC', 'Common Name': `${tag} broken`},
        {'smiles': chain(26, 'C(=O)[O-].[Na+]'), 'Common Name': `${tag} salt of row 1`},
      ]);
      expectArray(statuses(report), ['success', 'existing', 'failed', 'existing']);
      const ids = grokIds(report);
      expectArray(ids.slice(1), [acidId, '', ids[0]]);
      expectArray(messages(report).slice(1), [`This structure is already registered as ${acidId}`,
        'Cannot read the structure', 'Same compound as row 1']);
      const hash = (await Chemistry.identify(chain(26), MolTrackTables.matchRule())).hash_mol;
      expect((await moltrackEmsDb.compounds.query({filter: MolTrackTables.within('hash_mol', [hash])})).length, 1);
      const compound = await compoundOf(ids[0]);
      expect(compound.common_name, `${tag} new`, 'the first row of a structure carries the compound');
      expect(compound.structure, chain(26));
      expect((await compoundOf(acidId)).common_name, `${tag} acid`, 'a registered compound is not written');
    });

  test('a bad row fails alone; with reject all nothing is written', async () => {
    const rows = [{'smiles': 'C1CC', 'Common Name': `${tag}-3`}, {'smiles': chain(27, 'N'), 'Common Name': `${tag}-4`}];
    const strict = await register(rows, true);
    expectArray(statuses(strict), ['failed', 'not_processed']);
    expect(await moltrackEmsDb.compounds.count({filter: `common_name = "${tag}-4"`}), 0, 'nothing was written');
    const lenient = await register(rows);
    expectArray(statuses(lenient), ['failed', 'success']);
    expect(messages(lenient)[0], 'Cannot read the structure');
  });

  test('reject all is one transaction: a compound the server refuses leaves no other compound behind', async () => {
    const name = `${tag} atomic`;
    const report = await register([
      {'smiles': chain(28), 'Common Name': name, 'MolLogP': ''},
      {'smiles': chain(29), 'Common Name': `${tag} refused`, 'MolLogP': 'high'},
    ], true);
    expectArray(statuses(report), ['not_processed', 'failed']);
    expect(/MolLogP/.test(messages(report)[1]), true, `the refusal names the column: ${messages(report)[1]}`);
    expect(await moltrackEmsDb.compounds.count({filter: `common_name = "${name}"`}), 0);
  });

  test('a row with no structure fails alone', async () => {
    const report = await register([{'smiles': ''}, {'smiles': chain(30)}]);
    expectArray(statuses(report), ['failed', 'success']);
    expect(messages(report)[0], 'Structure is required');
  });

  test('every rule\'s hash is stored, so a stricter rule today and a looser one later still match', async () => {
    const l = `C[C@H](N)${'C'.repeat(20)}C(=O)O`;
    const d = `C[C@@H](N)${'C'.repeat(20)}C(=O)O`;
    const report = await register([{smiles: l}]);
    expectArray(statuses(report), ['success']);
    const looser = await Chemistry.identify(d, 'STEREO_INSENSITIVE_LAYERS');
    const found = await MolTrackTables.registered('hash_stereo_insensitive', [looser.hash_mol]);
    expect(found.size, 1, 'the enantiomer is the same compound once stereo is ignored');
    const stricter = await Chemistry.identify(d, 'ALL_LAYERS');
    expect((await MolTrackTables.registered('hash_all_layers', [stricter.hash_mol])).size, 0);
  });

  test('a corporate compound ID opens its context panel', async () => {
    const report = await register([{'smiles': chain(31), 'Common Name': `${tag} panel`}]);
    const panel = await MolTrackEms.compoundPanel(report.col(Registration.COMPOUND_ID)!.get(0));
    expect(panel.root.textContent!.includes(`${tag} panel`), true);
    expect((await MolTrackEms.compoundPanel('DG-999999')).root.textContent,
      'DG-999999 is not registered in MolTrack EMS');
  });

  test('the structure conditions find a registered compound, and a list filtered by one holds it', async () => {
    await register([{smiles: chain(32)}]);
    const id = (await MolTrackTables.registered('hash_mol',
      [(await Chemistry.identify(chain(32), MolTrackTables.matchRule())).hash_mol])).values().next().value!.id;
    const matching = (operator: string, value: string, options?: Record<string, unknown>) =>
      StructureSearch.matching({...Filters.cond('structure', operator, value), ...(options ? {options} : {})});
    expectArray(await matching(StructureSearch.EXACT, chain(32, 'C(=O)[O-].[Na+]')), [id]);
    expect((await matching(StructureSearch.SUBSTRUCTURE, chain(32))).includes(id), true);
    expect((await matching(StructureSearch.SIMILAR, chain(32), {threshold: 0.9})).includes(id), true);
    expect((await matching(StructureSearch.SUPERSTRUCTURE, 'CC(=O)O')).includes(id), false);
    const source = (await MolTrackTables.compounds()).source(
      {query: `structure ${StructureSearch.SUBSTRUCTURE} "${chain(32)}"`});
    try {
      await awaitCheck(() => source.state.value === 'ready', 'the filtered list never loaded', 30000);
      expectArray(source.rows.items.value.map((r) => r.id), [id]);
    } finally {
      source.dispose();
    }
  });

  test('the preview says what each row becomes and writes nothing', async () => {
    await register([{smiles: chain(33)}]);
    const df = DG.DataFrame.fromObjects([{smiles: chain(34)}, {smiles: chain(33, 'C(=O)[O-].[Na+]')},
      {smiles: 'C1CC('}, {smiles: chain(34)}])!;
    const preview = await Registration.preview(df, Registration.autoMap(df, await MolTrackTables.compounds()));
    expectArray(preview.map((p) => p.verdict), ['new', 'existing', 'error', 'new']);
    expect(preview[3].message, 'Same compound as row 1');
    expect(preview[1].message.includes('DG-'), true);
    expect((await MolTrackTables.registered('hash_mol',
      [(await Chemistry.identify(chain(34), MolTrackTables.matchRule())).hash_mol])).size, 0);
  });

  test('a new compound drawn as a registered structure is refused on its Structure, naming the compound', async () => {
    const report = await register([{smiles: chain(36)}]);
    const corporateId = report.col(Registration.COMPOUND_ID)!.get(0);
    const view = await MolTrackEms.view('compound');
    grok.shell.addView(view);
    const app = DomainApp.of(view)!;
    try {
      await app.create();
      await awaitCheck(() => app.entitySource.value?.currentRow.value != null, 'no new compound', 10000);
      app.entitySource.value!.currentRow.value!.original_structure = chain(36, 'C(=O)[O-].[Na+]');
      const message = MolTrackTables.alreadyRegistered(corporateId);
      await awaitCheck(() => app.form.value?.input('original_structure')?.validity.value === message,
        'the Structure input never said the structure is registered', 15000);
      const section = view.root.querySelector<HTMLElement>('.moltrack-ems-compound-section')!;
      await awaitCheck(() => !section.hidden && section.textContent!.includes(corporateId),
        'the page never showed the registered compound', 10000);
      expect(await app.session.save(), false);
      const problem = app.form.value!.problem.value ?? '';
      expect(problem.includes(message) && !/hash/i.test(problem), true, `the refusal reads "${problem}"`);
    } finally {
      app.session.discard();
      view.close();
    }
  });

  test('a new batch on a picked compound saves on it, the compound card shown, and stays on it', async () => {
    const report = await register([{'smiles': chain(35, 'O'), 'Common Name': `${tag} picked`}]);
    const corporateId = grokIds(report)[0];
    const compound = await compoundOf(corporateId);
    const view = await MolTrackEms.view('compound_batch');
    grok.shell.addView(view);
    const app = DomainApp.of(view)!;
    const section = () => view.root.querySelector<HTMLElement>('.moltrack-ems-compound-section')!;
    try {
      await app.create();
      await awaitCheck(() => app.entitySource.value?.currentRow.value != null && section() !== null, 'no new batch',
        10000);
      expect(section().hidden, true, 'no compound is shown before one is picked');
      expect(app.form.value!.input('compound_id') !== undefined, true, 'a new batch picks its compound');
      const batch = app.entitySource.value!.currentRow.value!;
      batch.compound_id = compound.id;
      batch.epa_batch_id = `${tag}-picked`;
      await awaitCheck(() => !section().hidden && section().textContent!.includes(corporateId),
        'the page never showed the picked compound', 10000);
      expect(await app.session.save(), true);
      const saved = app.entitySource.value!.currentRow.value!.id;
      created.batches.add(saved);
      expect((await moltrackEmsDb.compoundBatches.get(saved)).compound_id, compound.id);
      await awaitCheck(() => app.form.value?.input('compound_id') === undefined && !section().hidden &&
        section().textContent!.includes(corporateId), 'the saved batch does not read as its compound', 10000);
      const other = await register([{smiles: chain(37)}]);
      let refusal = '';
      await moltrackEmsDb.compoundBatches.update(saved, {compound_id: (await compoundOf(grokIds(other)[0])).id})
        .catch((e) => refusal = String((e as Error).message ?? e));
      expect(/write-once/.test(refusal), true, `a batch was moved to another compound: ${refusal}`);
    } finally {
      app.session.discard();
      view.close();
    }
  });

  const amide = (n: number) => chain(n, 'C(=O)N');

  test('with reject all, a row left unregistered names the row that stopped the run', async () => {
    const report = await register([{'smiles': amide(30)}, {'smiles': 'C1CC'}], true);
    expectArray(statuses(report), ['not_processed', 'failed']);
    expectArray(messages(report), ['Not registered: row 2 has a problem', 'Cannot read the structure']);
  });

  test('the preview tells a structure the file repeats from a new compound', async () => {
    const df = DG.DataFrame.fromObjects([{smiles: amide(31)}, {smiles: amide(32)}, {smiles: amide(31)}])!;
    const preview = await Registration.preview(df, Registration.autoMap(df, await MolTrackTables.compounds()));
    expectArray(preview.map((p) => p.sameAs), [undefined, undefined, 0]);
  });

  /** The bulk wizard over `rows`, open in the shell as the table a registrar would pick. */
  const bulkWizard = async (rows: Record<string, unknown>[]) => {
    const df = DG.DataFrame.fromObjects(rows)!;
    df.name = `${tag}-bulk`;
    const table = grok.shell.addTableView(df);
    const running = RegisterDialog.open(df);
    await awaitCheck(() => document.querySelector('.moltrack-ems-register') !== null, 'the bulk wizard never opened',
      30000);
    return {running, table, wizard: document.querySelector<HTMLElement>('.moltrack-ems-register')!};
  };
  const press = (wizard: HTMLElement, label: string) => [...wizard.querySelectorAll('button')]
    .find((b) => b.textContent!.trim().toUpperCase() === label)!.click();

  test('bulk registration rejects all by default, and refuses two file columns mapped onto one field', async () => {
    // a frame's column names ignore case: "CAS" and "cas_" both match the compound's CAS
    const {running, table, wizard} = await bulkWizard([{'smiles': amide(33), 'CAS': '50-78-2', 'cas_': '50-78-2'}]);
    try {
      const onError = [...wizard.querySelectorAll('select')].find((s) => s.textContent!.includes('Reject all'))!;
      expect(onError.selectedOptions[0].textContent, 'Reject all');
      press(wizard, 'NEXT');
      const reason = () => wizard.querySelector('.u2-wizard-reason')?.textContent ?? '';
      await awaitCheck(() => reason().includes('Two columns of the file are mapped to CAS'),
        'the mapping of two columns onto one field was not refused', 10000);
      press(wizard, 'NEXT');
      expect(wizard.querySelector('.moltrack-ems-register-preview .u2-data-table') === null, true,
        'the preview opened over a mapping that claims one field twice');
    } finally {
      press(wizard, 'CANCEL');
      table.close();
    }
    expect(await running, null);
  });

  test('bulk preview and report name each row by its own columns, and the report opens the compounds registered',
    async () => {
      const {running, table, wizard} = await bulkWizard([
        {'smiles': amide(34), 'Common Name': `${tag} bulk one`},
        {'smiles': amide(35), 'Common Name': `${tag} bulk two`}]);
      const heads = (host: Element) => [...host.querySelectorAll('.u2-data-table-head')].map((h) => h.textContent);
      let view: DG.ViewBase | undefined;
      try {
        press(wizard, 'NEXT');
        press(wizard, 'NEXT');
        const preview = wizard.querySelector('.moltrack-ems-register-preview')!;
        await awaitCheck(() => preview.querySelectorAll('.u2-data-table-row').length === 2,
          'the preview never listed the rows', 30000);
        expectArray(heads(preview), ['#', 'Structure', 'Common Name', 'Becomes', '']);
        expect(preview.textContent!.includes(`${tag} bulk one`), true);
        expect(preview.textContent!.includes('2 rows: 2 new compounds'), true, preview.textContent!);
        press(wizard, 'REGISTER');
        const report = wizard.querySelector('.moltrack-ems-register-report')!;
        await awaitCheck(() => report.textContent!.includes('Registered 2 of 2 rows'), 'the rows were not registered',
          60000);
        expectArray(heads(report), ['#', 'Common Name', 'Result', 'Grok ID', '']);
        const lines = [...report.querySelectorAll('.u2-data-table-row')].map((r) => r.textContent!);
        expect(lines.length === 2 && lines.every((l) => /Registered/.test(l) && /DG-\d{6}/.test(l)), true,
          lines.join(' / '));
        [...report.querySelectorAll<HTMLElement>('a, .d4-link-action')]
          .find((a) => a.textContent === 'Show the 2 compounds registered')!.click();
        const registered = (await running)!;
        await remember(registered);
        const molregnos = grokIds(registered).map((id) => Number(DomainPick.digits(id)));
        const shown = () => [...DomainApp.live].find((a) => a.base === MolTrackEms.baseOf(VIEWS[0]) &&
          String(a.listSource.query.value).startsWith('molregno in'));
        await awaitCheck(() => shown() !== undefined, 'the compounds never opened', 30000);
        view = grok.shell.v;
        const app = shown()!;
        await awaitCheck(() => app.listSource.state.value === 'ready' &&
          app.listSource.rows.items.value.map((r) => r.molregno).sort().join() === molregnos.sort().join(),
        'the compounds list does not hold exactly the registered compounds', 30000);
      } finally {
        if (document.contains(wizard))
          press(wizard, 'CANCEL');
        view?.close();
        table.close();
      }
    }, {timeout: 120000});

  test('Batches import names each batch\'s compound by its Grok ID; a Grok ID no compound has fails alone',
    async () => {
      const report = await register([{'smiles': amide(36), 'Common Name': `${tag} imported`}]);
      const grokId = grokIds(report)[0];
      const df = DG.DataFrame.fromObjects([{'Grok ID': grokId, 'EPA Batch ID': `${tag}-imp-1`, 'Purity': '99'},
        {'Grok ID': 'DG-999999', 'EPA Batch ID': `${tag}-imp-2`, 'Purity': '98'}])!;
      const view = grok.shell.addView(await MolTrackEms.view('compound_batch'));
      const app = DomainApp.of(view)!;
      const named = (name: string) => Control.forElement(document.querySelector(`[data-u2-name="${name}"]`)) as
        unknown as {value: {value: unknown}} | undefined;
      const button = (label: string) => [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog button')]
        .find((b) => b.textContent === label)!;
      const text = () => document.querySelector('.u2-domain-import-report')?.textContent ?? '';
      try {
        expect(app.importRows === DomainApp.prototype.importRows, true, 'Batches import is the standard import');
        const running = domains.import(app.table, {source: df});
        await awaitCheck(() => named('allOrNothing') !== undefined, 'the import never opened', 15000);
        named('allOrNothing')!.value.value = false;
        button('NEXT').click();
        await awaitCheck(() => named('Grok ID') !== undefined, 'no mapping step', 10000);
        named('Grok ID')!.value.value = 'compound_id';
        await delay(100);
        button('NEXT').click();
        await delay(500);
        button('NEXT').click();
        await awaitCheck(() => text().includes('inserted'), `the import never reported: ${text()}`, 30000);
        expect(text().includes('1 inserted'), true, text());
        expect(text().includes('No match for "DG-999999"'), true, text());
        button('CLOSE').click();
        await running;
        const batches = await moltrackEmsDb.compoundBatches.query({filter: `epa_batch_id starts "${tag}-imp"`});
        for (const batch of batches)
          created.batches.add(batch.id);
        expectArray(batches.map((b) => b.epa_batch_id), [`${tag}-imp-1`]);
        expect(batches[0].compound_id, (await compoundOf(grokId)).id);
      } finally {
        if (document.querySelector('.u2-domain-import') !== null)
          button('CANCEL')?.click();
        view.close();
      }
    }, {timeout: 120000});

  test('every demo structure is the compound its row names, and every demo CAS is one the form takes', async () => {
    const demo = await DemoData.table('batches');
    const rows = Array.from({length: demo.rowCount}, (_, i) => i);
    const smilesOf = (name: string) =>
      demo.get('smiles', rows.find((i) => demo.get('Common Name', i).trim() === name)!);
    for (const [name, formula] of [['Bisphenol A', 'C15H16O2'], ['Cymoxanil', 'C7H10N4O3'],
      ['4,4\'-Sulfonylbis[2-(prop-2-en-1-yl)phenol]', 'C18H18O4S'], ['Bisphenol AF', 'C15H10F6O2']])
      expect((await Chemistry.identify(smilesOf(name), 'ALL_LAYERS')).formula, formula, name);
    const compounds = await MolTrackTables.compounds();
    const refused = rows.map((i) => demo.get('CAS', i) as string)
      .filter((cas) => compounds.validators.check('cas', cas, {} as never) !== null);
    expectArray(refused, []);
  });

  test('a stored date is read as its UTC day, as the lists read it', async () => {
    expectArray([DemoData.day(dayjs('2024-02-02T23:30:00Z')), DemoData.day(new Date('2024-02-03T00:30:00Z')),
      DemoData.day('2/3/2024'), DemoData.day('2024-02-03')], ['2024-02-02', '2024-02-03', '2024-02-03', '2024-02-03']);
  });

  test('the demo assay runs are dated at UTC midnight of the day their names give', async () => {
    const runs = await MolTrackTables.allPages((page) =>
      moltrackEmsDb.assayRuns.query({columns: ['name', 'run_date'], sort: 'id', ...page}));
    const demo = runs.filter((r) => /^Hepatocyte Stability, .+, \d{4}-\d{2}-\d{2}$/.test(r.name ?? ''));
    expect(demo.length > 0, true, 'no demo assay run on the stand');
    const off = demo.filter((r) => dayjs(r.run_date).toISOString() !== `${r.name!.slice(-10)}T00:00:00.000Z`)
      .map((r) => `${r.name}: ${dayjs(r.run_date).toISOString()}`);
    expectArray(off, []);
  });
});

category('MolTrack EMS: app', () => {
  test('the compounds view lists compounds beside the filter builder', async () => {
    const view = await MolTrackEms.view('compound');
    grok.shell.addView(view);
    try {
      await awaitCheck(() => view.root.querySelector('[data-u2="domain-list"] [data-u2-row]') !== null,
        'the list never showed a compound', 10000);
      expect(view.root.querySelector('.u2-domain-app-filters-pane [data-u2="filter-builder"]') !== null, true);
    } finally {
      view.close();
    }
  });

  test('every table of the schema opens as a view', async () => {
    for (const table of ['compound_batch', 'addition', 'assay', 'assay_run', 'assay_result']) {
      const view = grok.shell.addView(await MolTrackEms.view(table));
      try {
        const app = DomainApp.of(view)!;
        await awaitCheck(() => app.listSource.state.value === 'ready' &&
          view.root.querySelector('[data-u2="domain-list"] [data-u2-row], .u2-domain-list-empty') !== null,
        `${table} never rendered its list`, 15000);
        expect(app.listSource.error.value === undefined, true, String(app.listSource.error.value));
      } finally {
        view.close();
      }
    }
  });

  test('New on a compound\'s Batches tab opens the batch create page on that compound', async () => {
    const compounds = grok.shell.addView(await MolTrackEms.view('compound'));
    const app = DomainApp.of(compounds)!;
    let batches: DG.ViewBase | undefined;
    try {
      await awaitCheck(() => app.listSource.rows.items.value.length > 0, 'the list never loaded', 15000);
      const compound = app.listSource.rows.items.value[0];
      await app.goTo('entity', compound.id);
      const button = () => compounds.root.querySelector<HTMLButtonElement>(
        '[data-u2-part="panes"] .u2-domain-children-toolbar [data-u2="new-button"]');
      await awaitCheck(() => button() !== null && !button()!.disabled, 'the Batches tab has no New', 15000);
      button()!.click();
      const base = MolTrackEms.baseOf(VIEWS[1]);
      await awaitCheck(() => [...DomainApp.live].some((a) => a.base === base && a.entity.value === DomainApp.NEW &&
        a.entitySource.value?.currentRow.value?.compound_id === compound.id),
      'no batch create page on the compound', 15000);
      const batchesApp = [...DomainApp.live].find((a) => a.base === base)!;
      batches = DomainApp.openAt(base, batchesApp.path.peek());
      expect(app.page.value, 'entity', 'the compound stays where it was');
      batchesApp.session.discard();
    } finally {
      batches?.close();
      compounds.close();
    }
  });

  test('the tree reuses the open view: a record, then the list, then the create page', async () => {
    const base = MolTrackEms.baseOf(VIEWS[0]);
    const apps = () => [...DomainApp.live].filter((a) => a.base === base);
    const view = await MolTrackEms.show('Compounds');
    try {
      const app = apps()[0];
      await awaitCheck(() => app.listSource.rows.items.value.length > 0, 'the list never loaded', 15000);
      await app.goTo('entity', app.listSource.rows.items.value[0].id);
      expect((await MolTrackEms.show('Compounds')).dart === view.dart, true);
      await awaitCheck(() => app.page.value === 'list', 'Search › Compounds left the record open', 5000);
      await MolTrackEms.show('Compounds?entity=new');
      await awaitCheck(() => app.entity.value === DomainApp.NEW, 'Register › Compound never opened the create page',
        5000);
      expect(apps().length, 1, 'views of the compounds');
    } finally {
      view.close();
    }
  });

  test('compounds and batches list newest first, and a row moved to shows its compound in the context panel',
    async () => {
      for (const [table, key] of [['compound_batch', 'batch_regno'], ['compound', 'molregno']]) {
        const view = await MolTrackEms.view(table);
        grok.shell.addView(view);
        const app = DomainApp.of(view)!;
        try {
          await awaitCheck(() => app.listSource.rows.items.value.length > 1, `the ${table} list never loaded`, 15000);
          const [first, second] = app.listSource.rows.items.value;
          expect((first[key] as number) > (second[key] as number), true,
            `${table}: ${first[key]} before ${second[key]}`);
          app.list.root.querySelector<HTMLElement>('[tabindex]')!.focus();
          app.listSource.currentRow.value = second;
          const corporateId = table === 'compound' ?
            MolTrackTables.corporateId(await MolTrackTables.compounds(), 'molregno', second.molregno) :
            String(second[Rows.caption('compound_id')]);
          await awaitCheck(() => grok.shell.o instanceof HTMLElement && grok.shell.o.textContent!.includes(corporateId),
            `the context panel never showed the compound of the ${table}`, 15000);
        } finally {
          view.close();
        }
      }
    });

  test('a batch page opens its compound: from the compound card and from the reference', async () => {
    const batches = await MolTrackEms.show('Batches');
    const app = DomainApp.of(batches)!;
    const base = MolTrackEms.baseOf(VIEWS[0]);
    const compounds = () => [...DomainApp.live].find((a) => a.base === base);
    const close = () => [...grok.shell.views]
      .find((v) => compounds() !== undefined && v.root.contains(compounds()!.root))?.close();
    try {
      await awaitCheck(() => app.listSource.rows.items.value.length > 0, 'the batches never loaded', 15000);
      const batch = app.listSource.rows.items.value[0];
      await app.goTo('entity', batch.id);
      // the card with the compounds open behind the batch: the click must not land back on the batch view; the
      // reference with none open: MolTrack's opener docks one
      for (const [where, target, behind] of [['card', '.moltrack-ems-compound-section .moltrack-ems-compound a', true],
        ['reference', '[data-u2-name="compound_id"] [data-u2-part="open"], [data-u2-name="compound_id"] a', false]]) {
        close();
        if (behind)
          await MolTrackEms.show('Compounds');
        grok.shell.v = batches;
        await awaitCheck(() => batches.root.querySelector(target as string) !== null,
          `the ${where} has nothing to click`, 15000);
        batches.root.querySelector<HTMLElement>(target as string)!.click();
        await awaitCheck(() => compounds() !== undefined && grok.shell.v?.root.contains(compounds()!.root) === true &&
          compounds()!.entity.value === batch.compound_id,
        `the ${where} did not bring the compound forward`, 15000);
      }
      // the compound a batch was registered on is a link to it once saved, never a picker
      expect(batches.root.querySelector('[data-u2-name="compound_id"] input'), null);
    } finally {
      close();
      batches.close();
    }
  });

  test('the tree registers compounds and batches, and no assay', async () => {
    const root = ui.tree();
    MolTrackEms.treeBrowser(root);
    expectArray(root.getOrCreateGroup('Register').items.map((node) => node.text), ['Compound', 'Batch', 'Bulk']);
  });

  test('a compound page reads structure and names first, the registration internals folded under them', async () => {
    const view = await MolTrackEms.view('compound');
    grok.shell.addView(view);
    const app = DomainApp.of(view)!;
    try {
      await awaitCheck(() => app.listSource.rows.items.value.length > 0, 'the list never loaded', 15000);
      await app.goTo('entity', app.listSource.rows.items.value[0].id);
      await awaitCheck(() => app.form.value?.form != null, 'the record never showed', 10000);
      const main = app.form.value!.form!.properties.map((p) => p.name);
      expectArray(main, CompoundsApp.FIELDS);
      expect(app.form.value!.input('original_structure') === undefined, true, 'a saved structure is not an input');
      const details = view.root.querySelector<HTMLElement>(
        '[data-u2-part="entity-page"] .u2-section:not(.moltrack-ems-compound-section)')!;
      expect(details.querySelector('.u2-section-title')!.textContent, 'Registration details');
      await awaitCheck(() => details.querySelector('[data-u2-name="hash_mol"]') !== null,
        'the registration hash is not under Registration details', 10000);
      expect(details.querySelector<HTMLElement>('[data-u2-part="body"]')!.hidden ||
        getComputedStyle(details.querySelector<HTMLElement>('[data-u2-part="body"]')!).display === 'none', true);
      expect((app.ribbon().tools[0] as DomainSearch).tooltipText, CompoundsApp.SEARCH_HINT);
    } finally {
      view.close();
    }
  });

  test('list lines name what a row belongs to: an addition its role, a batch its compound, a result both', async () => {
    const tables = await Promise.all(['addition', 'compound_batch', 'assay_result']
      .map((t) => MolTrackTables.table(t)));
    const sources = [tables[0].source({query: 'role = "SALT"', pageSize: 1}), tables[1].source({pageSize: 1}),
      tables[2].source({pageSize: 1})];
    const line = (at: number) => tables[at].renderer.listItem!(sources[at].rows.items.value[0]);
    try {
      await awaitCheck(() => sources.every((s) => s.rows.items.value.length > 0), 'a list never loaded', 15000);
      expect(line(0).textContent!.includes('Salt'), true, line(0).textContent!);
      // a reference's caption is rewritten in its target's display format once that is known
      await awaitCheck(() => /^DGB-\d{6}DG-\d{6}/.test(line(1).textContent!),
        `the batch line names no compound: ${line(1).textContent}`, 10000);
      await awaitCheck(() => /^DGB-\d{6}$/.test(String(sources[2].rows.items.value[0][Rows.caption('batch_id')])),
        'the result\'s batch caption was never formatted', 10000);
      const result = line(2);
      await awaitCheck(() => /^DG-\d{6}DGB-\d{6} · /.test(result.textContent!),
        `the result line names no compound and batch: ${result.textContent}`, 10000);
    } finally {
      for (const source of sources)
        source.dispose();
    }
  });

  test('a column of Grok IDs is detected without the app having been opened', async () => {
    const detect = DG.Func.find({package: 'Moltrack', name: 'detectGrokId'})[0];
    expect(detect !== undefined, true, 'no Grok ID detector in detectors.js');
    expect(await detect.apply({col: DG.Column.fromStrings('id', ['DG-000001', 'DG-000012', ''])}), 'Grok ID');
    expect(await detect.apply({col: DG.Column.fromStrings('id', ['DGB-000001'])}), null);
  });

  test('a draft derives its identity from the structure as submitted, and loses it when cleared', async () => {
    const compounds = await MolTrackTables.compounds();
    expect(compounds.derived.targets.includes('hash_mol'), true);
    const draft = compounds.draft();
    try {
      await awaitCheck(() => draft.currentRow.value !== null, 'no draft row', 5000);
      const row = draft.currentRow.value!;
      row.original_structure = 'CC(=O)O';
      await awaitCheck(() => row.structure === 'CC(=O)O', 'the structure was never derived', 10000);
      expect(row.formula, 'C2H4O2');
      expect(row.match_rule, MolTrackTables.matchRule());
      row.original_structure = '';
      await awaitCheck(() => row.hash_mol === null || row.hash_mol === undefined,
        'the hash outlived its structure', 5000);
    } finally {
      draft.discard();
      draft.dispose();
    }
  });
});
