import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {DomainApp, DomainPick, DomainSearch} from '@datagrok-libraries/ems/src/dg';
import {Filters, Rows} from '@datagrok-libraries/u2';
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

  const remember = async (report: DG.DataFrame) => {
    const track = async (column: string, client: DG.DomainTableClient, field: string, ids: Set<string>) => {
      for (const id of report.col(column)?.toList() ?? []) {
        const row = id ? await client.first({filter: `${field} = ${Number(DomainPick.digits(id))}`}) : null;
        if (row !== null)
          ids.add(row.id);
      }
    };
    await track(Registration.COMPOUND_ID, moltrackEmsDb.compounds, 'molregno', created.compounds);
    await track(Registration.BATCH_ID, moltrackEmsDb.compoundBatches, 'batch_regno', created.batches);
  };
  const register = async (rows: Record<string, string>[], batches: boolean, rejectAll = false) => {
    const report = await MolTrackEms.register(DG.DataFrame.fromObjects(rows)!, batches, rejectAll);
    await remember(report);
    return report;
  };
  const statuses = (report: DG.DataFrame) => report.col(Registration.STATUS)!.toList();
  const messages = (report: DG.DataFrame) => report.col(Registration.MESSAGE)!.toList();

  after(async () => {
    const failed: string[] = [];
    for (const id of created.batches)
      await moltrackEmsDb.compoundBatches.delete(id).catch((e) => failed.push(`batch ${id}: ${e}`));
    for (const id of created.compounds)
      await moltrackEmsDb.compounds.delete(id).catch((e) => failed.push(`compound ${id}: ${e}`));
    if (failed.length > 0)
      throw new Error(`Cleanup failed: ${failed.join('; ')}`);
  });

  test('a parent and its salt register as one compound with two batches', async () => {
    const report = await register([
      {'smiles': chain(25), 'EPA Batch ID': `${tag}-1`, 'Purity': '99.5', 'Common Name': 'Pentacosanoic acid'},
      {'smiles': chain(25, 'C(=O)[O-].[Na+]'), 'EPA Batch ID': `${tag}-2`, 'Purity': '97.25',
        'Common Name': 'Sodium pentacosanoate'},
    ], true);
    expectArray(statuses(report), ['success', 'success']);
    const ids = report.col(Registration.COMPOUND_ID)!.toList();
    expect(ids[0], ids[1]);
    expect(/^DG-\d{6}$/.test(ids[0]), true, `corporate compound id ${ids[0]}`);
    expect(/^DGB-\d{6}$/.test(report.col(Registration.BATCH_ID)!.get(1)), true);
    const batch = await moltrackEmsDb.compoundBatches.first({filter: `epa_batch_id = "${tag}-2"`});
    expect(batch!.purity, 97.25);
    const compound = await moltrackEmsDb.compounds.get(batch!.compound_id);
    expect(compound.common_name, 'Pentacosanoic acid', 'the first row of a structure carries the compound');
    expect(compound.structure, chain(25));
  });

  test('a registered structure registered again as a compound is refused, naming the compound', async () => {
    const first = await register([{smiles: chain(26)}], false);
    expectArray(statuses(first), ['success']);
    const again = await register([{smiles: chain(26)}], false);
    expectArray(statuses(again), ['failed']);
    const corporateId = first.col(Registration.COMPOUND_ID)!.get(0);
    expect(messages(again)[0], `This structure is already registered as ${corporateId}`);
  });

  test('a bad row fails alone; with reject all nothing is written', async () => {
    const strict = await register([{'smiles': 'C1CC', 'EPA Batch ID': `${tag}-3`},
      {'smiles': chain(27, 'N'), 'EPA Batch ID': `${tag}-4`}], true, true);
    expectArray(statuses(strict), ['failed', 'not_processed']);
    expect(await moltrackEmsDb.compoundBatches.count({filter: `epa_batch_id = "${tag}-4"`}), 0, 'nothing was written');
    const lenient = await register([{'smiles': 'C1CC', 'EPA Batch ID': `${tag}-3`},
      {'smiles': chain(27, 'N'), 'EPA Batch ID': `${tag}-4`}], true);
    expectArray(statuses(lenient), ['failed', 'success']);
    expect(messages(lenient)[0], 'Cannot read the structure');
  });

  test('reject all is one transaction: a batch the server refuses leaves no new compound behind', async () => {
    const name = `${tag} atomic`;
    const report = await register([
      {'smiles': chain(28), 'Common Name': name, 'EPA Batch ID': `${tag}-5`, 'Purity': '50'},
      {'smiles': chain(29), 'EPA Batch ID': `${tag}-6`, 'Purity': '150'},
    ], true, true);
    expectArray(statuses(report), ['not_processed', 'failed']);
    expect(/Purity/.test(messages(report)[1]), true, `the refusal names the column: ${messages(report)[1]}`);
    expect(await moltrackEmsDb.compounds.count({filter: `common_name = "${name}"`}), 0);
    expect(await moltrackEmsDb.compoundBatches.count({filter: `epa_batch_id = "${tag}-5"`}), 0);
  });

  test('a row with no structure fails alone', async () => {
    const report = await register([{'smiles': '', 'EPA Batch ID': `${tag}-7`},
      {'smiles': chain(30), 'EPA Batch ID': `${tag}-8`}], true);
    expectArray(statuses(report), ['failed', 'success']);
    expect(messages(report)[0], 'Structure is required');
  });

  test('every rule\'s hash is stored, so a stricter rule today and a looser one later still match', async () => {
    const l = `C[C@H](N)${'C'.repeat(20)}C(=O)O`;
    const d = `C[C@@H](N)${'C'.repeat(20)}C(=O)O`;
    const report = await register([{smiles: l}], false);
    expectArray(statuses(report), ['success']);
    const looser = await Chemistry.identify(d, 'STEREO_INSENSITIVE_LAYERS');
    const found = await MolTrackTables.registered('hash_stereo_insensitive', [looser.hash_mol]);
    expect(found.size, 1, 'the enantiomer is the same compound once stereo is ignored');
    const stricter = await Chemistry.identify(d, 'ALL_LAYERS');
    expect((await MolTrackTables.registered('hash_all_layers', [stricter.hash_mol])).size, 0);
  });

  test('a corporate compound ID opens its context panel', async () => {
    const report = await register([{'smiles': chain(31), 'Common Name': `${tag} panel`}], false);
    const panel = await MolTrackEms.compoundPanel(report.col(Registration.COMPOUND_ID)!.get(0));
    expect(panel.root.textContent!.includes(`${tag} panel`), true);
    expect((await MolTrackEms.compoundPanel('DG-999999')).root.textContent,
      'DG-999999 is not registered in MolTrack EMS');
  });

  test('the structure conditions find a registered compound, and a list filtered by one holds it', async () => {
    await register([{smiles: chain(32)}], false);
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
    await register([{smiles: chain(33)}], false);
    const df = DG.DataFrame.fromObjects([{smiles: chain(34)}, {smiles: chain(33, 'C(=O)[O-].[Na+]')},
      {smiles: 'C1CC('}, {smiles: chain(34)}])!;
    const mapping = Registration.autoMap(df, await MolTrackTables.compounds(), await MolTrackTables.batches());
    const preview = await Registration.preview(df, mapping);
    expectArray(preview.map((p) => p.verdict), ['new', 'existing', 'error', 'new']);
    expect(preview[1].message.includes('DG-'), true);
    expect((await MolTrackTables.registered('hash_mol',
      [(await Chemistry.identify(chain(34), MolTrackTables.matchRule())).hash_mol])).size, 0);
  });

  test('a new compound drawn as a registered structure is refused on its Structure, naming the compound', async () => {
    const report = await register([{smiles: chain(36)}], false);
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

  test('a new batch registers its new compound in the same save, and a second batch joins it', async () => {
    const view = await MolTrackEms.view('compound_batch');
    grok.shell.addView(view);
    const app = DomainApp.of(view)!;
    const draw = async () => {
      await app.create();
      await awaitCheck(() => app.inline.value.get('compound_id')?.currentRow.value != null, 'no new compound', 10000);
      const compound = app.inline.value.get('compound_id')!.currentRow.value!;
      compound.original_structure = chain(35, 'O');
      return compound;
    };
    try {
      const compound = await draw();
      await awaitCheck(() => app.entitySource.value?.currentRow.value?.compound_id === compound.id,
        'the batch never pointed at the new compound', 10000);
      app.entitySource.value!.currentRow.value!.epa_batch_id = `${tag}-inline-1`;
      expect(await app.session.save(), true);
      created.batches.add(app.entitySource.value!.currentRow.value!.id);
      expect([...document.querySelectorAll('.d4-balloon')].some((b) => b.textContent!.includes('not available while')),
        false, 'the save discards nothing while it runs');
      const hash = (await Chemistry.identify(chain(35, 'O'), MolTrackTables.matchRule())).hash_mol;
      const registered = (await MolTrackTables.registered('hash_mol', [hash])).get(hash)!;
      created.compounds.add(registered.id);
      await draw();
      await awaitCheck(() => app.entitySource.value?.currentRow.value?.compound_id === registered.id,
        'the second batch was not added to the registered compound', 20000);
      await awaitCheck(() => app.inline.value.get('compound_id')?.currentRow.value != null,
        'no new compound after', 5000);
      expect(app.inline.value.get('compound_id')!.isDirty.value, false, 'the dropped compound is not a change');
      const corporateId = MolTrackTables.corporateId(await MolTrackTables.compounds(), 'molregno', registered.molregno);
      const section = () => view.root.querySelector<HTMLElement>('.moltrack-ems-compound-section')!;
      const matched = () => section().querySelector<HTMLElement>('.moltrack-ems-matched')!;
      await awaitCheck(() => !section().hidden && section().textContent!.includes(corporateId) && !matched().hidden,
        'the page never showed the compound the new batch joined', 10000);
      expect(matched().textContent!.includes(`This structure is registered as ${corporateId}`), true);
      expect(app.form.value!.root.closest('.u2-section')?.querySelector('.u2-section-title')?.textContent, 'Batch');
      // "Use a different structure" takes the batch off the compound, and drawing it again attaches it again
      matched().querySelector<HTMLElement>('a')!.click();
      await awaitCheck(() => app.entitySource.value!.currentRow.value!.compound_id == null && section().hidden === true,
        'the batch stayed on the registered compound', 5000);
      await awaitCheck(() => app.inline.value.get('compound_id')?.currentRow.value != null, 'no new compound', 5000);
      app.inline.value.get('compound_id')!.currentRow.value!.original_structure = chain(35, 'O');
      await awaitCheck(() => app.entitySource.value?.currentRow.value?.compound_id === registered.id,
        'drawn again, the structure was not attached again', 20000);
      app.entitySource.value!.currentRow.value!.epa_batch_id = `${tag}-inline-2`;
      expect(await app.session.save(), true);
      const saved = app.entitySource.value!.currentRow.value!.id;
      created.batches.add(saved);
      const batches = await moltrackEmsDb.compoundBatches.query(
        {filter: MolTrackTables.within('compound_id', [registered.id])});
      expectArray(batches.map((b) => b.epa_batch_id).sort(), [`${tag}-inline-1`, `${tag}-inline-2`]);
      // the saved batch keeps showing its compound, and nothing moves it to another one
      await awaitCheck(() => matched().hidden === true, 'the registration note outlived the save', 5000);
      expect(!section().hidden && section().textContent!.includes(corporateId), true);
      await register([{smiles: chain(37)}], false);
      const otherHash = (await Chemistry.identify(chain(37), MolTrackTables.matchRule())).hash_mol;
      const other = (await MolTrackTables.registered('hash_mol', [otherHash])).get(otherHash)!;
      let refusal = '';
      await moltrackEmsDb.compoundBatches.update(saved, {compound_id: other.id})
        .catch((e) => refusal = String((e as Error).message ?? e));
      expect(/write-once/.test(refusal), true, `a batch was moved to another compound: ${refusal}`);
    } finally {
      view.close();
    }
  });

  const amide = (n: number) => chain(n, 'C(=O)N');

  test('with reject all, a row left unregistered names the row that stopped the run', async () => {
    const report = await register([{'smiles': amide(30), 'EPA Batch ID': `${tag}-9`},
      {'smiles': 'C1CC', 'EPA Batch ID': `${tag}-10`}], true, true);
    expectArray(statuses(report), ['not_processed', 'failed']);
    expectArray(messages(report), ['Not registered: row 2 has a problem', 'Cannot read the structure']);
  });

  test('the preview tells a structure the file repeats from a new compound', async () => {
    const df = DG.DataFrame.fromObjects([{smiles: amide(31)}, {smiles: amide(32)}, {smiles: amide(31)}])!;
    const preview = await Registration.preview(df,
      Registration.autoMap(df, await MolTrackTables.compounds(), await MolTrackTables.batches()));
    expectArray(preview.map((p) => p.sameAs), [undefined, undefined, 0]);
  });

  /** The bulk wizard over `rows`, open in the shell as the table a registrar would pick. */
  const bulkWizard = async (rows: Record<string, unknown>[]) => {
    const df = DG.DataFrame.fromObjects(rows)!;
    df.name = `${tag}-bulk`;
    const table = grok.shell.addTableView(df);
    const running = RegisterDialog.open(df, 'Batches');
    await awaitCheck(() => document.querySelector('.moltrack-ems-register') !== null, 'the bulk wizard never opened',
      30000);
    return {running, table, wizard: document.querySelector<HTMLElement>('.moltrack-ems-register')!};
  };
  const press = (wizard: HTMLElement, label: string) => [...wizard.querySelectorAll('button')]
    .find((b) => b.textContent!.trim().toUpperCase() === label)!.click();

  test('bulk registration rejects all by default, and refuses two file columns mapped onto one field', async () => {
    // a frame's column names ignore case: "Project" and "project_" both match the batch's Project
    const {running, table, wizard} = await bulkWizard([{'smiles': amide(33), 'Project': 'P1', 'project_': 'P2'}]);
    try {
      const onError = [...wizard.querySelectorAll('select')].find((s) => s.textContent!.includes('Reject all'))!;
      expect(onError.selectedOptions[0].textContent, 'Reject all');
      press(wizard, 'NEXT');
      const reason = () => wizard.querySelector('.u2-wizard-reason')?.textContent ?? '';
      await awaitCheck(() => reason().includes('Two columns of the file are mapped to Batch · Project'),
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

  test('bulk preview and report name each row by its own columns, and the report opens the batches registered',
    async () => {
      const {running, table, wizard} = await bulkWizard([
        {'smiles': amide(34), 'Common Name': `${tag} bulk one`, 'EPA Batch ID': `${tag}-b1`},
        {'smiles': amide(35), 'Common Name': `${tag} bulk two`, 'EPA Batch ID': `${tag}-b2`}]);
      const heads = (host: Element) => [...host.querySelectorAll('.u2-data-table-head')].map((h) => h.textContent);
      let view: DG.ViewBase | undefined;
      try {
        press(wizard, 'NEXT');
        press(wizard, 'NEXT');
        const preview = wizard.querySelector('.moltrack-ems-register-preview')!;
        await awaitCheck(() => preview.querySelectorAll('.u2-data-table-row').length === 2,
          'the preview never listed the rows', 30000);
        expectArray(heads(preview), ['#', 'Structure', 'Common Name', 'EPA Batch ID', 'Becomes', '']);
        expect(preview.textContent!.includes(`${tag} bulk one`) && preview.textContent!.includes(`${tag}-b2`), true);
        expect(preview.textContent!.includes('2 rows: 2 new compounds, 2 batches'), true, preview.textContent!);
        press(wizard, 'REGISTER');
        const report = wizard.querySelector('.moltrack-ems-register-report')!;
        await awaitCheck(() => report.textContent!.includes('Registered 2 of 2 rows'), 'the rows were not registered',
          60000);
        expectArray(heads(report), ['#', 'Common Name', 'EPA Batch ID', 'Result', 'Grok ID', 'Grok Batch ID', '']);
        const lines = [...report.querySelectorAll('.u2-data-table-row')].map((r) => r.textContent!);
        expect(lines.length === 2 && lines.every((l) => /Registered/.test(l) && /DG-\d{6}/.test(l) &&
          /DGB-\d{6}/.test(l)), true, lines.join(' / '));
        [...report.querySelectorAll<HTMLElement>('a, .d4-link-action')]
          .find((a) => a.textContent === 'Show the 2 batches registered')!.click();
        const registered = (await running)!;
        await remember(registered);
        const regnos = registered.col(Registration.BATCH_ID)!.toList().map((id) => Number(DomainPick.digits(id)));
        const shown = () => [...DomainApp.live].find((a) => a.base === MolTrackEms.baseOf(VIEWS[1]) &&
          String(a.listSource.query.value).startsWith('batch_regno in'));
        await awaitCheck(() => shown() !== undefined, 'the batches never opened', 30000);
        view = grok.shell.v;
        const app = shown()!;
        await awaitCheck(() => app.listSource.state.value === 'ready' &&
          app.listSource.rows.items.value.map((r) => r.batch_regno).sort().join() === regnos.sort().join(),
        'the batches list does not hold exactly the registered batches', 30000);
      } finally {
        if (document.contains(wizard))
          press(wizard, 'CANCEL');
        view?.close();
        table.close();
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
        expect((first[key] as number) > (second[key] as number), true, `${table}: ${first[key]} before ${second[key]}`);
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
