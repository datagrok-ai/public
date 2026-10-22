import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, awaitCheck, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {domains, DomainApp, DomainPick, IMPORT_COLUMNS} from '@datagrok-libraries/ems/src/dg';
import {Filters} from '@datagrok-libraries/u2';
import {moltrackEmsDb} from '../generated/db';
import {Chemistry, MatchRule} from '../ems/chemistry';
import {MolTrackEms, VIEWS} from '../ems/app';
import {MolTrackTables, SCHEMA} from '../ems/tables';
import {StructureSearch} from '../ems/structure-search';

// Expected values are RDKit 2026.03.6's (rdkit.Chem.RegistrationHash, MolToInchiKey, CalcMolFormula)
// for the same inputs — what Crux reproduces.
category('MolTrack EMS: chemistry', () => {
  test('aspirin and its sodium salt: one parent and hash; the rule decides on stereo and tautomers', async () => {
    const acid = await Chemistry.identify('CC(=O)Oc1ccccc1C(=O)O', 'ALL_LAYERS');
    const salt = await Chemistry.identify('CC(=O)Oc1ccccc1C(=O)[O-].[Na+]', 'ALL_LAYERS');
    expect(acid.structure, 'CC(=O)Oc1ccccc1C(=O)O');
    expect(acid.hash_mol, '60ec019cc6977893a0e2fe68e05fb6f163dc98a3');
    expect(salt.hash_mol, acid.hash_mol);
    expectArray([acid.inchikey, acid.formula, salt.inchikey, salt.formula],
      ['BSYNRYMUTXBXSQ-UHFFFAOYSA-N', 'C9H8O4', 'JZLOKWGVGHYBKD-UHFFFAOYSA-M', 'C9H7NaO4']);
    const same = async (a: string, b: string, rule: MatchRule) =>
      (await Chemistry.identify(a, rule)).hash_mol === (await Chemistry.identify(b, rule)).hash_mol;
    expectArray([await same('C[C@H](N)C(=O)O', 'C[C@@H](N)C(=O)O', 'ALL_LAYERS'),
      await same('C[C@H](N)C(=O)O', 'C[C@@H](N)C(=O)O', 'STEREO_INSENSITIVE_LAYERS'),
      await same('O=c1cccc[nH]1', 'Oc1ccccn1', 'ALL_LAYERS'),
      await same('O=c1cccc[nH]1', 'Oc1ccccn1', 'TAUTOMER_INSENSITIVE_LAYERS')], [false, true, false, true]);
  });
});

// Real rows in moltrack_ems, soft-deleted in after(); each test registers its own unusual structures.
category('MolTrack EMS: registration', () => {
  const tag = `T${Date.now()}`;
  const chain = (n: number, end = 'C(=O)O') => `${'C'.repeat(n)}${end}`;
  const created = {compounds: new Set<string>(), batches: new Set<string>()};
  const compoundOf = async (grokId: string) =>
    (await moltrackEmsDb.compounds.first({filter: `molregno = ${Number(DomainPick.digits(grokId))}`}))!;
  const results = (report: DG.DataFrame) => report.col(IMPORT_COLUMNS.result)!.toList();
  const messages = (report: DG.DataFrame) => report.col(IMPORT_COLUMNS.message)!.toList();
  // an empty cell is a hole in the list
  const grokIds = (report: DG.DataFrame) => Array.from(report.col('Grok ID')!.toList(), (id) => id ?? '');
  /** The compounds a run added — not the ones it found registered already, which may be the demo's — are deleted. */
  const register = async (rows: Record<string, string>[], rejectAll = false) => {
    const report: DG.DataFrame = await grok.functions.call('Moltrack:registerEmsBulk',
      {table: DG.DataFrame.fromObjects(rows)!, rejectAll});
    for (const [i, result] of results(report).entries()) {
      if (result === 'Added')
        created.compounds.add((await compoundOf(grokIds(report)[i])).id);
    }
    return report;
  };

  after(async () => {
    const failed: string[] = [];
    for (const id of created.batches)
      await moltrackEmsDb.compoundBatches.delete(id).catch((e) => failed.push(`batch ${id}: ${e}`));
    for (const id of created.compounds)
      await moltrackEmsDb.compounds.delete(id).catch((e) => failed.push(`compound ${id}: ${e}`));
    if (failed.length > 0)
      throw new Error(`Cleanup failed: ${failed.join('; ')}`);
  });

  test('bulk: a new structure is added, a registered one skipped by its Grok ID, a broken one and one the server ' +
    'refuses fail alone; with reject all nothing is written', async () => {
    const acidId = grokIds(await register([{'smiles': chain(25), 'Common Name': `${tag} acid`}]))[0];
    expect(/^DG-\d{6}$/.test(acidId), true, `corporate compound id ${acidId}`);
    expect((await MolTrackEms.compoundPanel(acidId)).root.textContent!.includes(`${tag} acid`), true);
    const rows = [{'smiles': chain(26), 'Common Name': `${tag} new`, 'MolLogP': ''},
      {'smiles': chain(25, 'C(=O)[O-].[Na+]'), 'Common Name': `${tag} salt`, 'MolLogP': ''},
      {'smiles': 'C1CC', 'Common Name': `${tag} broken`, 'MolLogP': ''},
      {'smiles': chain(27), 'Common Name': `${tag} refused`, 'MolLogP': 'high'}];
    expectArray(results(await register(rows, true)), ['Not imported', 'Skipped', 'Error', 'Not imported']);
    expect(await moltrackEmsDb.compounds.count({filter: `common_name = "${tag} new"`}), 0, 'reject all wrote nothing');
    const report = await register(rows);
    expectArray(results(report), ['Added', 'Skipped', 'Error', 'Error']);
    const [, registered, broken, refused] = messages(report);
    expect(registered.endsWith(MolTrackTables.alreadyRegistered(acidId)), true, registered);
    expect(broken.includes('Cannot read the structure'), true, broken);
    expect(/MolLogP/.test(refused), true, `the refusal names the column: ${refused}`);
    expect((await compoundOf(grokIds(report)[0])).common_name, `${tag} new`);
  });

  test('every rule\'s hash is stored, so a stricter rule today and a looser one later still match', async () => {
    const l = `C[C@H](N)${'C'.repeat(20)}C(=O)O`;
    const d = `C[C@@H](N)${'C'.repeat(20)}C(=O)O`;
    expectArray(results(await register([{smiles: l}])), ['Added']);
    const looser = await Chemistry.identify(d, 'STEREO_INSENSITIVE_LAYERS');
    const found = await MolTrackTables.registered('hash_stereo_insensitive', [looser.hash_mol]);
    expect(found.size, 1, 'the enantiomer is the same compound once stereo is ignored');
    const stricter = await Chemistry.identify(d, 'ALL_LAYERS');
    expect((await MolTrackTables.registered('hash_all_layers', [stricter.hash_mol])).size, 0);
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

  test('a new compound drawn as a registered structure is refused on its Structure, naming the compound, ' +
    'and cannot be saved', async () => {
    const corporateId = grokIds(await register([{smiles: chain(36)}]))[0];
    const view = grok.shell.addView(await MolTrackEms.views.view(`${SCHEMA}.compound`));
    const app = DomainApp.of(view)!;
    try {
      await app.create();
      await awaitCheck(() => app.entitySource.value?.currentRow.value != null, 'no new compound', 10000);
      app.entitySource.value!.currentRow.value!.original_structure = chain(36, 'C(=O)[O-].[Na+]');
      const message = MolTrackTables.alreadyRegistered(corporateId);
      await awaitCheck(() => app.form.value?.input('original_structure')?.validity.value === message,
        'the Structure input never said the structure is registered', 15000);
      await awaitCheck(() => (app.session.primaryButton as HTMLButtonElement | undefined)?.disabled === true,
        'Save stayed enabled under the verdict', 5000);
      expect(await app.session.save(), false);
      const problem = app.form.value!.problem.value ?? '';
      expect(problem.includes(message) && !/hash/i.test(problem), true, `the refusal reads "${problem}"`);
    } finally {
      app.session.discard();
      view.close();
    }
  });

  test('a batch import names each batch\'s compound by its Grok ID, one no compound has fails alone; ' +
    'the server keeps a batch on its compound', async () => {
    const grokId = grokIds(await register([{'smiles': chain(36, 'C(=O)N'), 'Common Name': `${tag} imported`}]))[0];
    const df = DG.DataFrame.fromObjects([{'Grok ID': grokId, 'EPA Batch ID': `${tag}-imp-1`, 'Purity': '99'},
      {'Grok ID': 'DG-999999', 'EPA Batch ID': `${tag}-imp-2`, 'Purity': '98'}])!;
    const report = await domains.importRows(await MolTrackTables.batches(), df,
      {allOrNothing: false, mapping: {compound_id: 'Grok ID'}});
    const batches = await moltrackEmsDb.compoundBatches.query({filter: `epa_batch_id starts "${tag}-imp"`});
    for (const batch of batches)
      created.batches.add(batch.id);
    expectArray(results(report), ['Added', 'Error']);
    expect(messages(report)[1].includes('No match for "DG-999999"'), true, messages(report)[1]);
    expectArray(batches.map((b) => b.epa_batch_id), [`${tag}-imp-1`]);
    expect(batches[0].compound_id, (await compoundOf(grokId)).id);
    let refusal = '';
    await moltrackEmsDb.compoundBatches.update(batches[0].id,
      {compound_id: (await compoundOf(grokIds(await register([{smiles: chain(37)}]))[0])).id})
      .catch((e) => refusal = String((e as Error).message ?? e));
    expect(/write-once/.test(refusal), true, `a batch was moved to another compound: ${refusal}`);
  }, {timeout: 120000});
});

category('MolTrack EMS: app', () => {
  test('every table of the schema opens as a view, its list beside the filter builder', async () => {
    for (const entry of VIEWS) {
      const view = grok.shell.addView(await MolTrackEms.views.view(`${SCHEMA}.${entry.table}`));
      try {
        const app = DomainApp.of(view)!;
        await awaitCheck(() => app.listSource.state.value === 'ready' &&
          view.root.querySelector('[data-u2="domain-list"] [data-u2-row], .u2-domain-list-empty') !== null,
        `${entry.table} never rendered its list`, 15000);
        expect(app.listSource.error.value === undefined, true, String(app.listSource.error.value));
        expect(view.root.querySelector('.u2-domain-app-filters-pane [data-u2="filter-builder"]') !== null, true,
          `${entry.table} has no filter builder`);
      } finally {
        view.close();
      }
    }
  });

  test('New on a compound\'s Batches tab opens the batch create page on that compound', async () => {
    const compounds = grok.shell.addView(await MolTrackEms.views.view(`${SCHEMA}.compound`));
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
      const base = MolTrackEms.views.baseOf(`${SCHEMA}.compound_batch`);
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
