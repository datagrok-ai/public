import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
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
  test('aspirin and its sodium salt: one parent and hash; InChIKey and formula of the structure as given', async () => {
    const acid = await Chemistry.identify('CC(=O)Oc1ccccc1C(=O)O', 'ALL_LAYERS');
    const salt = await Chemistry.identify('CC(=O)Oc1ccccc1C(=O)[O-].[Na+]', 'ALL_LAYERS');
    expect(acid.structure, 'CC(=O)Oc1ccccc1C(=O)O');
    expect(acid.hash_mol, '60ec019cc6977893a0e2fe68e05fb6f163dc98a3');
    expect(salt.hash_mol, acid.hash_mol);
    expectArray([acid.inchikey, acid.formula, salt.inchikey, salt.formula],
      ['BSYNRYMUTXBXSQ-UHFFFAOYSA-N', 'C9H8O4', 'JZLOKWGVGHYBKD-UHFFFAOYSA-M', 'C9H7NaO4']);
  });

  test('the matching rule decides what is the same: stereo and tautomers', async () => {
    const same = async (a: string, b: string, rule: MatchRule) =>
      (await Chemistry.identify(a, rule)).hash_mol === (await Chemistry.identify(b, rule)).hash_mol;
    expectArray([await same('C[C@H](N)C(=O)O', 'C[C@@H](N)C(=O)O', 'ALL_LAYERS'),
      await same('C[C@H](N)C(=O)O', 'C[C@@H](N)C(=O)O', 'STEREO_INSENSITIVE_LAYERS'),
      await same('O=c1cccc[nH]1', 'Oc1ccccn1', 'ALL_LAYERS'),
      await same('O=c1cccc[nH]1', 'Oc1ccccn1', 'TAUTOMER_INSENSITIVE_LAYERS')], [false, true, false, true]);
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
});

// Real rows in moltrack_ems, soft-deleted in after(); each test registers its own unusual structures.
category('MolTrack EMS: registration', () => {
  const tag = `T${Date.now()}`;
  const chain = (n: number, end = 'C(=O)O') => `${'C'.repeat(n)}${end}`;
  const amide = (n: number) => chain(n, 'C(=O)N');
  const created = {compounds: new Set<string>(), batches: new Set<string>()};
  const compoundOf = async (grokId: string) =>
    (await moltrackEmsDb.compounds.first({filter: `molregno = ${Number(DomainPick.digits(grokId))}`}))!;
  /** The compounds a run registered — not the ones it found registered already, which may be the demo's. */
  const remember = async (report: DG.DataFrame) => {
    for (let i = 0; i < report.rowCount; i++) {
      if (report.get(Registration.STATUS, i) === 'success')
        created.compounds.add((await compoundOf(report.get(Registration.COMPOUND_ID, i))).id);
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

  after(async () => {
    const failed: string[] = [];
    for (const id of created.batches)
      await moltrackEmsDb.compoundBatches.delete(id).catch((e) => failed.push(`batch ${id}: ${e}`));
    for (const id of created.compounds)
      await moltrackEmsDb.compounds.delete(id).catch((e) => failed.push(`compound ${id}: ${e}`));
    if (failed.length > 0)
      throw new Error(`Cleanup failed: ${failed.join('; ')}`);
  });

  test('bulk: a structure registers once, its salt and a registered one are reported by Grok IDs that open a panel',
    async () => {
      const acid = await register([{'smiles': chain(25), 'Common Name': `${tag} acid`}]);
      expectArray(statuses(acid), ['success']);
      const acidId = grokIds(acid)[0];
      expect(/^DG-\d{6}$/.test(acidId), true, `corporate compound id ${acidId}`);
      expect((await MolTrackEms.compoundPanel(acidId)).root.textContent!.includes(`${tag} acid`), true);
      expect((await MolTrackEms.compoundPanel('DG-999999')).root.textContent,
        'DG-999999 is not registered in MolTrack EMS');
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
      expect((await moltrackEmsDb.compounds.query({filter: DG.cond('hash_mol', '=', hash)})).length, 1);
      const compound = await compoundOf(ids[0]);
      expect(compound.common_name, `${tag} new`, 'the first row of a structure carries the compound');
      expect(compound.structure, chain(26));
      expect((await compoundOf(acidId)).common_name, `${tag} acid`, 'a registered compound is not written');
    });

  test('a bad row fails alone; with reject all nothing is written and the others name the row that stopped it',
    async () => {
      const rows = [{'smiles': 'C1CC', 'Common Name': `${tag}-3`},
        {'smiles': chain(27, 'N'), 'Common Name': `${tag}-4`}, {'smiles': '', 'Common Name': `${tag}-5`}];
      const strict = await register(rows, true);
      expectArray(statuses(strict), ['failed', 'not_processed', 'failed']);
      expectArray(messages(strict),
        ['Cannot read the structure', 'Not registered: row 1 has a problem', 'Structure is required']);
      expect(await moltrackEmsDb.compounds.count({filter: `common_name = "${tag}-4"`}), 0, 'nothing was written');
      expectArray(statuses(await register(rows)), ['failed', 'success', 'failed']);
    });

  test('a compound the server refuses: with reject all nothing is written, with Skip the row it fails alone',
    async () => {
      const kept = `${tag} kept`;
      const rows = [{'smiles': chain(28), 'Common Name': kept, 'MolLogP': ''},
        {'smiles': chain(29), 'Common Name': `${tag} refused`, 'MolLogP': 'high'}];
      const strict = await register(rows, true);
      expectArray(statuses(strict), ['not_processed', 'failed']);
      expect(/MolLogP/.test(messages(strict)[1]), true, `the refusal names the column: ${messages(strict)[1]}`);
      expect(await moltrackEmsDb.compounds.count({filter: `common_name = "${kept}"`}), 0);
      const lenient = await register(rows);
      expectArray(statuses(lenient), ['success', 'failed']);
      expect(/MolLogP/.test(messages(lenient)[1]), true, `the refusal names the column: ${messages(lenient)[1]}`);
      expect((await compoundOf(grokIds(lenient)[0])).common_name, kept);
    });

  test('every rule\'s hash is stored, so a stricter rule today and a looser one later still match', async () => {
    const l = `C[C@H](N)${'C'.repeat(20)}C(=O)O`;
    const d = `C[C@@H](N)${'C'.repeat(20)}C(=O)O`;
    expectArray(statuses(await register([{smiles: l}])), ['success']);
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

  test('the preview says what each row becomes, tells a repeated structure, and writes nothing', async () => {
    await register([{smiles: chain(33)}]);
    const df = DG.DataFrame.fromObjects([{smiles: chain(34)}, {smiles: chain(33, 'C(=O)[O-].[Na+]')},
      {smiles: 'C1CC('}, {smiles: chain(34)}])!;
    const preview = await Registration.preview(df, Registration.autoMap(df, await MolTrackTables.compounds()));
    expectArray(preview.map((p) => p.verdict), ['new', 'existing', 'error', 'new']);
    expectArray(preview.map((p) => p.sameAs), [undefined, undefined, undefined, 0]);
    expect(preview[3].message, 'Same compound as row 1');
    expect(preview[1].message.includes('DG-'), true);
    expect((await MolTrackTables.registered('hash_mol',
      [(await Chemistry.identify(chain(34), MolTrackTables.matchRule())).hash_mol])).size, 0);
  });

  test('a new compound drawn as a registered structure is refused on its Structure, naming the compound, ' +
    'and cannot be saved', async () => {
    const corporateId = grokIds(await register([{smiles: chain(36)}]))[0];
    const view = grok.shell.addView(await MolTrackEms.view('compound'));
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

  test('a new batch on a picked compound shows it in the context panel, saves on it, stays on it and opens it',
    async () => {
      const corporateId = grokIds(await register([{'smiles': chain(35, 'O'), 'Common Name': `${tag} picked`}]))[0];
      const compound = await compoundOf(corporateId);
      const view = grok.shell.addView(await MolTrackEms.view('compound_batch'));
      const app = DomainApp.of(view)!;
      const compounds = () => [...DomainApp.live].find((a) => a.base === MolTrackEms.baseOf(VIEWS[0]));
      const link = () => view.root.querySelector<HTMLElement>(
        '[data-u2-name="compound_id"] [data-u2-part="open"], [data-u2-name="compound_id"] a');
      try {
        await app.create();
        await awaitCheck(() => app.entitySource.value?.currentRow.value != null, 'no new batch', 10000);
        expect(app.form.value!.input('compound_id') !== undefined, true, 'a new batch picks its compound');
        const batch = app.entitySource.value!.currentRow.value!;
        batch.compound_id = compound.id;
        batch.epa_batch_id = `${tag}-picked`;
        await awaitCheck(() => grok.shell.o instanceof HTMLElement && grok.shell.o.textContent!.includes(corporateId),
          'the context panel never showed the picked compound', 10000);
        expect(await app.session.save(), true);
        const saved = app.entitySource.value!.currentRow.value!.id;
        created.batches.add(saved);
        expect((await moltrackEmsDb.compoundBatches.get(saved)).compound_id, compound.id);
        await awaitCheck(() => app.form.value?.input('compound_id') === undefined && link() !== null,
          'the saved batch\'s compound is not a link to it', 10000);
        link()!.click();
        await awaitCheck(() => compounds()?.entity.value === compound.id &&
          grok.shell.v?.root.contains(compounds()!.root) === true, 'the link never brought the compound forward',
        15000);
        let refusal = '';
        await moltrackEmsDb.compoundBatches.update(saved,
          {compound_id: (await compoundOf(grokIds(await register([{smiles: chain(37)}]))[0])).id})
          .catch((e) => refusal = String((e as Error).message ?? e));
        expect(/write-once/.test(refusal), true, `a batch was moved to another compound: ${refusal}`);
      } finally {
        app.session.discard();
        [...grok.shell.views].find((v) => compounds() !== undefined && v.root.contains(compounds()!.root))?.close();
        view.close();
      }
    });

  test('bulk preview and report name each row by its own columns, and the report opens the compounds registered',
    async () => {
      const df = DG.DataFrame.fromObjects([{'smiles': amide(34), 'Common Name': `${tag} bulk one`},
        {'smiles': amide(35), 'Common Name': `${tag} bulk two`}])!;
      const table = grok.shell.addTableView(df);
      const running = RegisterDialog.open(df);
      await awaitCheck(() => document.querySelector('.moltrack-ems-register') !== null, 'no bulk wizard', 30000);
      const wizard = document.querySelector<HTMLElement>('.moltrack-ems-register')!;
      const press = (label: string) => [...wizard.querySelectorAll('button')]
        .find((b) => b.textContent!.trim().toUpperCase() === label)!.click();
      const heads = (host: Element) => [...host.querySelectorAll('.u2-data-table-head')].map((h) => h.textContent);
      let view: DG.ViewBase | undefined;
      try {
        press('NEXT');
        press('NEXT');
        const preview = wizard.querySelector('.moltrack-ems-register-preview')!;
        await awaitCheck(() => preview.querySelectorAll('.u2-data-table-row').length === 2,
          'the preview never listed the rows', 30000);
        expectArray(heads(preview), ['#', 'Structure', 'Common Name', 'Becomes', '']);
        expect(preview.textContent!.includes('2 rows: 2 new compounds'), true, preview.textContent!);
        press('REGISTER');
        const report = wizard.querySelector('.moltrack-ems-register-report')!;
        await awaitCheck(() => report.textContent!.includes('Registered 2 of 2 rows'), 'the rows were not registered',
          60000);
        expectArray(heads(report), ['#', 'Common Name', 'Result', 'Grok ID', '']);
        [...report.querySelectorAll<HTMLElement>('a, .d4-link-action')]
          .find((a) => a.textContent === 'Show the 2 compounds registered')!.click();
        const registered = (await running)!;
        await remember(registered);
        const molregnos = grokIds(registered).map((id) => Number(DomainPick.digits(id))).sort().join();
        const shown = () => [...DomainApp.live].find((a) => a.base === MolTrackEms.baseOf(VIEWS[0]) &&
          String(a.listSource.query.value).startsWith('molregno in'));
        await awaitCheck(() => shown() !== undefined, 'the compounds never opened', 30000);
        view = grok.shell.v;
        await awaitCheck(() => shown()!.listSource.state.value === 'ready' &&
          shown()!.listSource.rows.items.value.map((r) => r.molregno).sort().join() === molregnos,
        'the compounds list does not hold exactly the registered compounds', 30000);
      } finally {
        if (document.contains(wizard))
          press('CANCEL');
        view?.close();
        table.close();
      }
    }, {timeout: 120000});

  test('Batches import names each batch\'s compound by its Grok ID; a Grok ID no compound has fails alone',
    async () => {
      const grokId = grokIds(await register([{'smiles': amide(36), 'Common Name': `${tag} imported`}]))[0];
      const df = DG.DataFrame.fromObjects([{'Grok ID': grokId, 'EPA Batch ID': `${tag}-imp-1`, 'Purity': '99'},
        {'Grok ID': 'DG-999999', 'EPA Batch ID': `${tag}-imp-2`, 'Purity': '98'}])!;
      const named = (name: string) => Control.forElement(document.querySelector(`[data-u2-name="${name}"]`)) as
        unknown as {value: {value: unknown}} | undefined;
      const button = (label: string) => [...document.querySelectorAll<HTMLButtonElement>('.u2-dialog button')]
        .find((b) => b.textContent === label)!;
      const text = () => document.querySelector('.u2-domain-import-report')?.textContent ?? '';
      try {
        const running = domains.import(await MolTrackTables.batches(), {source: df});
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
});

category('MolTrack EMS: app', () => {
  test('every table of the schema opens as a view, its list beside the filter builder', async () => {
    for (const entry of VIEWS) {
      const view = grok.shell.addView(await MolTrackEms.view(entry.table));
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
        const view = grok.shell.addView(await MolTrackEms.view(table));
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

  test('a compound page reads structure and names first, the registration identity folded in its form', async () => {
    const view = grok.shell.addView(await MolTrackEms.view('compound'));
    const app = DomainApp.of(view)!;
    try {
      await awaitCheck(() => app.listSource.rows.items.value.length > 0, 'the list never loaded', 15000);
      await app.goTo('entity', app.listSource.rows.items.value[0].id);
      await awaitCheck(() => app.form.value?.form != null, 'the record never showed', 10000);
      expect(app.form.value!.input('original_structure') === undefined, true, 'a saved structure is not an input');
      const details = [...app.form.value!.root.querySelectorAll<HTMLElement>('.u2-section')]
        .find((s) => s.querySelector('.u2-section-title')?.textContent === 'Registration details');
      expect(details !== undefined, true, 'the form has no Registration details');
      const inDetails = (columns: string[]) => columns.filter((c) => details!.querySelector(`[data-u2-name="${c}"]`));
      expectArray(inDetails(CompoundsApp.DETAILS), CompoundsApp.DETAILS);
      expectArray(inDetails(CompoundsApp.FIELDS), []);
      const body = details!.querySelector<HTMLElement>('[data-u2-part="body"]')!;
      expect(body.hidden || getComputedStyle(body).display === 'none', true, 'Registration details are not folded');
      expect((app.ribbon().tools[0] as DomainSearch).tooltipText, CompoundsApp.SEARCH_HINT);
    } finally {
      view.close();
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
