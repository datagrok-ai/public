import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {awaitCheck, before, category, delay, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {CruxSubstructureService} from '../crux-service/crux-substructure-service';
import {Fingerprint} from '../utils/chem-common';
import * as chemRdKit from '../utils/chem-common-rdkit';
import {getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {getSearchProgressEventName, getSearchQueryAndType, getTerminateEventName,
  SubstructureSearchType} from '../constants';
import {_package} from '../package-test';
import {SubstructureFilter} from '../widgets/chem-substructure-filter';

function columnOf(molecules: string[]): DG.Column {
  const column = DG.Column.fromList(DG.TYPE.STRING, 'smiles', molecules);
  DG.DataFrame.fromColumns([column]).name = 'Crux test';
  return column;
}

function indexes(bits: BitArray): number[] {
  const hits: number[] = [];
  for (let i = bits.findNext(-1); i !== -1; i = bits.findNext(i))
    hits.push(i);
  return hits;
}

function search(service: CruxSubstructureService, column: DG.Column, query: string, mask?: BitArray) {
  return service.search(column, query, query, Fingerprint.Morgan, 0.8, true, mask);
}

category('Crux substructure', () => {
  before(async () => {
    if (!chemRdKit.moduleInitialized) {
      chemRdKit.setRdKitWebRoot(_package.webRoot);
      await chemRdKit.initRdKitModuleLocal();
    }
  });

  test('original rows survive malformed and empty cells', async () => {
    const service = new CruxSubstructureService();
    try {
      const column = columnOf(['CCO', '', 'invalid', 'c1ccccc1', 'bad\ncell', 'Oc1ccccc1', 'CCN']);
      expectArray(indexes(await search(service, column, 'c1ccccc1')), [3, 5]);
      expectArray(indexes(await search(service, column, 'N')), [6]);
      expectArray(indexes(await search(service, column, '[Si]')), []);
    } finally {
      service.dispose();
    }
  });

  test('Contains agrees with RDKit on molecule and SMARTS queries', async () => {
    const service = new CruxSubstructureService();
    const molecules = ['CCO', 'CCN', 'CC(=O)N', 'CC(=O)O', 'c1ccccc1', 'Oc1ccccc1',
      'c1ccncc1', 'C1CCCCC1', 'C1CC1', '[NH4+]', '[Na+].[Cl-]', 'N#CC', 'CCOC', 'O=C=O'];
    const column = columnOf(molecules);
    const rdkit = chemRdKit.getRdKitModule();
    try {
      for (const query of ['CC', 'C1=CC=CC=C1', 'c1ccccc1', '[#7,#8]', '[N;H2]', '[R]', 'C(=O)N', '[N;+]', '[N+]']) {
        const qmol = getQueryMolSafe(query, query, rdkit)!;
        const expected: number[] = [];
        try {
          molecules.forEach((smiles, index) => {
            const mol = rdkit.get_mol(smiles)!;
            try {
              if (mol.get_substruct_match(qmol) !== '{}')
                expected.push(index);
            } finally {
              mol.delete();
            }
          });
        } finally {
          qmol.delete();
        }
        const actual = indexes(await search(service, column, query));
        expect(JSON.stringify(actual), JSON.stringify(expected),
          `${query}: Crux ${JSON.stringify(actual)}, RDKit ${JSON.stringify(expected)}`);
      }
    } finally {
      service.dispose();
    }
  });

  test('column edits and row changes invalidate the index', async () => {
    const service = new CruxSubstructureService();
    const column = columnOf(['CCO', 'CCN', 'CCC']);
    try {
      expectArray(indexes(await search(service, column, 'N')), [1]);
      column.set(0, 'CN');
      expectArray(indexes(await search(service, column, 'N')), [0, 1]);
      column.dataFrame.rows.removeAt(0, 1);
      expectArray(indexes(await search(service, column, 'N')), [0]);
      column.dataFrame.rows.addNew(['N']);
      expectArray(indexes(await search(service, column, 'N')), [0, 2]);
      const mask = new BitArray(column.length);
      mask.setBit(2, true);
      expectArray(indexes(await search(service, column, 'N', mask)), [2]);
    } finally {
      service.dispose();
    }
  });

  test('molblock targets and query', async () => {
    const service = new CruxSubstructureService();
    const rdkit = chemRdKit.getRdKitModule();
    const benzene = rdkit.get_mol('c1ccccc1')!;
    try {
      const block = benzene.get_molblock();
      const column = columnOf(['CCO', block, 'Oc1ccccc1']);
      expectArray(indexes(await service.search(column, block, benzene.get_smarts(),
        Fingerprint.Morgan, 0.8, true)), [1, 2]);
    } finally {
      benzene.delete();
      service.dispose();
    }
  });

  test('cancellation suppresses stale progress and results', async () => {
    const service = new CruxSubstructureService();
    const column = columnOf(Array.from({length: 50_001}, (_, i) => i % 2 ? 'CCN' : 'CCO'));
    let progress = 0;
    const subscription = grok.events.onCustomEvent(getSearchProgressEventName('Crux test', 'smiles'))
      .subscribe(() => progress++);
    try {
      const old = await service.search(column, 'O', 'O', Fingerprint.Morgan, 0.8);
      // Let index construction begin, then cancel while changing the query.
      await delay(10);
      grok.events.fireCustomEvent(getTerminateEventName('Crux test', 'smiles'),
        getSearchQueryAndType('O', SubstructureSearchType.CONTAINS, Fingerprint.Morgan, 0.8));
      const progressAtCancel = progress;
      const result = await search(service, column, 'N');
      expect(result.countBits(true), 25_000);
      expect(old.countBits(true), 0);
      expect(progress, progressAtCancel);
    } finally {
      subscription.unsubscribe();
      service.dispose();
    }
  });

  test('independent filters keep independent indexes', async () => {
    const first = new CruxSubstructureService();
    const second = new CruxSubstructureService();
    try {
      const results = await Promise.all([
        search(first, columnOf(['CCN', 'CCO']), 'N'),
        search(second, columnOf(['CCO', 'CCO', 'CCN']), 'N'),
      ]);
      expectArray(indexes(results[0]), [0]);
      expectArray(indexes(results[1]), [2]);
    } finally {
      first.dispose();
      second.dispose();
    }
  });

  test('million-row progress is bounded and complete', async () => {
    const source = await _package.files.readCsv('tests/smiles_50K.csv');
    const smiles = source.columns.byName('smiles').toList();
    const column = columnOf(Array.from({length: 1_000_000}, (_, i) => smiles[i % smiles.length]));
    const service = new CruxSubstructureService();
    const progress: number[] = [];
    const progressSubscription = grok.events.onCustomEvent(getSearchProgressEventName('Crux test', 'smiles'))
      .subscribe((value: number) => progress.push(value));
    let finishSubscription: {unsubscribe(): void} | undefined;
    const finished = new Promise<void>((resolve) => {
      finishSubscription = grok.events.onCustomEvent(getTerminateEventName('Crux test', 'smiles'))
        .subscribe(() => resolve());
    });
    try {
      const started = performance.now();
      const matches = await service.search(column, 'c1ccccc1', 'c1ccccc1', Fingerprint.Morgan, 0.8);
      await finished;
      const firstSearchMs = performance.now() - started;
      const warmStarted = performance.now();
      const warmMatches = await search(service, column, 'c1ccncc1');
      console.log(JSON.stringify({firstSearchMs, warmSearchMs: performance.now() - warmStarted,
        progress, hits: matches.countBits(true), warmHits: warmMatches.countBits(true)}));
      expect(matches.countBits(true), 903_980);
      expect(warmMatches.countBits(true), 150_220);
      expect(progress.length <= 7, true, `Expected at most seven updates, received ${progress.length}`);
      expect(progress[progress.length - 1], 100);
      expect(progress.filter((value) => value === 100).length, 1);
      expect(progress.every((value, i) => i === 0 || value > progress[i - 1]), true);
    } finally {
      progressSubscription.unsubscribe();
      finishSubscription?.unsubscribe();
      service.dispose();
    }
  }, {timeout: 120_000});

  test('filter applies consecutive cell edits and clears its query', async () => {
    const column = columnOf(['CCO', 'CCN', 'CCC']);
    column.semType = DG.SEMTYPE.MOLECULE;
    column.meta.units = DG.chem.Notation.Smiles;
    const filter = new SubstructureFilter();
    filter.attach(column.dataFrame);
    filter.applyState({columnName: column.name});
    const dialog = ui.dialog().add(filter.root).show();
    try {
      await awaitCheck(() => filter.sketcher.sketcher?.isInitialized === true, 'Sketcher initialization', 10_000);
      filter.sketcher.setSmiles('N');
      await awaitCheck(() => column.dataFrame.filter.trueCount === 1 && !filter.calculating,
        'Initial Crux filter', 10_000);
      column.set(0, 'CN');
      column.set(2, 'CCN');
      await awaitCheck(() => column.dataFrame.filter.trueCount === 3, 'Both edited rows must be rechecked', 10_000);
      filter.sketcher.setMolFile(DG.WHITE_MOLBLOCK);
      await awaitCheck(() => !filter.isFiltering && !filter.calculating && column.dataFrame.filter.trueCount === 3,
        'Cleared filter', 10_000);
    } finally {
      filter.detach();
      dialog.close();
    }
  });
});
