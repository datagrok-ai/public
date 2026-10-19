import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {_fetchENASequence, complement, countSubsequenceTableAugment, enaSequence, getOrders,
  openTable1, openTable2, openTable3} from '../package';

category('Sequence', () => {
  test('complement', async () => {
    expect(complement('fasta: ATTGCA'), 'fasta: TAACGT');
    expect(complement('gttc ca'), 'caag gt');
  });

  test('nucleotide renderer', async () => {
    const df = DG.DataFrame.fromCsv('sequence\nfasta: ATGCATGC\nfasta: GGCCTAAT');
    df.col('sequence')!.semType = 'dna_nucleotide';
    const tv = grok.shell.addTableView(df);
    // The grid starts with the string renderer and switches once the package renderer loads.
    // It returns a proxy whose name is the cell type, so check the cell type it resolved.
    await awaitCheck(() => tv.grid.col('sequence')!.renderer.cellType === 'dna_nucleotide',
      'Nucleotide renderer was not picked up', 5000);
    tv.close();
  });

  test('countSubsequenceTableAugment', async () => {
    const df = DG.DataFrame.fromCsv('Sequence\naccacc\ngttc\naccc');
    await countSubsequenceTableAugment(df, df.col('Sequence')!, 'acc');
    const counts = df.col('N(acc)')!;
    expect([counts.get(0), counts.get(1), counts.get(2)].join(','), '2,0,1');
  }, {timeout: 120000});
});

category('ENA', () => {
  test('enaSequence', async () => {
    const widget = await enaSequence('MT079845.1');
    expect(widget.root.textContent!.includes('>ENA|MT079845|MT079845.1'), true);
    expect(widget.root.querySelector('textarea')!.value.startsWith('CTACAAGAGATCGAAAGTTGG'), true);
  }, {timeout: 60000});

  test('_fetchENASequence', async () => {
    const df = await _fetchENASequence('coronavirus', 3, 20);
    expect(df.rowCount, 3);
    expect(df.columns.names().join(','), 'ID,Sequence');
    expect(df.col('Sequence')!.semType, 'dna_nucleotide');
    for (let i = 0; i < df.rowCount; i++) {
      expect(/^[A-Z]{2}\d{6}$/.test(df.get('ID', i)), true);
      expect(df.get('Sequence', i).length <= 20, true);
    }
  }, {timeout: 60000});
});

category('Queries', () => {
  test('getOrders', async () => {
    const df = await getOrders();
    expect(df.rowCount, 13);
    expect(df.columns.names().join(','), 'customerid,sum');
  });
});

category('Files', () => {
  test('openTable variants', async () => {
    const tables = [
      await openTable1('bio/sars-cov-2.csv'),
      await openTable2('bio/sars-cov-2.csv'),
      await openTable3('bio/sars-cov-2.csv'),
    ];
    for (const df of tables) {
      expect(df.rowCount > 0, true);
      expect(df.rowCount, tables[0].rowCount);
      expect(df.columns.names().join(','), 'GenBank,Sequence');
      await grok.data.detectSemanticTypes(df);
      expect(df.col('Sequence')!.semType, 'dna_nucleotide');
      expect(df.col('GenBank')!.semType, 'EnaID');
    }
  });

  test('addTables', async () => {
    const before = Array.from(grok.shell.tables).length;
    // Called through the platform: in the test bundle the imported _package is not initialized,
    // so its files root would be System:AppData itself.
    await grok.functions.call('GuillermoSequence:addTables');
    // cars.csv, demog.csv, test-ena-ids.csv and test-sequences.csv; aspirin.mol is not a CSV.
    expect(Array.from(grok.shell.tables).length - before, 4);
  });
}, {clear: true});
