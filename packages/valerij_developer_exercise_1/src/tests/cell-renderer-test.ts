import {category, expect, test} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {NucleotideBoxCellRenderer, nucleotideBoxCellRenderer} from '../package';

category('NucleotideBoxCellRenderer', () => {
  test('cellType is dna_nucleotide', async () => {
    const renderer = new NucleotideBoxCellRenderer();
    expect(renderer.cellType, 'dna_nucleotide');
  });

  test('name is set', async () => {
    const renderer = new NucleotideBoxCellRenderer();
    expect(renderer.name, 'Nucleotide cell renderer');
  });

  test('factory function returns renderer instance', async () => {
    const renderer = nucleotideBoxCellRenderer();
    expect(renderer instanceof DG.GridCellRenderer, true);
  });

  test('renderer is registered for dna_nucleotide', async () => {
    const funcs = DG.Func.find({tags: ['cellRenderer']});
    const found = funcs.find((f: DG.Func) => f.options['cellType'] === 'dna_nucleotide');
    expect(found != null, true);
  });


  test('renders on a grid without errors', async () => {
    const df = DG.DataFrame.fromCsv(
      `sequence
GATTACA
ATTCGGA`
    );
    df.col('sequence')!.semType = 'dna_nucleotide';

    const tv = grok.shell.addTableView(df);
    try {
      // Give the grid time to render
      await new Promise((resolve) => setTimeout(resolve, 500));

      const grid = tv.grid;
      expect(grid != null, true);

      // Verify the grid's canvas exists and has content
      const canvas = grid.canvas;
      expect(canvas != null, true);
      expect(canvas.width > 0, true);
      expect(canvas.height > 0, true);
    } finally {
      tv.close();
    }
  });
});
