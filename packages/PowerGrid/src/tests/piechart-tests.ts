import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {PieChartCellRenderer} from '../sparklines/piechart';

category('PieChart', () => {
  test('Legacy subsector without functionType', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [0.5, 2])]);
    const tv = grok.shell.addTableView(df);
    try {
      const gc = tv.grid.columns.add({gridColumnName: 'pie', cellType: 'piechart'});
      gc.settings = {
        columnNames: ['x'],
        sectors: {lowerBound: 0, upperBound: 1, values: '', sectors: [{name: 'Group', sectorColor: '#ff0000',
          subsectors: [{name: 'x', weight: 1, line: [[0, 0], [1, 1], [3, 1]], min: 0, max: 3}]}]},
      };
      const g = ui.canvas(80, 80).getContext('2d')!;
      const cell = tv.grid.cell('pie', 0);
      new PieChartCellRenderer().render(g, 0, 0, 80, 80, cell, cell.style);
      expect(gc.settings.sectors.sectors[0].subsectors[0].functionType, 'numerical');
    } finally {
      tv.close();
    }
  });
});
