import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/utils/src/test';

category('Widgets', () => {
  test('Legend', async () => {
    const df = grok.data.demo.demog();
    const column = df.getCol('disease');
    const legend = DG.Legend.create(column);
    expect(legend instanceof DG.Legend, true);
    expect(legend.column, column);
    expect(legend.root instanceof HTMLElement, true);
    expect(legend.root.classList.contains('d4-legend'), true);
    const newColumn = df.getCol('site');
    legend.column = newColumn;
    expect(legend.column, newColumn);
    legend.showNulls = false;
    expect(legend.showNulls, false);
    legend.position = DG.LEGEND_POSITION.BOTTOM;
    expect(legend.position, DG.LEGEND_POSITION.BOTTOM);
    expect(legend.root.classList.contains(`d4-legend-${DG.LEGEND_POSITION.BOTTOM}`), true);
  });
}, {owner: 'agolovko@datagrok.ai'});
