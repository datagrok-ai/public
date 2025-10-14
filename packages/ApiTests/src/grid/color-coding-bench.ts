import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/utils/src/test';


category('Grid: Color Coding Benchmarks', () => {
  test('Rendering', async () => {
    const t = grok.data.demo.demog(DG.Test.isInBenchmark ? 100: 1000);
    t.col('age')!.meta.colors.setConditional({'<30': DG.Color.green, '30-70': '#ff0000'});
    t.col('sex')!.meta.colors.setCategorical({'M': 0xFF0000FF, 'F': 0xFF800080});

    const v = grok.shell.addTableView(t);

    for (let i = 0; i < t.rowCount; i++)
      v.grid.invalidate();
  }, { benchmark: true, owner: 'dkovalyov@datagrok.ai'});
});
