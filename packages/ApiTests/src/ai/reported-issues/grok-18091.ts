import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, findProp, look} from '../helpers';

// Regression coverage for GROK-18091: PC plot transformation persistence across filter-panel toggles.
category('AI: GROK-18091: PC plot transformation persistence', () => {
  test('transformation property exists and is reachable on PC plot', async () => {
    const v = DG.Viewer.pcPlot(demog(), {columnNames: ['age', 'height']});
    expect(v.type, DG.VIEWER.PC_PLOT);
    expect(findProp(v, 'transformation') != null, true);
  });

  test('empty-string transformation round-trip is non-throwing', async () => {
    const v = DG.Viewer.pcPlot(demog(20), {columnNames: ['age', 'height']});
    expectNoThrow(() => v.setOptions({transformation: ''}));
    const t = look(v)['transformation'];
    expect(t === '' || t == null, true);
  });
});
