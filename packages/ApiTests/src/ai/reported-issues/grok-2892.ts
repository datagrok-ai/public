import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, withTableView} from '../helpers';

// Regression coverage for GROK-2892: Pie chart onClick filter semantics (property round-trip).
category('AI: GROK-2892: Pie chart onClick filter semantics', () => {
  test('factory + setOptions: onClick=Filter round-trips through props and look', async () => {
    await withTableView(demog(), async (tv) => {
      const v = DG.Viewer.pieChart(demog(), {category: 'race', onClick: 'Filter'});
      tv.addViewer(v);
      expect(v.props.onClick, 'Filter');
      v.setOptions({onClick: 'Filter'});
      expect(v.props.onClick, 'Filter');
      expectLook(v, {onClick: 'Filter'});
    });
  });

  test('RowGroupAction: all three values round-trip via setOptions + props + look', async () => {
    await withTableView(demog(), async (tv) => {
      const v = DG.Viewer.pieChart(demog(), {category: 'race'});
      tv.addViewer(v);
      for (const action of ['Select', 'Filter', 'None']) {
        v.setOptions({onClick: action});
        expect(v.props.onClick, action);
        expectLook(v, {onClick: action});
      }
    });
  });

  test('getProperties(): onClick is exposed (best-effort description check)', async () => {
    await withTableView(demog(), async (tv) => {
      const v = DG.Viewer.pieChart(demog(), {category: 'race', onClick: 'Filter'});
      tv.addViewer(v);
      const props = v.getProperties();
      const onClickProp = props.find((p) => p.name === 'onClick');
      expect(onClickProp != null, true);
      if (onClickProp != null && onClickProp.description != null)
        expect(typeof onClickProp.description === 'string' && onClickProp.description.length > 0, true);
    });
  });
});
