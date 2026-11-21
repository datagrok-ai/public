import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {waitUntil, category, expect, test} from '@datagrok-libraries/test/src/test';

category('FormsViewer', () => {
  test('Property change after the view is closed', async () => {
    const tv = grok.shell.addTableView(grok.data.demo.demog(100));
    const viewer = tv.addViewer('Forms');
    await waitUntil(() => viewer.dataFrame != null, 'Forms viewer has not been attached to the table', 3000);
    tv.close();
    expect(viewer.dataFrame == null, true);
    viewer.props.rendererSize = 'large';
    viewer.props.showCurrentRow = false;
  });

  test('Colour-coded fields use the contrast text colour', async () => {
    const df = grok.data.demo.demog(100);
    df.col('race')!.meta.colors.setCategorical();
    const tv = grok.shell.addTableView(df);
    try {
      const viewer = tv.addViewer('Forms');
      const row = df.col('race')!.toList().indexOf('Other');
      df.currentRowIdx = row;
      let field: HTMLInputElement | undefined;
      await waitUntil(() => (field = Array.from(viewer.root.querySelectorAll('input'))
        .find((i) => i.value === 'Other' && i.style.backgroundColor !== '')) != null, 'colour-coded field not rendered', 5000);
      const expected = DG.Color.toHtml(DG.Color.getContrastColor(df.col('race')!.meta.colors.getColor(row)));
      expect(field!.style.getPropertyValue('color'), expected === '#ffffff' ? 'rgb(255, 255, 255)' : 'rgb(0, 0, 0)');
    }
    finally {
      tv.close();
    }
  });
});
