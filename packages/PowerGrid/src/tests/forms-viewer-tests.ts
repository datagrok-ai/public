import * as grok from 'datagrok-api/grok';

import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';

category('FormsViewer', () => {
  test('Property change after the view is closed', async () => {
    const tv = grok.shell.addTableView(grok.data.demo.demog(100));
    const viewer = tv.addViewer('Forms');
    await awaitCheck(() => viewer.dataFrame != null, 'Forms viewer has not been attached to the table', 3000);
    tv.close();
    expect(viewer.dataFrame == null, true);
    viewer.props.rendererSize = 'large';
    viewer.props.showCurrentRow = false;
  });
});
