import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/utils/src/test';

category('TreeView', () => {
  test('Group expanding', async () => {
    const tree = ui.tree();
    const g = tree.group('Group 1', null, false);
    expect(g.expanded, false);
    g.expanded = true;
    expect(g.expanded, true);
  }, {owner: 'dkovalyov@datagrok.ai'});
});
