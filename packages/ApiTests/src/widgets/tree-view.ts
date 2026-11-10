import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('TreeView', () => {
  test('Group expanding', async () => {
    const tree = ui.tree();
    const g = tree.group('Group 1', null, false);
    expect(g.expanded, false);
    g.expanded = true;
    expect(g.expanded, true);
  }, {owner: 'dkovalyov@datagrok.ai'});

  test('loadSources loads the query of the source it is given', async () => {
    const funcs = grok.dapi.functions;
    const source = funcs.filter('package.shortName = "ApiTests"').by(3);
    await funcs.filter('name like "Echo%"').list();
    const group = ui.tree().group('Functions', null, false);
    await group.loadSources(source);
    const loaded = group.children.map((n) => DG.toJs(n.value)?.id).filter((id) => id != null).sort();
    const expected = (await source.list()).map((f) => f.id).sort();
    expect(expected.length, 3);
    expect(loaded.join(), expected.join());
  }, {owner: 'drizhinashvili@datagrok.ai'});
});
