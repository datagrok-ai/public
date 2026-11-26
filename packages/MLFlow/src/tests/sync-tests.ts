import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {ListedModel, OwnRow, syncDiff} from '../sync';

category('Sync', () => {
  test('syncDiff registers new models, updates changed ones and deletes the gone ones', async () => {
    const listed = (connectionId: string, name: string, version: string, error = ''): ListedModel => ({connectionId,
      connection: connectionId, name, version, description: '', error, inputs: '[{"name": "a", "type": "double"}]',
      outputs: '[{"name": "y", "type": "float"}]'});
    const row = (id: string, connectionId: string, name: string, version: string, rowName = name): OwnRow =>
      ({id, name: rowName, description: '', target: 'y', features: 'a',
        options: {connectionId, name, version, type: 'double'}});
    const {register, remove} = syncDiff(
      [listed('c1', 'same', '1'), listed('c1', 'newer', '2'), listed('c1', 'new', '1'), listed('c2', '', '', 'down'),
        listed('c1', 'unsigned', '', 'no signature')],
      [row('r1', 'c1', 'same', '1', 'renamed by the user'), row('r2', 'c1', 'newer', '1'), row('r3', 'c1', 'gone', '1'),
        row('r4', 'c2', 'kept', '1'), row('r5', 'c3', 'unlisted', '1'), row('r6', 'c1', 'same', '1')]);
    expect(register.map((r) => `${r.model.name}:${r.row?.id ?? 'new'}`).join(), 'newer:r2,new:new');
    expect(remove.map((r) => r.id).sort().join(), 'r3,r5,r6');
  });
});
