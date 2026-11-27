import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {syncDiff} from '../sync';
import {ListedModel, OwnRow} from '../types';

category('Sync', () => {
  test('syncDiff registers new models, updates changed ones and deletes the gone ones', async () => {
    const listed = (connection: string, name?: string, version?: string, error?: string): ListedModel =>
      ({connection_id: connection, connection, name, version, error, inputs: '[{"name": "a", "type": "double"}]',
        outputs: '[{"name": "y", "type": "float"}]'});
    const row = (id: string, connectionId: string, name: string, version: string, rowName = name): OwnRow =>
      ({Id: id, Name: rowName, Target: 'y', Features: 'a',
        Options: JSON.stringify({connectionId, name, version, type: 'double'})});
    const {register, remove} = syncDiff(
      [listed('c1', 'same', '1'), listed('c1', 'newer', '2'), listed('c1', 'new', '1'),
        listed('c2', undefined, undefined, 'down'), listed('c1', 'unsigned', undefined, 'no signature')],
      [row('r1', 'c1', 'same', '1', 'renamed by the user'), row('r2', 'c1', 'newer', '1'), row('r3', 'c1', 'gone', '1'),
        row('r4', 'c2', 'kept', '1'), row('r5', 'c3', 'unlisted', '1'), row('r6', 'c1', 'same', '1')]);
    expect(register.map((r) => `${r.model.name}:${r.row?.Id ?? 'new'}`).join(), 'newer:r2,new:new');
    expect(remove.map((r) => r.Id).sort().join(), 'r3,r5,r6');
  });
});
