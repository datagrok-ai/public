import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, test, expect, before, after} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';

let group: DG.Group;

category('Packages: Properties', () => {
  before(async () => {
    if (!group)
      group = (await grok.dapi.groups.getGroupsLookup('All users'))[0];
  });

  test('INT', async () => {
    await changeProp('INT', 10);
  });

  test('FLOAT', async () => {
    await changeProp('FLOAT', 10.567);
  });

  test('BOOL', async () => {
    await changeProp('BOOL', false);
  });

  test('STRING', async () => {
    await changeProp('STRING', 'value2');
  });

  after(async () => {
    await _package.setProperties({INT: 1, FLOAT: 1.234, BOOL: true, STRING: 'value1'});
  });
});

async function changeProp(name: string, value: any): Promise<void> {
  //@ts-ignore
  await _package.setSettings({[name]: value}, group);
  const props = await _package.getSettings();
  //@ts-ignore
  expect(props[name], value);
}
