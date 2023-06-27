import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('Dapi: entities', () => { 
  let group: DG.Group;

  before(async () => {
    group = DG.Group.create('js-api-test-group1');
    group = await grok.dapi.groups.save(group);
    const properties = {
      entityId: group.id,
      property: 'myProp',
      value: 'value',
    };
    await group.setProperties(properties);
  });

  test('getProperties', async () => {
    const props = await group.getProperties();
    expect(typeof props === 'object', true);
    expect(Object.keys(props).length, 3);
  });

  test('setProperties', async () => {
    await group.setProperties({testProp1: 'prop1', testProp2: 'prop2'});
    expect(Object.keys(await group.getProperties()).length, 5);
  });

  after(async () => {
    await grok.dapi.groups.delete(group);
  });
});
