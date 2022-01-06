// import {after, before, category, expect, test} from "../test";
// import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
// import * as DG from 'datagrok-api/dg';
//
// category('Dapi: properties', () => {
//
//     test('Dapi: properties - save, get, delete', async () => {
//         let group = DG.Group.create('js-api-test-group1');
//         group = await grok.dapi.groups.save(group);
//
//         let properties = {
//             'entityId': group.id,
//             'property': 'myProp',
//             'value': 'value'
//         };
//
//         await grok.dapi.entities.saveProperties([properties]);
//         await grok.dapi.entities.getProperties(group);
//         await grok.dapi.entities.deleteProperties([properties]);
//     });
//
// });
