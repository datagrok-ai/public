import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package';

export async function initTestStickyMeta() {
    let schemas: DG.Schema[] = await grok.dapi.stickyMeta.getSchemas();
    let schema: DG.Schema | undefined = schemas.filter((schema) => schema.name == 'Autotests').at(0);
    if (schema !== undefined)
        return;
    await grok.dapi.stickyMeta.createSchema('Autotests', [{
        name: 'autotest', matchBy: 'semtype=test'
    }], [{
        name: 'tickets', type: DG.TYPE.STRING 
    }]);
}