import {after, before, category, expect, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('Dapi: projects', () => {

    test('Dapi: projects - open', async () => {
        await grok.dapi.projects.open('demog');
    });

});
