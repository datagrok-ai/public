import {after, before, category, expect, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('Dapi: layouts', () => {

    test('Dapi: layouts - get applicable', async () => {
        let layouts = await grok.dapi.layouts.getApplicable(grok.data.demo.demog());
        if (layouts.length === 0)
            throw "Can't find any applicable layouts"
    });

});
