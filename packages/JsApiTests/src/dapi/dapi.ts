import {after, before, category, expect, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('Dapi', () => {

    test('Dapi: all data sources', async () => {
        await grok.dapi.queries.first();
        await grok.dapi.connections.first();
        await grok.dapi.credentials.first();
        await grok.dapi.jobs.first();
        await grok.dapi.notebooks.first();
        await grok.dapi.models.first();
        await grok.dapi.packages.first();
        await grok.dapi.layouts.first();
        await grok.dapi.tables.first();
        await grok.dapi.users.first();
        await grok.dapi.groups.first();
        await grok.dapi.scripts.first();
        await grok.dapi.projects.first();
        await grok.dapi.environments.first();
    });

});
