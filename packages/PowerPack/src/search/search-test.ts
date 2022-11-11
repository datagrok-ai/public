import {after, before, category, expect, test} from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {initTemplates} from "./templates-search";

category('Search', () => {

  test('InitTemplates', async () => await initTemplates(), {skipReason : 'GROK-11445'});

});
