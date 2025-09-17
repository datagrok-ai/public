import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {awaitCheck, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
//import { funcs } from '../package-api';
import { openRevvityNode } from '../view-utils';
import { Operators } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { funcs } from '../package-api';


category('revvity signals', () => {

  before(async () => {
  });
  
  test('app tests', async () => {
    const view = await funcs.revvitySignalsLinkApp();
    //check that statistics view has been created
    await awaitCheck(() => view.root.getElementsByTagName('table') != null, 'Initial statistics hasn\'t been loaded', 30000);
    grok.shell.closeAll();
  });

  test('open compounds|assets node', async () => {
    const node = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Chem').getOrCreateGroup('Revvity Signals');
    node.expanded = true;
    openRevvityNode(node, ['Compounds'], 'asset', 'Compounds', 'asset');
    await awaitCheck(() => grok.shell.tv?.dataFrame.rowCount === 20, 'Compounds -> Assests node hasn\'t been loaded', 30000);
    grok.shell.closeAll();
  });

    test('open compounds|assets with pre-defined query', async () => {
        const query = {
            "logicalOperator": Operators.Logical.and,
            "conditions": [
                {
                    "field": "createdAt",
                    "operator": "before",
                    "value": "2025-09-06T21:00:00.000Z"
                }
            ]
        }
    const node = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Chem').getOrCreateGroup('Revvity Signals');
    node.expanded = true;
    openRevvityNode(node, ['Compounds'], 'asset', 'Compounds', 'asset', query);
    await awaitCheck(() => grok.shell.tv?.dataFrame.rowCount === 20, 'Compounds -> Assests node with pre-defined query hasn\'t been loaded', 30000);
    grok.shell.closeAll();
  });
});


