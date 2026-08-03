import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import { awaitCheck, category, expect, test } from '@datagrok-libraries/test/src/test';
import { openRevvityNode } from '../view-utils';
import { Operators } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { funcs } from '../package-api';
import { MOL_COL_NAME } from '../constants';
import { getLibrariesWithEntityTypes } from '../libraries';

category('revvity signals app', () => {

  test('app initial statistics', async () => {
    const view = await funcs.revvitySignalsLinkApp();
    //check that statistics view has been created
    await awaitCheck(() => view.root.getElementsByTagName('table').length > 0, 'Initial statistics hasn\'t been loaded', 30000);
  });

  test('open compounds|assets node', async () => {
    const node = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Chem').getOrCreateGroup('Revvity Signals');
    node.expanded = true;
    openRevvityNode(node, ['Compounds'], 'Assets', 'Compounds', 'asset');
    await awaitCheck(() => grok.shell.tv?.dataFrame.rowCount === 100, 'Compounds -> Assets node hasn\'t been loaded', 30000);
  }, {timeout: 60000});

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
    openRevvityNode(node, ['Compounds'], 'Assets', 'Compounds', 'asset', query);
    await awaitCheck(() => grok.shell.tv?.dataFrame.rowCount === 100, 'Compounds -> Assets node with pre-defined query hasn\'t been loaded', 30000);
  }, {timeout: 60000});
});


category('revvity signals functions', () => {

  test('searchEntities', async () => {
    const query = {
      "query": {
        "$match": {
          "field": "type",
          "value": "batch",
          "mode": "keyword"
        }
      }
    };
    const df = await funcs.searchEntities(JSON.stringify(query), '{}', 'assetType:686ecf60e3c7095c954bd94f', 'batch');
    expect(df.rowCount > 0, true, 'Returned empty dataframe');
  }, {timeout: 60000});

  test('searchEntitiesWithStructures', async () => {
    const query = {
      "query": {
        "$match": {
          "field": "type",
          "value": "asset",
          "mode": "keyword"
        }
      }
    };
    const df = await funcs.searchEntitiesWithStructures(JSON.stringify(query), '{}', 'assetType:686ecf60e3c7095c954bd94f', 'asset');
    expect(df.rowCount > 0, true, 'Returned empty dataframe');
    await awaitCheck(() => !df.col(MOL_COL_NAME)!.isEmpty, 'Molecules column is not filling with values', 30000);
  }, {timeout: 60000});

  test('getLibrariesWithEntityTypes', async () => {
    const libs = await getLibrariesWithEntityTypes();
    expect(libs.length > 0, true, 'Returned empty libraries list');
  }, {timeout: 60000});

  test('getUsers', async () => {
    const users = JSON.parse(await funcs.getUsers());
    expect(users.length > 0, true, 'Returned empty users list');
  }, {timeout: 60000});

  test('getTags', async () => {
    const tags = JSON.parse(await funcs.getTags('batch', 'assetType:686ecf60e3c7095c954bd94f'));
    expect(Object.keys(tags).length > 0, true, 'Returned empty tags list');
  }, {timeout: 60000});

  test('getTerms', async () => {
    const terms = await funcs.getTermsForField('materials.Batch Chemical Name', 'batch', 'assetType:686ecf60e3c7095c954bd94f', true);
    expect(terms.length > 0, true, 'Returned empty terms list');
  }, {timeout: 60000});

  test('getWidget', async () => {
    const widget: DG.Widget = await funcs.entityTreeWidget(DG.SemanticValue.fromValueType('DGS-0000009-001', 'revvity-id'));
    await awaitCheck(() => widget.root.querySelector('table') != null, 'Widget hasn\'t been created', 30000);
  }, {timeout: 60000});

});
