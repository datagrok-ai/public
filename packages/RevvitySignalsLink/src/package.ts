/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { u2 } from "@datagrok-libraries/utils/src/u2";
import '../css/revvity-signals-styles.css';
import { SignalsSearchParams, SignalsSearchQuery } from './signals-search-query';
import { queryEntities, queryLibraries, queryMaterialById, queryTags, queryTerms, queryUsers, RevvityApiResponse, RevvityData, RevvityUser } from './revvity-api';
import { dataFrameFromObjects, reorderColummns, transformData, createRevvityResponseWidget, getViewNameByCompoundType } from './utils';
import { addMoleculeStructures, assetsQuery, MOL_COL_NAME, retrieveQueriesMap } from './compounds';
import { getRevvityUsersWithMapping } from './users';
import { createInitialSatistics, getRevvityLibraries, RevvityLibrary, RevvityType } from './libraries';
import { createViewFromPreDefinedQuery, handleInitialURL } from './view-utils';
import { SAVED_SEARCH_STORAGE } from './search-utils';


export const _package = new DG.Package();

export type CurrentRevvityLibrary = {
  libId: string;
  libName: string;
  type?: string;
}

export type RevvityConfig = {
  libraries?: RevvityLibrary[];
}

let config: RevvityConfig = { libraries: undefined };

//tags: app
//name: Revvity Signals
//input: string path { meta.url: true; optional: true }
//output: view v
//meta.browsePath: Chem
export async function revvitySignalsLinkApp(path?: string): Promise<DG.ViewBase> {

  console.log(path);
  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/img/revvity.png',
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/RevvitySignalsLink/README.md',
    description: '- Integrate with your Revvity account.\n' +
      '- Browse the tenant content.\n' +
      '- Perfrom searches through you tenant.' +
      '- Find contextual information on entities like assets, batches etc.\n' +
      '- Analyze assay data.'
  });

  const statsDiv = ui.div('', {style: {position: 'relative'}});
  const view = DG.View.fromRoot(ui.divV([appHeader, statsDiv]));
  view.name = 'Revvity';

  if (path) {
    const cddNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Chem').getOrCreateGroup('Revvity Signals');
    cddNode.expanded = true;
    handleInitialURL(cddNode, path);
  } else
    createInitialSatistics(statsDiv);

  return view;

}

//input: dynamic treeNode
//input: view browseView
export async function revvitySignalsLinkAppTreeBrowser(treeNode: DG.TreeViewGroup, browseView: DG.View) {
  const loadingNode = treeNode.item('');
  loadingNode.captionLabel.classList.add('revvity-signals-starting-app');
  ui.setUpdateIndicator(loadingNode.captionLabel, true);
  let libs;
  try {
    await getRevvityUsersWithMapping();
    libs = await getRevvityLibraries();
  } catch(e: any) {
    loadingNode.remove();
    grok.shell.error(`Revvity libraries haven't been loaded: ${e?.message ?? e}`);
    return;
  }
  loadingNode.remove();
  for (const lib of libs) {
    const libNode = treeNode.group(lib.name);
    for (const libType of lib.types) {
      const viewName = getViewNameByCompoundType(libType.name);
      const typeNode = libNode.item(`${viewName.charAt(0).toUpperCase()}${viewName.slice(1)}`);
      typeNode.onSelected.subscribe(async () => {

        // //need woraround with nodeToDeselect to deselect node which was selected via routing (openRevvityNode function)
        // const nodeToDeselect = treeNode.items
        //   .find((node) => node.text.toLowerCase() !== viewName.toLowerCase() && node.root.classList.contains('d4-tree-view-node-selected'));
        // nodeToDeselect?.root.classList.remove('d4-tree-view-node-selected');

        await createViewFromPreDefinedQuery(treeNode, [lib.name, getViewNameByCompoundType(libType.name)], lib.name, libType.name);
      });
    }
  }
  const savedSearchesNode = treeNode.group('Saved searches');
  for (const lib of libs) {
    const libNode = savedSearchesNode.group(lib.name);
    for (const libType of lib.types) {
      const typeName = getViewNameByCompoundType(libType.name);
      const typeNode = libNode.group(`${typeName.charAt(0).toUpperCase()}${typeName.slice(1)}`);

      const storageKey = `${lib.name}|${libType.name}`;
      const savedSearchesStr = grok.userSettings.getValue(SAVED_SEARCH_STORAGE, storageKey) || '{}';
      const savedSearches: { [key: string]: string } = JSON.parse(savedSearchesStr);
      for (let key of Object.keys(savedSearches)) {
        const savedSearchNode = typeNode.item(key);
        savedSearchNode.onSelected.subscribe(async () => {
          await createViewFromPreDefinedQuery(treeNode,
            ['saved searches', lib.name, getViewNameByCompoundType(libType.name), key], lib.name, libType.name,
            JSON.parse(savedSearches[key]), true);
        });
      }
    }
  }

}

//name: Search Entities
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: string query
//input: string params
//output: dataframe df
export async function searchEntities(query: string, params: string): Promise<DG.DataFrame> {
  let df = DG.DataFrame.create();
  try {
    const queryJson: SignalsSearchQuery = JSON.parse(query);
    const paramsJson: SignalsSearchParams = JSON.parse(params);
    const response = await queryEntities(queryJson, Object.keys(paramsJson).length ? paramsJson : undefined);
    const data: Record<string, any>[] = !Array.isArray(response.data) ? [response.data!] : response.data!;
    df = dataFrameFromObjects(data);
    await grok.data.detectSemanticTypes(df);
  } catch (e: any) {
    grok.shell.error(e?.message ?? e);
  }
  return df;
}

//name: Search Entities With Structures
//input: string query
//input: string params
//output: dataframe df
export async function searchEntitiesWithStructures(query: string, params: string): Promise<DG.DataFrame> {
  let df = DG.DataFrame.create();
  try {
    const queryJson: SignalsSearchQuery = JSON.parse(query);
    const paramsJson: SignalsSearchParams = JSON.parse(params);
    const response = await queryEntities(queryJson, Object.keys(paramsJson).length ? paramsJson : undefined);
    if (!response.data || (Array.isArray(response.data) && response.data.length === 0))
      return DG.DataFrame.create();
    const data: Record<string, any>[] = !Array.isArray(response.data) ? [response.data!] : response.data!;

    const rows = await transformData(data);
    const moleculeIds = data.map((it) => it.id);
    df = DG.DataFrame.fromObjects(rows)!;
    const moleculeColumn = DG.Column.fromStrings(MOL_COL_NAME, new Array<string>(moleculeIds.length).fill(''));
    moleculeColumn.semType = DG.SEMTYPE.MOLECULE;
    moleculeColumn.meta.units = DG.UNITS.Molecule.MOLBLOCK;
    df.columns.add(moleculeColumn);
    reorderColummns(df);
    addMoleculeStructures(moleculeIds, moleculeColumn);
  } catch (e: any) {
    grok.shell.error(e?.message ?? e);
  }
  return df;
}

//name: Get Users
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//output: string users
export async function getUsers(): Promise<string> {
  const users: { [key: string]: RevvityUser } = {};
  const response = await queryUsers();
  if (!response.data || (Array.isArray(response.data) && response.data.length === 0))
    return '{}';
  const data: Record<string, any>[] = !Array.isArray(response.data) ? [response.data!] : response.data!;
  for (const user of data)
    users[user.id] = Object.assign({}, user.attributes || {});
  return JSON.stringify(users);
}


//name: Get Libraries
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//output: string libraries
export async function getLibraries(): Promise<string> {
  const response = await queryLibraries();
  const libraries: RevvityLibrary[] = [];
  if (!response.data || (Array.isArray(response.data) && response.data.length === 0))
    return '[]';
  const data: Record<string, any>[] = !Array.isArray(response.data) ? [response.data!] : response.data!;
  for (const lib of data) {
    if (lib.attributes.name && lib.id) {
      let types: RevvityType[] = [];
      const query = {
        "query": {
          "$and": [
            {
              "$match": {
                "field": "assetTypeEid",
                "value": `${lib.type}:${lib.id}`,
              }
            },
            {
              "$not": [
                {
                  "$match": {
                    "field": "type",
                    "value": "assetType"
                  }
                }
              ]
            }
          ]
        },
        "field": "type"
      }
      const typesResponse = await queryTerms(query);
      if (typesResponse.data) {
        const typesData: Record<string, any>[] = !Array.isArray(typesResponse.data) ? [typesResponse.data!] : typesResponse.data!;
        types = typesData.map((it) => {
          return {name: it.id, count: it.attributes?.count};
        });
      }
      libraries.push({ name: lib.attributes.name, id: `${lib.type}:${lib.id}`, types: types });
    }
  }
  return JSON.stringify(libraries);
}


//name: Get Tags
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: string type
//input: string assetTypeId
//output: string fields
export async function getTags(type: string, assetTypeId: string): Promise<string> {
  const query = {
    "query": {
      "$and": [
        {
          "$match": {
            "field": "assetTypeEid",
            "value": assetTypeId,
          }
        },
        {
          "$match": {
            "field": "type",
            "value": type,
            "mode": "keyword"
          }
        },
        {
          "$match": {
            "field": "isTemplate",
            "value": false
          }
        }
      ]
    }
  };
  const response = await queryTags(query);
  if (!response.data || (Array.isArray(response.data) && response.data.length === 0))
    return '{}';
  const data: Record<string, any>[] = !Array.isArray(response.data) ? [response.data!] : response.data!;
  const tags: { [key: string]: string } = {};
  for (const tag of data) {
    tags[tag.id] = tag.attributes.types[0].type;
  }
  return JSON.stringify(tags);
}


//name: Get Terms
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: string fieldName
//input: string type
//input: string assetTypeId
//input: bool isMaterial
//output: string terms
export async function getTerms(fieldName: string, type: string, assetTypeId: string, isMaterial: boolean): Promise<string> {
  const innerAndConditions: any[] = [
    {
      "$match": {
        "field": "assetTypeEid",
        "value": assetTypeId,
      }
    },
    {
      "$match": {
        "field": "type",
        "value": type,
        "mode": "keyword"
      }
    },
  ];
  if (isMaterial) {
    innerAndConditions.push({
      "$and": [
        {
          "$match": {
            "field": "isMaterial",
            "value": true
          }
        },
        {
          "$not": [
            {
              "$match": {
                "field": "type",
                "value": "assetType"
              }
            }
          ]
        }
      ]
    })
  }
  const query: SignalsSearchQuery = {
    "query": {
      "$and": innerAndConditions
    },
    field: fieldName,
    in: "tags"
  };
  const response = await queryTerms(query);
  if (!response.data || (Array.isArray(response.data) && response.data.length === 0))
    return '{}';
  const data: RevvityData[] = !Array.isArray(response.data) ? [response.data!] : response.data!;
  const tags = data.map((it) => it.id).filter((tag) => tag != undefined);
  return JSON.stringify(tags);
}


//name: Revvity Signals
//tags: panel, widgets
//input: semantic_value id { semType: RevvitySignalsId }
//output: widget result
export async function entityTreeWidget(idSemValue: DG.SemanticValue<string>): Promise<DG.Widget> {
  const obj = (await queryMaterialById(idSemValue.value)) as RevvityApiResponse;
  const div = createRevvityResponseWidget(obj, idSemValue);
  return new DG.Widget(div);
}
