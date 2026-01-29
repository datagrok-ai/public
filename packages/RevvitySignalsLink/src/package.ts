/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/revvity-signals-styles.css';
import { SignalsSearchParams, SignalsSearchQuery } from './signals-search-query';
import { queryLibraries, queryStructureById, queryTags, queryTerms, queryUsers, RevvityData, RevvityUser, search } from './revvity-api';
import { reorderColumns, transformData, getViewNameByCompoundType, createRevvityWidgetByCorporateId, createWidgetByRevvityLabel } from './utils';
import { addMoleculeStructures, getConditionForLibAndType } from './compounds';
import { createInitialSatistics, getRevvityLibraries, refreshStats, resetRevvityLibraries, RevvityLibrary } from './libraries';
import { createViewForExpandabelNode, createViewFromPreDefinedQuery, handleInitialURL } from './view-utils';
import { createSavedSearchesSatistics, SAVED_SEARCH_STORAGE } from './search-utils';
import { funcs } from './package-api';
import { HIDDEN_ID_COL_NAME, ID_COL_NAME, MOL_COL_NAME, REVVITY_LABEL_SEM_TYPE, REVVITY_SEARCH_RES_TOTAL_COUNT, REVVVITY_LABEL_FIELDS, USER_FIELDS } from './constants';
import { getRevvityUsers, getUsersAllowed, updateRevvityUsers } from './users';
import { convertIdentifierFormatToRegexp } from './detectors';


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


//meta.role: init
export function init() {
  registerRevvityIdsFormats();
}

//meta.role: app
//name: Revvity Signals
//input: string path { meta.url: true; optional: true }
//output: view v
//meta.browsePath: Chem
export async function revvitySignalsLinkApp(path?: string): Promise<DG.ViewBase> {

  console.log(path);
  const initViewDiv = ui.divV([], {style: {flexGrow: '0'}});
  const view = DG.View.fromRoot(initViewDiv);
  view.name = 'Revvity';

  if (path) {
    const cddNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').getOrCreateGroup('Chem').getOrCreateGroup('Revvity Signals');
    cddNode.expanded = true;
    handleInitialURL(cddNode, path);
  } else {
    createInitialSatistics(initViewDiv).then((stats: HTMLDivElement) => {
      const refreshButton = ui.button('Refresh', async () => {
        const func = DG.Func.find({ package: 'RevvitySignalsLink', name: "getLibraries" });
        if (func.length)
          await grok.functions.clientCache.clear(func[0].id);
        resetRevvityLibraries();
        ui.empty(stats);
        refreshStats(stats);
      });
      refreshButton.style.width = '72px';
      initViewDiv.append(refreshButton);
    });
  }

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
    await getRevvityUsers();
    libs = await getRevvityLibraries();
  } catch(e: any) {
    loadingNode.remove();
    grok.shell.error(`Revvity libraries haven't been loaded: ${e?.message ?? e}`);
    return;
  }
  loadingNode.remove();

  for (const lib of libs) {
    const libNode = treeNode.group(lib.name);
    libNode.onSelected.subscribe(() => createViewForExpandabelNode(lib.name, createInitialSatistics, lib.name));
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

  const savedSearchesNode = treeNode.group('Saved Searches');
  savedSearchesNode.onSelected.subscribe(() => createViewForExpandabelNode('saved searches', createSavedSearchesSatistics));
  for (const lib of libs) {
    const libNode = savedSearchesNode.group(lib.name);
    libNode.onSelected.subscribe(() => createViewForExpandabelNode(lib.name, createSavedSearchesSatistics, lib.name));
    for (const libType of lib.types) {
      const typeName = getViewNameByCompoundType(libType.name);
      const typeNode = libNode.group(`${typeName.charAt(0).toUpperCase()}${typeName.slice(1)}`);
      typeNode.onSelected.subscribe(() => createViewForExpandabelNode(lib.name, createSavedSearchesSatistics, lib.name, libType.name));
      
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
//input: string query
//input: string params
//input: string libId
//input: string entityType
//output: dataframe df
export async function searchEntities(query: string, params: string, libId: string, entityType: string): Promise<DG.DataFrame> {
  const queryJson: SignalsSearchQuery = JSON.parse(query);
  const paramsJson: SignalsSearchParams = JSON.parse(params);
  const response = await search(queryJson, Object.keys(paramsJson).length ? paramsJson : undefined);
  const data: Record<string, any>[] = !Array.isArray(response.data) ? [response.data!] : response.data!;

  //in case /users endpoint is not allowed, extract users from included section
  if (!getUsersAllowed) {
    const newUsers = response.included?.filter((it) => it.type === 'user').map((it) => it.attributes);
    if (newUsers?.length)
      updateRevvityUsers(newUsers);
  }
  const df = await transformData(data, libId, entityType);

  //saving total result count
  df.setTag(REVVITY_SEARCH_RES_TOTAL_COUNT, response.meta ? response.meta!.total.toString() : '');

  USER_FIELDS.forEach((field) => {
    const col = df.col(field);
    if (col)
      col.semType = DG.TYPE.USER;
  });
  const idCol = df.col(ID_COL_NAME);
  if (idCol)
    idCol.name = HIDDEN_ID_COL_NAME;
  for (const colName of REVVVITY_LABEL_FIELDS) {
    const col = df.col(colName);
    if (col)
      col.semType = REVVITY_LABEL_SEM_TYPE;
  }
  return df;
}

//name: Search Entities With Structures
//input: string query
//input: string params
//input: string libId
//input: string entityType
//input: bool doNotAddStructures {optional: true}
//output: dataframe df
export async function searchEntitiesWithStructures(query: string, params: string,
  libId: string, entityType: string, doNotAddStructures?: boolean): Promise<DG.DataFrame> {
  let df = DG.DataFrame.create();
  try {
    df = await funcs.searchEntities(query, params, libId, entityType);
    let idCol = df.col(HIDDEN_ID_COL_NAME);
    if (!idCol)
      idCol = df.col(ID_COL_NAME);
    if (idCol) {
      const moleculeIds = idCol.toList();
      const moleculeColumn = DG.Column.fromStrings(MOL_COL_NAME, new Array<string>(moleculeIds.length).fill(''));
      moleculeColumn.semType = DG.SEMTYPE.MOLECULE;
      moleculeColumn.meta.units = DG.UNITS.Molecule.MOLBLOCK;
      df.columns.add(moleculeColumn);
      reorderColumns(df);
      if (!doNotAddStructures)
        addMoleculeStructures(moleculeIds, moleculeColumn);
    }
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
  const users: RevvityUser[] = [];
  const response = await queryUsers();
  if (!response.data || (Array.isArray(response.data) && response.data.length === 0))
    return '[]';
  const data: Record<string, any>[] = !Array.isArray(response.data) ? [response.data!] : response.data!;
  for (const user of data)
    users.push(user.attributes);
  return JSON.stringify(users);
}


//name: Get Libraries
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//output: string libraries
export async function getLibraries(): Promise<string> {
  const response = await queryLibraries();
  if (!response.data || (Array.isArray(response.data) && response.data.length === 0))
    return '[]';
  const data: Record<string, any>[] = !Array.isArray(response.data) ? [response.data!] : response.data!;
  return JSON.stringify(data);
}


//name Register Revvity Ids Formats
export async function registerRevvityIdsFormats() {
  const data = JSON.parse(await funcs.getLibraries());
  const formats: string[] = [];
  for (const lib of data) {
    if (lib.attributes?.assets?.numbering?.format) {
      formats.push(lib.attributes?.assets?.numbering?.format);
      if (lib.attributes?.batches?.numbering?.format)
        formats.push(`${formats[0]}-${lib.attributes?.batches?.numbering?.format}`);
    }
  }
  const regexp = convertIdentifierFormatToRegexp(Object.values(formats));
  if (regexp)
    DG.SemanticValue.registerRegExpDetector('revvity-id', regexp);
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
    if (tag?.attributes?.types && Array.isArray(tag?.attributes?.types) && tag?.attributes?.types.length > 0)
    tags[tag.id] = tag?.attributes?.types[0].type;
  }
  return JSON.stringify(tags);
}

//name: Search Terms
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: string query
//output: string terms
export async function searchTerms(query: string): Promise<string> {
  const response = await queryTerms(JSON.parse(query));
  if (!response.data || (Array.isArray(response.data) && response.data.length === 0))
    return '[]';
  const data: RevvityData[] = !Array.isArray(response.data) ? [response.data!] : response.data!;
  return JSON.stringify(data);
}


//name: Get Terms For Field
//input: string fieldName
//input: string type
//input: string assetTypeId
//input: bool isMaterial
//output: list<string> terms
export async function getTermsForField(fieldName: string, type: string, assetTypeId: string, isMaterial: boolean): Promise<string[]> {
  const innerAndConditions = getConditionForLibAndType(type, assetTypeId, isMaterial);
  const query: SignalsSearchQuery = {
    "query": {
      "$and": innerAndConditions
    },
    field: fieldName,
    in: "tags"
  };
  const data: RevvityData[] = JSON.parse(await funcs.searchTerms(JSON.stringify(query)));
  const terms = data.map((it) => it.id).filter((term) => term != undefined);
  return terms;
}


//name: Get Structure By Id
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: string id
//output: string structure
export async function getStructureById(id: string): Promise<string> {
  const molecule = await queryStructureById(id) as string;
  return molecule;
}


//name: Revvity Signals
//meta.role: widgets,panel
//input: semantic_value id { semType: revvity-id }
//output: widget result
export async function entityTreeWidget(idSemValue: DG.SemanticValue<string>): Promise<DG.Widget> {
  const query = {
    "query": {
      "$match": {
        "field": "name",
        "value": idSemValue.value,
        "mode": "keyword"
      }
    }
  }
  const obj = await search(query);
  const div = createRevvityWidgetByCorporateId(obj, idSemValue);
  return new DG.Widget(div);
}

//name: Revvity Signals
//meta.role: widgets,panel
//input: semantic_value id { semType: revvity-label }
//output: widget result
export async function revvityLabelWidget(idSemValue: DG.SemanticValue<string>): Promise<DG.Widget> {
  const div = await createWidgetByRevvityLabel(idSemValue);
  return new DG.Widget(div);
}