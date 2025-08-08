/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from "@datagrok-libraries/utils/src/u2";
import { buildOperatorUI, createDefaultOperator, signalsSearchBuilderUI } from './signals-search-builder';
import '../css/revvity-signals-styles.css';
import { SignalsSearchParams, SignalsSearchQuery } from './signals-search-query';
import { queryEntities, queryEntityById, queryLibraries, queryMaterialById, queryStructureById, queryTerms, queryUsers, RevvityApiResponse, RevvityData, RevvityUser } from './revvity-api';
import { dataFrameFromObjects, reorderColummns, transformData, widgetFromObject, createRevvityResponseWidget } from './utils';
import { addMoleculeStructures, assetsQuery, batchesQuery, MOL_COL_NAME } from './compounds';
import { RevvityFilters } from './filters';
import { getDefaultProperties } from './properties';
import { buildPropertyBasedQueryBuilder, QueryBuilder } from './query-builder';
import { testFilterCondition } from './const';


export const _package = new DG.Package();
export type RevvityLibrary = {
  id: string;
  name: string;
  types: string[];
}
export type RevvityConfig = {
  libraries?: RevvityLibrary[];
}
let openedView: DG.View | null = null;
let config: RevvityConfig = {libraries: undefined};

//tags: app
//name: Revvity Signals
//output: view v
//meta.browsePath: Chem
export async function revvitySignalsLinkApp(): Promise<DG.ViewBase> {

  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/images/benchling.png',
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/RevvitySignalsLink/README.md',
    description: '- Integrate with your Revvity account.\n' +
      '- Analyze assay data.\n' +
      '- Browse the tenant content.\n'
  });

  const view = DG.View.fromRoot(appHeader);
  view.name = 'Revvity';
  return view;
}


//input: dynamic treeNode
//input: view browseView
export async function revvitySignalsLinkAppTreeBrowser(treeNode: DG.TreeViewGroup) {
  const search = treeNode.item('Search');
  search.onSelected.subscribe(() => {
    const v = DG.View.create('Search');
    const queryBuilder = signalsSearchBuilderUI();
    v.append(queryBuilder);
    grok.shell.addPreview(v);
  });

  const search2 = treeNode.item('Search 2');
  search2.onSelected.subscribe(() => {
    const v = DG.View.create('Search 2');
    const queryBuilder = buildPropertyBasedQueryBuilder(getDefaultProperties(), JSON.parse(JSON.stringify(testFilterCondition)));
    v.append(queryBuilder);
    grok.shell.addPreview(v);
  });

  const createRadio = (form: HTMLElement, types: string[], searchObj: { library: string, type: string | undefined }) => {
    const previousRadio = form.querySelector('.revvity-seach-type-radio');
    if (previousRadio)
      previousRadio.remove();
    if (types.length) {
      const radio = ui.input.radio('Type', {
        value: types[0],
        items: types,
        onValueChanged: () => {
          searchObj.type = radio.value!;
        }
      })
      radio.classList.add('revvity-seach-type-radio');
      form.append(radio.root);
    }
  };

  const search3 = treeNode.item('Materials Search');
  search3.onSelected.subscribe(async () => {
    const v = DG.View.create({ name: 'Materials Search' });
    ui.setUpdateIndicator(v.root, true, 'Loading Revvity filters...');
    if (!config.libraries)
      getLibraries().then((libs: string) => {
        config.libraries = JSON.parse(libs);
        createMaterialsSearchUI(v);
      })
    else
      createMaterialsSearchUI(v);
    grok.shell.addPreview(v);
  });


  const createMaterialsSearchUI = (v: DG.View) => {
    const searchDiv = ui.divV([], 'revvity-signals-materials-search');
    if (config.libraries?.length) {
      const radioDiv = ui.div();
      const searchObj = { library: config.libraries[0].name, type: config.libraries[0].types.length ? config.libraries[0].types[0] : undefined};
      const library = ui.input.choice('Search', {
        value: config.libraries[0].name,
        items: config.libraries.map((it) => it.name),
        onValueChanged: () => {
          searchObj.library = library.value!;
          const selectedLib = config.libraries!.filter((it) => it.name === library.value!);
          searchObj.type = selectedLib[0].types.length ? selectedLib[0].types[0] : undefined;
          createRadio(radioDiv, selectedLib[0].types, searchObj);
        }
      });
      const libraryTypeForm = ui.form([library]);
      createRadio(libraryTypeForm, config.libraries[0].types, searchObj);
      searchDiv.append(libraryTypeForm);
    }
    const queryBuilder = new QueryBuilder(getDefaultProperties());
    searchDiv.append(queryBuilder.root);
    const runButton = ui.bigButton('RUN', () => {
      console.log(queryBuilder.condition);
    });
    runButton.classList.add('revvity-signals-materials-search-run');
    searchDiv.append(runButton);
    v.append(searchDiv);
    ui.setUpdateIndicator(v.root, false);
  }

  const createViewFromPreDefinedQuery = async (query: string, name: string) => {
    const df = await grok.functions.call('RevvitySignalsLink:searchEntitiesWithStructures', {
      query: query,
      params: '{}'
    });
    const tv = grok.shell.addTablePreview(df);
    tv.name = name;
    new RevvityFilters(tv);
  }

  const compounds = treeNode.group('Compounds');

  const assets = compounds.item('Assets');
  assets.onSelected.subscribe(async () => {
    await createViewFromPreDefinedQuery(JSON.stringify(assetsQuery), 'Assets');
  });

  const batches = compounds.item('Batches');
  batches.onSelected.subscribe(async () => {
    await createViewFromPreDefinedQuery(JSON.stringify(batchesQuery), 'Batches');
  });
}

//name: Search Entities
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
    if (!response.data)
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
//output: string users
export async function getUsers(): Promise<string> {
  const users: {[key: string]: RevvityUser} = {};
  const response = await queryUsers();
  if (!response.data)
    return '{}';
  const data: Record<string, any>[] = !Array.isArray(response.data) ? [response.data!] : response.data!;
  for (const user of data)
    users[user.id] = Object.assign({}, user.attributes || {});
  return JSON.stringify(users);
}


//name: Get Libraries
//output: string libraries
export async function getLibraries(): Promise<string> {
  const response = await queryLibraries();
  const libraries: RevvityLibrary[] = [];
  if (!response.data)
    return '[]';
  const data: Record<string, any>[] = !Array.isArray(response.data) ? [response.data!] : response.data!;
  for (const lib of data) {
    if (lib.attributes.name && lib.id) {
      let types: string[] = [];
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
        types = typesData.map((it) => it.id);
      }
      libraries.push({name: lib.attributes.name, id: `${lib.type}:${lib.id}`, types: types});
    }
  }
  return JSON.stringify(libraries);
}


//name: Revvity Signals
//tags: panel, widgets
//input: string id { semType: RevvitySignalsId }
//output: widget result
export async function entityTreeWidget(id: string): Promise<DG.Widget> {
  const obj = (await queryMaterialById(id)) as RevvityApiResponse;
  const div = createRevvityResponseWidget(obj);
  return new DG.Widget(div);
}