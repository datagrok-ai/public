/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { u2 } from "@datagrok-libraries/utils/src/u2";
import '../css/revvity-signals-styles.css';
import { convertComplexConditionToSignalsSearchQuery, SignalsSearchParams, SignalsSearchQuery } from './signals-search-query';
import { queryEntities, queryLibraries, queryMaterialById, queryTags, queryTerms, queryUsers, RevvityApiResponse, RevvityData, RevvityUser } from './revvity-api';
import { dataFrameFromObjects, reorderColummns, transformData, createRevvityResponseWidget } from './utils';
import { addMoleculeStructures, assetsQuery, batchesQuery, materialsCondition, MOL_COL_NAME } from './compounds';
import { getDefaultProperties, REVVITY_FIELD_TO_PROP_TYPE_MAPPING } from './properties';
import { ComplexCondition, Operators, QueryBuilder, QueryBuilderLayout, SUGGESTIONS_FUNCTION } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { getRevvityUsers } from './users';
import { createInitialSatistics, getRevvityLibraries, RevvityLibrary, RevvityType } from './libraries';



export const _package = new DG.Package();

export type CurrentRevvityLibrary = {
  libId: string;
  libName: string;
  type?: string;
}

export type RevvityConfig = {
  libraries?: RevvityLibrary[];
}

let openedView: DG.View | null = null;
let config: RevvityConfig = { libraries: undefined };

//tags: app
//name: Revvity Signals
//output: view v
//meta.browsePath: Chem
export async function revvitySignalsLinkApp(): Promise<DG.ViewBase> {

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
  createInitialSatistics(statsDiv);
  return view;

}

//input: dynamic treeNode
export async function revvitySignalsLinkAppTreeBrowser(treeNode: DG.TreeViewGroup) {
  getRevvityUsers();
  const libs = await getRevvityLibraries();

  const createViewFromPreDefinedQuery = async (query: string, name: string, libName: string, compoundType: string) => {
    //!!!!!!!!!!TODO: Change to .then()
    const df = await grok.functions.call('RevvitySignalsLink:searchEntitiesWithStructures', {
      query: query,
      params: '{}'
    });
    const tv = grok.shell.addTablePreview(df);
    tv.name = name;
    const filtersDiv = ui.div([]);
    let queryBuilder: QueryBuilder | null = null;

    const updateQueryBuilderLayout = (width: number) => {
      if (!queryBuilder) return;

      // Switch to narrow layout if width is less than 300px, otherwise use standard
      const newLayout = width < 300 ? QueryBuilderLayout.Narrow : QueryBuilderLayout.Standard;

      if (queryBuilder.getLayout() !== newLayout) {
        queryBuilder.setLayout(newLayout);
      }
    };

    ui.onSizeChanged(filtersDiv).subscribe(() => {
      updateQueryBuilderLayout(filtersDiv.clientWidth);
    });


    const initializeQueryBuilder = async (libId: string, compoundType: string) => {
      ui.setUpdateIndicator(filtersDiv, true, 'Loading filters...');
      const filterFields = getDefaultProperties();
      const tagsStr = await grok.functions.call('RevvitySignalsLink:getTags', {
        assetTypeId: libId,
        type: compoundType
      });
      const tags: {[key: string]: string} = JSON.parse(tagsStr);
      Object.keys(tags).forEach((tagName) => {
        const propOptions: {[key: string]: any} = {
          name: tagName,
          type: REVVITY_FIELD_TO_PROP_TYPE_MAPPING[tags[tagName]],
        };
        const nameArr = tagName.split('.');
        if (nameArr.length > 1)
          propOptions.friendlyName = nameArr[1];
        const prop = DG.Property.fromOptions(propOptions);
        prop.options[SUGGESTIONS_FUNCTION] = async (text: string) => {
          const termsStr =  await getTerms(tagName, compoundType, libId, true);
          const terms: string[] =  JSON.parse(termsStr);
          return terms.filter((it) => it.toLowerCase().includes(text.toLowerCase()));
        }
        filterFields.push(prop);
      });
      queryBuilder = new QueryBuilder(filterFields, undefined, QueryBuilderLayout.Narrow);
      const runSearchButton = ui.bigButton('Search', async () => {
        ui.setUpdateIndicator(tv.grid.root, true, 'Searching...');
        const resultDf = await runSearch(libId, compoundType, queryBuilder!.condition);
        tv.dataFrame = resultDf;
        ui.setUpdateIndicator(tv.grid.root, false);
      });
      ui.setUpdateIndicator(filtersDiv, false);
      filtersDiv.append(queryBuilder.root);
      filtersDiv.append(ui.div(runSearchButton, {style: {paddingLeft: '4px'}}));
    }

    const initializeFilters = async () => {
      const libs = await getRevvityLibraries();
      const selectedLib = libs.filter((l) => l.name === libName);
      if (selectedLib.length) {
        //create filters button
        const filtersButton = ui.button('Add filters', () => {
          tv.dockManager.dock(filtersDiv, 'left', null, 'Filters', 0.2);
          if (!queryBuilder)
            initializeQueryBuilder(selectedLib[0].id, compoundType);
        });
        tv.setRibbonPanels([[filtersButton]]);
      }
    }

    const runSearch = async (libId: string, compoundType: string,
      queryBuilderCondition: ComplexCondition): Promise<DG.DataFrame> => {
      const condition: ComplexCondition = {
        logicalOperator: Operators.Logical.and,
        conditions: [
          materialsCondition
        ]
      }
      condition.conditions.push(
        {
          field: "assetTypeEid",
          operator: Operators.EQ,
          value: libId
        }
      );
      condition.conditions.push(
        {
          field: "type",
          operator: Operators.EQ,
          value: compoundType
        }
      );
      condition.conditions.push(queryBuilderCondition);
      const signalsQuery: SignalsSearchQuery = convertComplexConditionToSignalsSearchQuery(condition);
      console.log(signalsQuery);
      const resultDf = await grok.functions.call('RevvitySignalsLink:searchEntitiesWithStructures', {
        query: JSON.stringify(signalsQuery),
        params: '{}'
      });
      return resultDf;
    }
    initializeFilters();
  }

  for (const lib of libs) {
    const libNode = treeNode.group(lib.name);
    for (const libType of lib.types) {
      const typeNode = libNode.item(`${libType.name.charAt(0).toUpperCase()}${libType.name.slice(1)}`);
      typeNode.onSelected.subscribe(async () => {
        await createViewFromPreDefinedQuery(JSON.stringify(assetsQuery), `${libType.name.charAt(0).toUpperCase()}${libType.name.slice(1)}`,
          lib.name, libType.name);
      });
    }
  }
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
//input: string id { semType: RevvitySignalsId }
//output: widget result
export async function entityTreeWidget(id: string): Promise<DG.Widget> {
  const obj = (await queryMaterialById(id)) as RevvityApiResponse;
  const div = createRevvityResponseWidget(obj);
  return new DG.Widget(div);
}
