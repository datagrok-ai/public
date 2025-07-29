/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from "@datagrok-libraries/utils/src/u2";
import { buildOperatorUI, createDefaultOperator, signalsSearchBuilderUI } from './signalsSearchBuilder';
import '../css/revvity-signals-styles.css';
import { SignalsSearchParams, SignalsSearchQuery } from './signalsSearchQuery';
import { queryEntities, queryEntityById, queryMaterialById, queryStructureById, queryUsers, RevvityApiResponse, RevvityData, RevvityUser } from './revvityApi';
import { dataFrameFromObjects, reorderColummns, transformData, widgetFromObject, createRevvityResponseWidget } from './utils';
import { addMoleculeStructures, assetsQuery, batchesQuery, MOL_COL_NAME } from './compounds';


export const _package = new DG.Package();
let openedView: DG.View | null = null;

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

  const createViewFromPreDefinedQuery = async (query: string, name: string) => {
    const df = await grok.functions.call('RevvitySignalsLink:searchEntitiesWithStructures', {
      query: query,
      params: '{}'
    });
    const tv = grok.shell.addTablePreview(df);
    tv.name = name;
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


//name: Revvity Signals
//tags: panel, widgets
//input: string id { semType: RevvitySignalsId }
//output: widget result
export async function entityTreeWidget(id: string): Promise<DG.Widget> {
  const obj = (await queryMaterialById(id)) as RevvityApiResponse;
  const div = createRevvityResponseWidget(obj);
  return new DG.Widget(div);
}