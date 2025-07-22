/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from "@datagrok-libraries/utils/src/u2";
import { buildOperatorUI, createDefaultOperator, signalsSearchBuilderUI } from './signalsSearchBuilder';
import '../css/revvity-signals-styles.css';
import { SignalsSearchParams, SignalsSearchQuery } from './signalsSearchQuery';
import { queryEntities } from './revvityApi';
import { dataFrameFromObjects } from './utils';
import { assetsQuery, batchesQuery } from './compounds';


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
    //const v = DG.View.create(name);
    const df = await grok.functions.call('RevvitySignalsLink:searchEntities', {
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
  } catch (e: any) {
    grok.shell.error(e?.message ?? e);
  }
  return df;
}