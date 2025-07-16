/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from "@datagrok-libraries/utils/src/u2";
import { buildOperatorUI, createDefaultOperator, signalsSearchBuilderUI } from './signalsSearchBuilder';
import '../css/revvity-signals-styles.css';


export const _package = new DG.Package();
const STORAGE_NAME = 'Revvity';
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
    const queryBuilder = signalsSearchBuilderUI(() => {
      console.log('**********on change');
    });
    v.append(queryBuilder);
    grok.shell.addPreview(v);
  })
}