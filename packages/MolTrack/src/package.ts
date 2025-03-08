/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from "@datagrok-libraries/utils/src/u2";

export const _package = new DG.Package();

//tags: app
//name: CDD Vault
//meta.icon: images/cdd-icon-small.png
//output: view v
//meta.browsePath: Chem
export async function molTrackApp(): Promise<DG.ViewBase> {

  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/images/cdd-icon-big.png',
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/MolTrack/README.md',
    description: '- Chemical compound registration system\n' +
      '- Analyze assay data\n' +
      '- Find contextual information on molecules.\n'
  });


  return DG.View.fromRoot(ui.divV([
    appHeader,
    ui.wait(async () => (await grok.functions.call('MolTrack:getCompounds') as DG.DataFrame).plot.grid().root)
  ]));
}


//input: dynamic treeNode
//input: view browseView
export async function cddVaultAppTreeBrowser(appNode: DG.TreeViewGroup, browseView: any) {
  appNode.group('Protocols');
  appNode.group('Plates');
  appNode.group('Assays');
}