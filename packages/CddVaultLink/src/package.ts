/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from "@datagrok-libraries/utils/src/u2";
import {getVaults, queryMolecules} from "./cdd-vault-api";

export const _package = new DG.Package();

//tags: app
//name: CDD Vault
//meta.icon: images/cdd-icon-small.png
//output: view v
//meta.browsePath: Chem
export async function cddVaultApp(): Promise<DG.ViewBase> {

  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/images/cdd-icon-big.png',
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/CddVaultLink/README.md',
    description: '- Integrate with your CDD Vault.\n' +
      '- Analyze assay data.\n' +
      '- Find contextual information on molecules.\n' +
      '- Browse the vault content.\n'
  });

  return DG.View.fromRoot(appHeader);
}

//input: dynamic treeNode
//input: view browseView
export async function cddVaultAppTreeBrowser(treeNode: DG.TreeViewGroup, browseView: any) {
  const vaults = await getVaults();
  console.log(vaults);

  for (const vault of vaults.data!) {
    const vaultNode = treeNode.group(vault.name);
    const moleculesNode = vaultNode.item('Molecules');
    moleculesNode.onSelected.subscribe(async (_) => {
      const molecules = await queryMolecules(vault.id, { page_size: 100 });
      console.log(molecules);
      const df = DG.DataFrame.fromObjects(molecules.data!.objects!)!;
      grok.shell.addTableView(df);
    });

    vaultNode.group('Protocols');
    vaultNode.group('Plates');
    vaultNode.group('Assays');
  }
}


//name: Databases | CDD Vault
//input: string mol {semType: Molecule}
//tags: panel
//output: widget result
export function molColumnPropertyPanel(molecule: string): DG.Widget {
  return DG.Widget.fromRoot(ui.wait(async () => {
    const vaults = await getVaults();
    const vaultId = vaults.data![0].id;
    const cddMols = await queryMolecules(vaultId, { structure: molecule, structure_search_type: "exact"});

    if (!cddMols.data?.objects?.length)
      return ui.divText('Not found');

    return ui.tableFromMap(cddMols.data.objects[0]);
  }))
}