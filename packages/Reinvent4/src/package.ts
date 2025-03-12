/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ReinventBaseEditor, TARGET_PATH } from './utils/reinvent-editor';
import { zipFolder } from './utils/utils';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: getFolders
//output: list<string> targetFiles
export async function getFolders(): Promise<string[]> {
  const targetsFiles: DG.FileInfo[] = await grok.dapi.files.list(TARGET_PATH, true);
  return targetsFiles.filter((dir) =>  dir.isDirectory).map((dir) => dir.name);
}

//name: ReinventEditor
//tags: editor
//input: funccall call
export function reinventEditor(call: DG.FuncCall): void {
  const funcEditor = new ReinventBaseEditor();
  ui.dialog({title: 'Reinvent'})
    .add(funcEditor.getEditor())
    .onOK(async () => {
      const params = funcEditor.getParams();
      call.func.prepare({
        ligand: params.ligand,
        optimize: params.optimize
      }).call(true);
    }).show();
}

//name: runReinvent
//meta.cache: all
//meta.cache.invalidateOn: 0 * * * *
//input: string ligand {semType: Molecule}
//input: string optimize
//output: dataframe result
export async function runReinvent(ligand: string, optimize: string): Promise<DG.DataFrame> {
  const container = await grok.dapi.docker.dockerContainers.filter('reinvent').first();
  const files = (await grok.dapi.files.list(`${TARGET_PATH}/${optimize}`));

  const zipBlob = await zipFolder(files);
  const formData = new FormData();
  formData.append('folder', zipBlob, 'folder.zip');
  [ligand].forEach(smiles => {
    formData.append('smiles', smiles);
  });

  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/run_reinvent', {
    method: 'POST',
    body: formData,
  });

  const resultDf = DG.DataFrame.fromJson(await response.text());
  return resultDf;
}

//top-menu: Chem | Generate molecules...
//name: Reinvent
//tags: HitDesignerFunction
//input: string ligand = "OC(CN1CCCC1)NC(CCC1)CC1Cl" {semType: Molecule}
//input: string optimize {choices: Reinvent4:getFolders}
//editor: Reinvent4: ReinventEditor
//output: dataframe result
export async function reinvent(
  ligand: string, optimize: string
): Promise<DG.DataFrame> {
  const resultDfPromise = grok.functions.call('Reinvent4:runReinvent', {
    ligand: ligand,
    optimize: optimize
  });
  const schemas = await grok.dapi.stickyMeta.getSchemas();
  const lineageSchema = schemas.find((s) => s.name === 'Lineage');

  const resultDf: DG.DataFrame = await resultDfPromise;

  if (lineageSchema) {
    const molCol = DG.Column.fromStrings('canonical_smiles', [ligand]);
    molCol.semType = DG.SEMTYPE.MOLECULE;

    const seedDf = generateStickyDf('seed ligand');
    const generatedDf = generateStickyDf('generated');

    const lineagePromise = grok.dapi.stickyMeta.setAllValues(lineageSchema, molCol, seedDf);

    const resultMolCol = resultDf.columns.byName('SMILES');
    const initialMolCol = resultDf.columns.byName('Input_SMILES');
    resultMolCol.semType = DG.SEMTYPE.MOLECULE;
    initialMolCol.semType = DG.SEMTYPE.MOLECULE;
    await grok.data.detectSemanticTypes(resultDf);

    if (resultMolCol) {
      const lineagePromises = resultMolCol.toList().map((smiles: string) => {
        const col = DG.Column.fromStrings('smiles', [smiles]);
        col.semType = DG.SEMTYPE.MOLECULE;
        return grok.dapi.stickyMeta.setAllValues(lineageSchema, col, generatedDf);
      });

      await Promise.all([lineagePromise, ...lineagePromises]);
    }
  }

  const scoreColumn = resultDf.col('Score');
  if (scoreColumn)
    scoreColumn.meta.colors.setLinear([DG.Color.red, DG.Color.green]);

  return resultDf;
}

function generateStickyDf(role: string): DG.DataFrame {
  const lineageDf = DG.DataFrame.create(1);
  const roleCol = lineageDf.columns.addNewString('Role');
  roleCol.set(0, role);
  const optimizedCol = lineageDf.columns.addNewString('Optimized parameters');
  optimizedCol.set(0, 'Binding affinity, hERG, BBB, CYP Inhibitors');
  const authorCol = lineageDf.columns.addNewString('User');
  authorCol.set(0, grok.shell.user.friendlyName);
  const dateCol = lineageDf.columns.addNewString('Date');
  dateCol.set(0, new Date().toLocaleString());
  return lineageDf;
}
