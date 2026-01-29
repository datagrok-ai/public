/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ReinventBaseEditor, TARGET_PATH } from './utils/reinvent-editor';
import { zipFolder } from './utils/utils';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions
{
  @grok.decorators.func()
  static info() {  
    grok.shell.info(_package.webRoot);
  }

  @grok.decorators.func()
  static async getFolders(): Promise<string[]> {  
    const targetsFiles: DG.FileInfo[] = await grok.dapi.files.list(TARGET_PATH, true);
    return targetsFiles.filter((dir) =>  dir.isDirectory).map((dir) => dir.name);
  }

  @grok.decorators.func({
    'name': 'ReinventEditor',
    meta: {role: 'editor'},
  })
  static reinventEditor(
    call: DG.FuncCall): void {
    const funcEditor = new ReinventBaseEditor();
    ui.dialog({title: 'Generate molecules'})
      .add(funcEditor.getEditor())
      .onOK(async () => {
        const params = funcEditor.getParams();
        call.func.prepare({
          ligand: params.ligand,
          optimize: params.optimize
        }).call(true);
      }).show();
  }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
  static async runReinvent(
    @grok.decorators.param({'options':{'semType':'Molecule'}})  ligand: string,
    optimize: string): Promise<string> {
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

    const resultText = await response.text();
    return resultText;
  }

  @grok.decorators.func({
    'name': 'Reinvent',
    'editor': 'Reinvent4:ReinventEditor',
    'meta': {'role': 'hitDesignerFunction'}
  })
  static async reinvent(
    @grok.decorators.param({'options':{'semType':'Molecule','initialValue':'\'OC(CN1CCCC1)NC(CCC1)CC1Cl\''}}) ligand: string,
    @grok.decorators.param({'options':{'choices':'Reinvent4:getFolders'}}) optimize: string
  ): Promise<DG.DataFrame> {  
    const resultDfPromise = grok.functions.call('Reinvent4:runReinvent', {
      ligand: ligand,
      optimize: optimize
    });
    const schemas = await grok.dapi.stickyMeta.getSchemas();
    const lineageSchema = schemas.find((s) => s.name === 'Lineage');

    const resultDf: DG.DataFrame = DG.DataFrame.fromJson(await resultDfPromise);

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

  @grok.decorators.func({
    'top-menu': 'Chem | Generate molecules...',
    'editor': 'Reinvent4:ReinventEditor'
  })
  static async reinventTopMenu(
    @grok.decorators.param({'options':{'semType':'Molecule','initialValue':'\'OC(CN1CCCC1)NC(CCC1)CC1Cl\''}})  ligand: string,
    @grok.decorators.param({'options':{'choices':'Reinvent4:getFolders'}})   optimize: string) {  
    const generatedDf = await PackageFunctions.reinvent(ligand, optimize);
    grok.shell.addTableView(generatedDf);
  }
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
