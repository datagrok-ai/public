import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {delay} from '@datagrok-libraries/utils/src/test';
import {CACHED_DOCKING, POSE_COL, TARGET_PATH, _package} from '../utils/constants';
import { BiostructureData } from '@datagrok-libraries/bio/src/pdb/types';
import { AutoDockDataType } from '../apps/auto-dock-app';
import { runAutodock5 } from '../package';

async function openMoleculeDataset(name: string): Promise<DG.TableView> {
  const table = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(name));
  grok.shell.windows.showProperties = true;
  return grok.shell.addTableView(table);
}

export async function _demoDocking(): Promise<void> {
  const demoScript = new DemoScript('Docking',
    '', undefined, {autoStartFirstStep: true});
  let table: DG.DataFrame;
  let tv: DG.TableView;
  demoScript
    .step('Load data', async () => {
      tv = await openMoleculeDataset('System:AppData/Docking/demo_files/demo_dataset.csv');
      table = tv.dataFrame;
    }, {description: 'Load dataset that contains compounds designed BACE1 Inhibition.\
      Each entry includes a compound ID and its corresponding SMILES representation.'})
    .step('Perform docking', async () => {
      const autodockResults = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText('System:AppData/Docking/demo_files/autodock_results.csv'));
      const targetsFiles: DG.FileInfo[] = await grok.dapi.files.list(TARGET_PATH, true);
      const receptor = targetsFiles.filter(file => file.isFile).find(file => file.name === 'BACE1.pdbqt');
      const target = targetsFiles.filter(file => file.isDirectory).find(file => file.name === 'BACE1');
      const receptorData: BiostructureData = {
        binary: false,
        data: (await grok.dapi.files.readAsText('System:AppData/Docking/targets/BACE1/BACE1.pdbqt')),
        ext: receptor!.extension,
        options: {name: receptor!.name,},
      };
      const data: AutoDockDataType = {
        ligandDf: table,
        ligandMolColName: 'SMILES',
        receptor: receptorData,
        gpfFile: (await grok.dapi.files.readAsText('System:AppData/Docking/targets/BACE1/BACE1.gpf')),
        confirmationNum: 10,
        ligandDfString: table.columns.byName('SMILES').toString(),
      };
      //@ts-ignore
      CACHED_DOCKING.K.push(data);
      //@ts-ignore
      CACHED_DOCKING.V.push(autodockResults);
      await runAutodock5(table, table.col('SMILES')!, target!.name, 10);
      tv.dataFrame.col(POSE_COL)!.semType = DG.SEMTYPE.MOLECULE3D;
      tv.dataFrame.col(POSE_COL)!.setTag(DG.TAGS.UNITS, 'pdbqt');
      tv.grid.invalidate();
      }, {description: 'Conducting docking analysis is executed through the **Chem > Autodock** module,\
        involving the configuration of a dataframe representing unique molecular setups, ligand identification,\
        and the specification of the target folder housing docking configurations and the macromolecule.\
        The analysis yields valuable results in the form of molecular poses and binding energies, providing\
        insights into the potential interactions between ligands and the macromolecular target.'})
      .step('Explore results', async () => {
        tv.dataFrame.currentCell = tv.dataFrame.cell(0, POSE_COL);
        await delay(100);
      }, {description: 'After clicking on a pose, users can view the pocket, position, and additional characteristics\
        like electrostatic properties and fixed ligand details for in-depth analysis.', delay: 100})
      .start();
  }