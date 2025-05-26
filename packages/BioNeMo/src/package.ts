/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { CONSTANTS, DiffDockModel, PosesJson } from './diffdock/diffdock-model';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

export async function getApiKey(): Promise<string> {
  //@ts-ignore
  const apiKey = (await _package.getSettings())['apiKey'];
  if (apiKey)
    return apiKey;
  throw new Error('Api key not set');
}

//name: MolMIMModel
//input: string algorithm = "CMA-ES"
//input: int num_molecules = 30
//input: string property_name = "QED"
//input: bool minimize = false
//input: double min_similarity = 0.3
//input: int particles = 30
//input: int iterations = 10
//input: string smi = "[H][C@@]12Cc3c[nH]c4cccc(C1=C[C@H](NC(=O)N(CC)CC)CN2C)c34" {semType: Molecule}
export async function molMIMModel(algorithm: string, num_molecules: number, property_name: string, minimize: boolean, min_similarity: number,
  particles: number, iterations: number, smi: string
) {
  const apiKey = await getApiKey();
  const results = await grok.functions.call('BioNeMo:MolMIMGenerate', { algorithm, num_molecules, property_name, minimize, min_similarity, particles, iterations, smi, apiKey });
}

//name: EsmFoldModel
//top-menu: Bio | Folding | EsmFold...
//input: dataframe df 
//input: column sequences {semType: Macromolecule}
export async function esmFoldModel(df: DG.DataFrame, sequences: DG.Column) {
  try {
    const apiKey = await getApiKey();
    const grid = grok.shell.getTableView(df.name).grid;
    const protein = DG.Column.fromType(DG.TYPE.STRING, 'Protein', sequences.length);
    for (let i = 0; i < sequences.length; ++i) {
      const colValue = sequences.get(i);
      const predictedValue = await grok.functions.call('BioNeMo:esmfold', { sequence: colValue, api_key: apiKey });
      protein.set(i, predictedValue);
    }
    protein.setTag(DG.TAGS.SEMTYPE, DG.SEMTYPE.MOLECULE3D);
    df.columns.add(protein);
    await grok.data.detectSemanticTypes(df);
    grid.invalidate();
  } catch (e: any) {
    console.error('Error running ESMFold model:', e.message);
    grok.shell.error(`Failed to run ESMFold model: ${e.message}`);
  }
}

//name: Bio | EsmFold
//input: semantic_value sequence {semType: Macromolecule}
//output: widget result
export async function esmFoldModelPanel(sequence: DG.SemanticValue): Promise<DG.Widget> {
  const result = new DG.Widget(ui.div());
  const loader = ui.loader();
  result.root.appendChild(loader);

  try {
    const apiKey = await getApiKey();
    const res = await grok.functions.call('BioNeMo:esmfold', { sequence: sequence.value, api_key: apiKey });
    const { success, error, pdb } = JSON.parse(res);

    result.root.removeChild(loader);

    if (error)
      result.root.appendChild(ui.divText(error));
    else {
      const molstarViewer = await sequence.cell.dataFrame.plot.fromType('Biostructure', { pdb });
      result.root.appendChild(molstarViewer.root);
    }
  } catch (e: any) {
    result.root.removeChild(loader);
    result.root.appendChild(ui.divText(e.message));
  }

  return result;
}

//name: getTargetFilses
//output: list<string> targetFiles
export async function getTargetFiles(): Promise<string[]> {
  const targetsFiles: DG.FileInfo[] = await grok.dapi.files.list(CONSTANTS.TARGET_PATH, true);
  return targetsFiles.filter(file => file.isDirectory).map(file => file.name);
}

//name: diffDockModelScript
//meta.cache: client
//meta.cache.invalidateOn: 0 * * * *
//input: string ligand
//input: string target
//input: int poses
//output: string result
export async function diffDockModelScript(ligand: string, target: string, poses: number): Promise<string | undefined> {
  try {
    const apiKey = await getApiKey();
    const encodedPoses = await grok.functions.call('Bionemo:diffdock', {
      protein: target,
      ligand: ligand,
      num_poses: poses,
      api_key: apiKey
    });
    return new TextDecoder().decode(encodedPoses.data);
  } catch (e: any) {
    console.error('Error running DiffDock model:', e.message);
    throw new Error(`Failed to run DiffDock model: ${e.message}`);
  }
}

//name: DiffDockModel
//top-menu: Chem | Docking | DiffDock...
//input: dataframe df
//input: column ligands {semType: Molecule}
//input: string target {choices: Bionemo: getTargetFiles}
//input: int poses = 5
export async function diffDockModel(df: DG.DataFrame, ligands: DG.Column, target: string, poses: number) {
  const receptorFile = (await grok.dapi.files.list(`${CONSTANTS.TARGET_PATH}/${target}`)).find((file) => file.extension === 'pdbqt')!;
  const receptor = await grok.dapi.files.readAsText(receptorFile);
  const diffDockModel = new DiffDockModel(df, ligands, receptor, receptorFile.name, poses);
  await diffDockModel.run();
}

//name: Biology | DiffDock
//tags: panel, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export async function diffDockPanel(smiles: DG.SemanticValue): Promise<DG.Widget> {
  const posesInput = ui.input.int('Poses', { value: 10 });
  const targetInput = ui.input.choice('Target', { value: (await getTargetFiles())[0], items: await getTargetFiles() });

  const resultsContainer = ui.div();
  const form = ui.form([targetInput, posesInput]);
  const panels = ui.divV([form, ui.button('Run', async () => {
    await handleRunClick(smiles, posesInput.value!, targetInput.value!, resultsContainer);
  }), resultsContainer]);

  return DG.Widget.fromRoot(panels);
}

async function handleRunClick(smiles: DG.SemanticValue, poses: number, target: string, resultsContainer: HTMLDivElement) {
  resultsContainer.innerHTML = '';
  const loader = ui.loader();
  resultsContainer.appendChild(loader);

  try {
    const { value: smilesValue, units: smilesUnits, cell } = smiles;
    const { dataFrame: table, column, rowIndex } = cell;

    const receptorFile = (await grok.dapi.files.list(`${CONSTANTS.TARGET_PATH}/${target}`)).find((file: DG.FileInfo) => file.extension === 'pdbqt')!;
    const receptor = await grok.dapi.files.readAsText(receptorFile);
    const { name: receptorName } = receptorFile;

    const diffDockModel = new DiffDockModel(table, column, receptor, receptorName, poses);
    const virtualColName = getVirtualPosesColumnName(receptorName, poses);

    let posesColumn = table.columns.byName(virtualColName);
    let posesJson: PosesJson;

    if (!posesColumn) {
      posesJson = await diffDockModel.getPosesJson(smilesValue, smilesUnits);
      posesColumn = await diffDockModel.createColumn(DG.TYPE.STRING, virtualColName, table.rowCount);
      table.columns.add(posesColumn);
      posesColumn.set(rowIndex, JSON.stringify(posesJson));
    } else {
      const existingValue = posesColumn.get(rowIndex);
      if (!existingValue) {
        posesJson = await diffDockModel.getPosesJson(smilesValue, smilesUnits);
        posesColumn.set(rowIndex, JSON.stringify(posesJson));
      } else {
        posesJson = JSON.parse(existingValue);
      }
    }

    diffDockModel.virtualPosesColumn = posesColumn;

    const { bestId, bestPose } = diffDockModel.findBestPose(posesJson);
    const viewer = await diffDockModel.createCombinedControl(posesJson, bestPose, bestId, false);

    resultsContainer.removeChild(loader);
    resultsContainer.append(viewer);
  } catch (e: any) {
    resultsContainer.removeChild(loader);
    const errorDiv = ui.divText(e.message)
    resultsContainer.append(errorDiv);
  }
}

function getVirtualPosesColumnName(target: string, poses: number): string {
  return `${CONSTANTS.VIRTUAL_POSES_COLUMN_NAME}_${target}_${poses}`;
}