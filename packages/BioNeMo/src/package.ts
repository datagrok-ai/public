/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: MolMIMModel
//top-menu: Chem | BioNeMo | MolMIM...
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
  const results = await grok.functions.call('BioNeMo:MolMIMGenerate', {algorithm, num_molecules, property_name, minimize, min_similarity, particles, iterations, smi});
}

//name: EsmFoldModel
//top-menu: Bio | BioNeMo | EsmFold...
//input: dataframe df 
//input: column sequences {semType: Macromolecule}
export async function esmFoldModel(df: DG.DataFrame, sequences: DG.Column) {
  const grid = grok.shell.getTableView(df.name).grid;
  const protein = DG.Column.fromType(DG.TYPE.STRING, 'Protein', sequences.length);
  //init not working
  /*protein.init((i: number) => {
    const value = sequences.get(i);
    grok.functions.call('BioNeMo:esmfold', {value}).then((res) => res);
  });*/
  for (let i = 0; i < sequences.length; ++i) {
    const colValue = sequences.get(i);
    const predictedValue = await grok.functions.call('BioNeMo:esmfold', {sequence: colValue});
    protein.set(i, predictedValue);
  }
  protein.setTag(DG.TAGS.SEMTYPE, DG.SEMTYPE.MOLECULE3D);
  df.columns.add(protein);
  await grok.data.detectSemanticTypes(df);
  grid.invalidate();
}

//name: EsmFoldModelPanel
//input: semantic_value sequence {semType: Macromolecule}
//output: widget result
export async function esmFoldModelPanel(sequence: DG.SemanticValue): Promise<DG.Widget> {
  const result = new DG.Widget(ui.div());
  const loader = ui.loader();
  result.root.appendChild(loader);
  grok.functions.call('BioNeMo:esmfold', {sequence: sequence.value}).then(async (res) => {
    result.root.removeChild(loader);
    const molstarViewer = await sequence.cell.dataFrame.plot.fromType('Biostructure', {pdb: res});
    result.root.appendChild(molstarViewer.root);
  });
  return result;
}

//name: DiffDockModel
//top-menu: Chem | BioNeMo | DiffDock...
//input: dataframe df
//input: column ligands {semType: Molecule}
//input: file target
//input: int poses = 20
export async function diffDockModel(df: DG.DataFrame, ligands: DG.Column, target: DG.FileInfo, poses: number) {
  const encodedPoses = await grok.functions.call('BioNeMo:diffdock', {protein: await target.readAsString(), ligand: ligands.get(0), num_poses: poses});
  const posesJson = new TextDecoder().decode(encodedPoses.data);
}