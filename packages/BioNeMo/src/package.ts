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
  console.log(results);
}

//name: EsmFoldModel
//top-menu: Chem | BioNeMo | EsmFold...
//input: string sequence = "MDILCEENTSLSSTTNSLMQLNDDTRLYSNDFNSGEANTSDAFNWTVDSENRTNLSCEGCLSPSCLSLLHLQEKNWSALLTAVVIILTIAGNILVIMAVSLEKKLQNATNYFLMSLAIADMLLGFLVMPVSMLTILYGYRWPLPSKLCAVWIYLDVLFSTASIMHLCAISLDRYVAIQNPIHHSRFNSRTKAFLKIIAVWTISVGISMPIPVFGLQDDSKVFKEGSCLLADDNFVLIGSFVSFFIPLTIMVITYFLTIKSLQKEATLCVSDLGTRAKLASFSFLPQSSLSSEKLFQRSIHREPGSYTGRRTMQSISNEQKACKVLGIVFFLFVVMWCPFFITNIMAVICKESCNEDVIGALLNVFVWIGYLSSAVNPLVYTLFNKTYRSAFSRYIQCQYKENKKPLQLILVNTIPALAYKSSQLQMGQKKNSKQDAKTTDNDCSMVALGKQHSEEASKDNSDGVNEKVSCV"
export async function esmFoldModel(sequence: string) {
  const currentView = grok.shell.tv;
  const structure = await grok.functions.call('BioNeMo:esmfold', {sequence});
  const viewer = await currentView.dataFrame.plot.fromType('Biostructure', {pdb: structure});
  currentView.dockManager.dock(viewer.root, DG.DOCK_TYPE.RIGHT);
}