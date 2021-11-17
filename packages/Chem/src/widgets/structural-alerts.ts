// The file is imported from a WebWorker. Don't use Datagrok imports
import { drawMoleculeToCanvas } from '../chem_common_rdkit';

let _structuralAlertsRdKitModule: any = null;
let _webRoot: string | null = null;
let _smartsMap: Map<string, any> = new Map();
let _data: string[] | null = null;

export function setStructuralAlertsRdKitModule(module: any, webRoot: string) {
  _structuralAlertsRdKitModule = module;
  _webRoot = webRoot;
}

export function loadAlertsCollection(smarts: string[]) {
  _data = smarts;
  for (let i = 0; i < smarts.length; i++) {
    const currentSmarts = smarts[i];
    _smartsMap.set(currentSmarts, _structuralAlertsRdKitModule.get_qmol(currentSmarts));
  }
}

export function getStructuralAlerts(smiles: string) {
  const alerts: number[] = [];
  const mol = _structuralAlertsRdKitModule.get_mol(smiles);
  //TODO: use SustructLibrary and count_matches instead. Currently throws an error on rule id 221
  // const lib = new _structuralAlertsRdKitModule.SubstructLibrary();
  // lib.add_smiles(smiles);
  for (let i = 0; i < _data!.length; i++) {
    const subMol = _smartsMap.get(_data![i]);
    // lib.count_matches(subMol);
    const matches = mol.get_substruct_matches(subMol);
    if (matches !== '{}') {
      alerts.push(i);
    }
  }
  return alerts;
}