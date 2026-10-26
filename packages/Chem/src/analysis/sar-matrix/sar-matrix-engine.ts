import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getRdKitModule} from '../../utils/chem-common-rdkit';
import {byEngineModule, withMol} from '../../engine/chem-engine';
import {hasNewLines} from '../../utils/chem-common';
import {MAX_SMILES_LENGTH} from '../../utils/chem-constants';

/**
 * `read(mol, module)` of a structure parsed on the main thread by the "Chemical engine" (src/engine/README.md): Crux's
 * module in Crux mode, the module the build's other structures come from (the workers' canonical SMILES, fragments and
 * links), and RDKit's where Crux reads none or throws; RDKit's in RDKit mode. Null where the structure does not read
 * (`get_mol`'s must be valid; `get_qmol`'s when `asQuery`), `read` gives null or throws on RDKit. `read` parses any
 * other structure it needs on `module`.
 */
export function withEngineMol<T>(call: string, smiles: string, read: (mol: RDMol, module: RDModule) => T | null,
  asQuery = false): T | null {
  return byEngineModule(call, getRdKitModule(), (module, crux) => {
    if (crux && !hasNewLines(smiles) && smiles.length > MAX_SMILES_LENGTH)
      return null; // RDKit's: a SMILES this long could run Crux's wasm out of memory
    return withMol(module, crux, smiles, (mol) => asQuery || mol.is_valid() ? read(mol, module) : null, asQuery);
  });
}
