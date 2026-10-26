/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import MolNotation = DG.chem.Notation;
// temporary polyfills for extended notations, remove these when js-api is updated. added to avoid local dependency
// @ts-ignore
MolNotation.CxSmarts = 'cxsmarts';
// @ts-ignore
MolNotation.CxSmiles = 'cxsmiles';


// enum MolNotation {
//   Smiles = 'smiles',
//   Smarts = 'smarts',
//   MolBlock = 'molblock', // molblock V2000
//   V3KMolBlock = 'v3Kmolblock', // molblock V3000
//   Unknown = 'unknown',
//   CxSmiles = 'cxsmiles', // extended smiles
//   CxSmarts = 'cxsmarts', // extended smarts
// }

// datagrok libraries dependencies
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {getMolSafe, isSmilesSyntax} from './mol-creation_rdkit';
import {byEngineModule} from '../engine/chem-engine';
import {noteCruxError} from '../engine/crux-module';
import {convertMolNotationByEngine} from '../engine/worker-jobs';
import {cruxReadsCxExtension} from '../crux/crux-service';


export const MALFORMED_MOL_V2000 = `
Malformed

  0  0  0  0  0  0  0  0  0  0999 V2000
M  END`;
const MALFORMED_MOL_V3000 = `
Malformed

  0  0  0  0  0  0            999 V3000
M  END`;

const cxWithoutCoordsJSON = JSON.stringify({CX_ALL_BUT_COORDS: true});

/**
 * Convert between the following notations: SMILES, SMARTS, Molfile V2000 and Molfile V3000
 *
 * @param {string} moleculeString  admissible notations: listed below
 * @param {string} sourceNotation  possible values: 'smiles', 'smarts',
 * 'molblock' (corresponding to Molfile V2000) and 'V3KMolBlock' (corresponding to
 * Molfile V3000)
 * @param {string} targetNotation  possible values: same as for sourceNotation
 * @param {RDModule} rdKitModule
 * @param {boolean} addHs
 * @return {string} the converted representation
 *
 * In Crux mode ("Chemical Engine", engine/README.md) the Crux module converts (RDKit.js 2026.03.6's strings), and
 * `rdKitModule` what Crux cannot read or refuses, and a CXSMILES whose query features Crux does not read.
 */
export function _convertMolNotation(
  moleculeString: string,
  sourceNotation: DG.chem.Notation,
  targetNotation: DG.chem.Notation,
  rdKitModule: RDModule,
  addHs: boolean = false,
): string {
  if (sourceNotation === targetNotation)
    throw new Error(`Convert molecule notation: source and target notations must differ: "${sourceNotation}"`);
  if (!cruxReadsCxExtension(moleculeString)) // query features Crux's reader skips: RDKit's in both modes
    return _convertMolNotationByRdKit(moleculeString, sourceNotation, targetNotation, rdKitModule, addHs);
  const none = notConverted(targetNotation);
  return byEngineModule('convertMolNotation', rdKitModule, (module, crux) => {
    const converted = convertWith(moleculeString, sourceNotation, targetNotation, module, addHs, crux);
    return crux || converted !== none ? converted : null; // (RDKit's `none` is no result: no fallback is logged)
  }, `RDKit converted a molecule Crux cannot read: ${moleculeString}`) ?? none;
}

/** What a conversion gives a molecule it cannot read or write in `targetNotation`. */
function notConverted(targetNotation: DG.chem.Notation): string {
  return (targetNotation === MolNotation.MolBlock) ? MALFORMED_MOL_V2000 :
    (targetNotation === MolNotation.V3KMolBlock) ? MALFORMED_MOL_V3000 : 'MALFORMED_INPUT_VALUE';
}

/** {@link _convertMolNotation} by RDKit in both modes: for code whose next step is RDKit's (reaction drawing, the
 * substructure filter's query names). */
export function _convertMolNotationByRdKit(moleculeString: string, sourceNotation: DG.chem.Notation,
  targetNotation: DG.chem.Notation, rdKitModule: RDModule, addHs: boolean = false): string {
  if (sourceNotation === targetNotation)
    throw new Error(`Convert molecule notation: source and target notations must differ: "${sourceNotation}"`);
  return convertWith(moleculeString, sourceNotation, targetNotation, rdKitModule, addHs, false)!;
}

/** {@link _convertMolNotation} by one module. Under Crux (`crux`): null for a molecule it cannot read, and what it
 * throws thrown, so that RDKit converts the molecule. */
function convertWith(moleculeString: string, sourceNotation: DG.chem.Notation, targetNotation: DG.chem.Notation,
  module: RDModule, addHs: boolean, crux: boolean): string | null {
  let result = notConverted(targetNotation);
  let mol: RDMol | null = null;
  try {
    mol = getMolSafe(moleculeString, {}, module, true, crux ? noteCruxError : undefined).mol;
    if (!mol && crux)
      return null;
    if (mol) {
      if (targetNotation === MolNotation.MolBlock) {
        //when converting from smiles set coordinates and rendering parameters
        if (sourceNotation === MolNotation.Smiles) {
          if (!mol.has_coords())
            mol!.set_new_coords();
          mol.normalize_depiction(1);
          mol.straighten_depiction(false);
          // SMILES stereo is absolute, and so is the molblock's (its chiral flag): without the flag, a reader that honours
          // it (Ketcher's Indigo) took the centres as relative, and the molecule cell editor wrote them back as an AND group
          if (moleculeString.includes('@'))
            mol.set_prop('_MolFileChiralFlag', '1');
        }
        if (addHs) {
          try {
            mol.add_hs_in_place();
          } catch (e) {}
        }
        result = mol.get_molblock();
      }
      if (targetNotation === MolNotation.Smiles)
        result = mol.get_smiles();
      if (targetNotation === MolNotation.V3KMolBlock)
        result = mol.get_v3Kmolblock();
      if (targetNotation === MolNotation.Smarts)
        result = mol.get_smarts(); // @ts-ignore
      if (targetNotation === MolNotation.CxSmiles) // @ts-ignore temporary until chem-meta update
        result = mol.get_cxsmiles(cxWithoutCoordsJSON); // @ts-ignore
      if (targetNotation === MolNotation.CxSmarts)
        result = mol.get_cxsmarts();
      // if (targetNotation === MolNotation.Inchi)
      //   return mol.get_inchi();
    }
  } catch (err) {
    if (crux)
      throw err;
    console.error(errorToConsole(err));
  } finally {
    mol?.delete();
  }
  return result;
}

export function molToMolblock(molStr: string, module: RDModule): string {
  return DG.chem.isMolBlock(molStr) ?
    molStr : _convertMolNotation(molStr, MolNotation.Unknown, MolNotation.MolBlock, module);
}

export function molToSmiles(molStr: string, module: RDModule): string {
  return _convertMolNotation(molStr, MolNotation.Unknown, MolNotation.Smiles, module);
}

/** A molecule string as SMILES, as Copy SMILES gives it: the string itself when it is written in SMILES syntax alone, else
 * RDKit's SMILES of it (a molblock, a SMARTS, a CXSMILES, whose extension a SMILES does not carry: one was copied as it
 * stood unless SMARTS syntax was read in its extension, `&1` or `$Q_e$`). */
export function asSmiles(molStr: string, module: RDModule): string {
  return isSmilesSyntax(molStr) ? molStr : molToSmiles(molStr, module);
}

export async function convertNotationForColumn(molecules: DG.Column<string>, targetNotation: MolNotation, kekulize = false): Promise<string[]> {
  return convertMolNotationByEngine(molecules.toList(), targetNotation, kekulize);
}
