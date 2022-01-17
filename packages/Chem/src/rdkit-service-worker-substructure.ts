import {RdKitServiceWorkerSimilarity} from './rdkit-service-worker-similarity';
import {rdKitFingerprintToBitArray} from './chem-common';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {isMolBlock} from './chem-utils';

export class RdkitServiceWorkerSubstructure extends RdKitServiceWorkerSimilarity {
  _patternFps: BitArray[] | null = null;
   readonly _patternFpLength = 2048;

   constructor(module: Object, webRoot: string) {
    super(module, webRoot);
  }

  initMoleculesStructures(dict: string[], usePatternFingerprints: boolean) {
    this.freeMoleculesStructures();
    this.freeMorganFingerprints();
    if (dict.length === 0)
      return;
    this._rdKitMols = [];
    if (usePatternFingerprints)
      this._patternFps = [];
    const hashToMolblock: {[_:string] : any} = {};
    const molIdxToHash = [];
    for (let i = 0; i < dict.length; ++i) {
      let item = dict[i];
      let mol = null;
      let fp: BitArray | null = null;
      if (usePatternFingerprints)
        fp = new BitArray(this._patternFpLength);
      try {
        mol = this._rdKitModule.get_mol(item);
        if (usePatternFingerprints) {
          const fpRdKit = mol.get_pattern_fp(this._patternFpLength);
          fp = rdKitFingerprintToBitArray(fpRdKit);
        }
        if (isMolBlock(item)) {
          item = mol.normalize_2d_molblock();
          mol.straighten_2d_layout();
          if (!hashToMolblock[item])
            hashToMolblock[item] = mol.get_molblock();
        }
      } catch (e) {
        console.error('Chem | Possibly a malformed molString: `' + item + '`');
        // preserving indices with a placeholder
        mol?.delete();
        mol = this._rdKitModule.get_mol('');
        // Won't rethrow
      }
      this._rdKitMols.push(mol);
      if (this._patternFps) {
        this._patternFps.push(fp!);
      }
      molIdxToHash.push(item);
    }
    return {molIdxToHash, hashToMolblock};
  }

  searchSubstructure(queryMolString: string, querySmarts: string) {
    const matches: number[] = [];
    if (this._rdKitMols) {
      try {
        let queryMol = null;
        try {
          try {
            queryMol = this._rdKitModule.get_mol(queryMolString, '{"mergeQueryHs":true}');
          } catch (e2) {
            queryMol?.delete();
            queryMol = null;
            if (querySmarts !== null && querySmarts !== '') {
              console.log('Chem | Cannot parse a MolBlock. Switching to SMARTS');
              queryMol = this._rdKitModule.get_qmol(querySmarts);
            } else
              throw new Error('Chem | SMARTS not set');
          }
          if (queryMol && queryMol.is_valid()) {
            if (this._patternFps) {
              const fpRdKit = queryMol.get_pattern_fp(this._patternFpLength);
              const queryMolFp = rdKitFingerprintToBitArray(fpRdKit);
              for (let i = 0; i < this._patternFps.length; ++i) {
                const crossedFp = BitArray.fromAnd(this._patternFps[i], queryMolFp);
                if (crossedFp.equals(queryMolFp))
                  if (this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}') // Is patternFP iff?
                    matches.push(i);
              }
            } else {
              for (let i = 0; i < this._rdKitMols!.length; ++i)
                if (this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}')
                  matches.push(i);
            }
          }
        } finally {
          queryMol?.delete();
        }
      } catch (e) {
        console.error(
          'Possibly a malformed query: `' + queryMolString + '`');
        // Won't rethrow
      }
    }
    return '[' + matches.join(', ') + ']';
  }

  freeMoleculesStructures() {
    if (this._rdKitMols !== null) {
      for (let mol of this._rdKitMols!)
        mol.delete();
      this._rdKitMols = null;
    }
    this._patternFps = null;
  }
}
