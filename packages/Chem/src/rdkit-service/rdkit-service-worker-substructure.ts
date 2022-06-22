import {RdKitServiceWorkerSimilarity} from './rdkit-service-worker-similarity';
import {rdKitFingerprintToBitArray} from '../utils/chem-common';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {isMolBlock} from '../utils/chem-utils';
import {RDModule, RDMol} from "../rdkit-api";

export class RdKitServiceWorkerSubstructure extends RdKitServiceWorkerSimilarity {
  _patternFps: BitArray[] | null = null;
   readonly _patternFpLength = 2048;

   constructor(module: RDModule, webRoot: string) {
    super(module, webRoot);
  }

  initMoleculesStructures(dict: string[], normalizeCoordinates: boolean, usePatternFingerprints: boolean) {
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
          const fpRdKit = mol.get_pattern_fp_as_uint8array(this._patternFpLength);
          fp = rdKitFingerprintToBitArray(fpRdKit);
        }
        if (normalizeCoordinates) {
          if (isMolBlock(item)) {
            mol.normalize_depiction();
            item = mol.compute_hash();
            mol.straighten_depiction();
            if (!hashToMolblock[item])
              hashToMolblock[item] = mol.get_molblock();
          }
        }
      } catch (e) {
        console.error('Chem | Possibly a malformed molString: `' + item + '`');
        // preserving indices with a placeholder
        mol?.delete();
        mol = this._rdKitModule.get_mol('');
        // Won't rethrow
      }
      this._rdKitMols.push(mol);
      if (this._patternFps)
        this._patternFps.push(fp!);
      if (normalizeCoordinates)
        molIdxToHash.push(item);
    }
    return {molIdxToHash, hashToMolblock};
  }

  getMol(molString: string) {
     return this._rdKitModule.get_mol(molString, '{"mergeQueryHs":true}');
  }

  getQMol(molString: string) {
     return this._rdKitModule.get_qmol(molString);
  }

  searchSubstructure(queryMolString: string, queryMolBlockFailover: string): string {
    const matches: number[] = [];
    if (this._rdKitMols) {
      let queryMol: RDMol | null = null;
      try {
        if (isMolBlock(queryMolString)) {
          queryMol = this.getMol(queryMolString);
        } else {
          let createdQMol = false;
          try {
            queryMol = this.getQMol(queryMolString);
            createdQMol = queryMol.is_valid();
          } catch (e) {
            queryMol?.delete();
            queryMol = null;
          }
          if (createdQMol) {
            let mol = null;
            let createdMol = false;
            try {
              mol = this.getMol(queryMolString);
              createdMol = mol.is_valid();
            } catch (e) {
              mol?.delete();
            }
            if (createdMol) { // check the qmol is proper
              let match = mol!.get_substruct_match(queryMol!);
              if (match === '{}') {
                queryMol!.delete();
                queryMol = mol;
              }
            } // else, this looks to be a real SMARTS
          } else { // failover to queryMolBlockFailover
            queryMol = this.getMol(queryMolBlockFailover);
          }
        }
        if (queryMol && queryMol.is_valid()) {
          if (this._patternFps) {
            const fpRdKit = queryMol.get_pattern_fp_as_uint8array(this._patternFpLength);
            const queryMolFp = rdKitFingerprintToBitArray(fpRdKit);
            for (let i = 0; i < this._patternFps.length; ++i) {
              const crossedFp = BitArray.fromAnd(this._patternFps[i], queryMolFp);
              if (crossedFp.equals(queryMolFp))
                if (this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}') // Is patternFP iff?
                  matches.push(i);
            }
          } else
            for (let i = 0; i < this._rdKitMols!.length; ++i)
              if (this._rdKitMols[i]!.get_substruct_match(queryMol) !== '{}')
                matches.push(i);
        } else
          throw new Error('Chem | Search pattern cannot be set');
      } catch (e) {
        console.error(
          'Chem | Possibly a malformed query: `' + queryMolString + '`');
        // Won't rethrow
      } finally {
        queryMol?.delete();
      }
    }
    return '[' + matches.join(', ') + ']';
  }
  
  freeMoleculesStructures(): void {
    if (this._rdKitMols !== null) {
      for (let mol of this._rdKitMols!)
        mol.delete();
      this._rdKitMols = null;
    }
    this._patternFps = null;
  }
}
