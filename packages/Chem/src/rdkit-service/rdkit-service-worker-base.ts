import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {LRUCache} from 'typescript-lru-cache';

export const MAX_MOL_CACHE_SIZE = 10000;

export class RdKitServiceWorkerBase {
  _rdKitModule: RDModule;
  _webRoot: string;
  _rdKitMols: (RDMol | null)[] | null = null;
  _molsCache: LRUCache<string, RDMol> | null = null;
  _cacheCounter = 0;
  _requestTerminated = false;
  /** Per-op cancellation set. Lets one op be cancelled without affecting concurrent work. */
  _terminatedOps: Set<string> = new Set<string>();
  _terminationCheckDelay = 50;

  isOpTerminated(opId?: string): boolean {
    return opId !== undefined && this._terminatedOps.has(opId);
  }

  setOpTerminate(opId: string, flag: boolean): void {
    if (flag) this._terminatedOps.add(opId);
    else this._terminatedOps.delete(opId);
  }

  constructor(module: RDModule, webRoot: string) {
    this._rdKitModule = module;
    this._webRoot = webRoot;
    this._molsCache = new LRUCache<string, RDMol>({
      maxSize: MAX_MOL_CACHE_SIZE,
      onEntryEvicted: ({key, value}) => {
        value?.delete();
        console.log(`${key} deleted`);
      },
    });
  }

  addToCache(rdMol: RDMol): boolean {
    let added = false;
    const key = rdMol.get_smiles();
    if (this._molsCache?.has(key))
      return true;
    if (this._cacheCounter < MAX_MOL_CACHE_SIZE && !this._molsCache?.has(key)) {
      this._molsCache?.set(key, rdMol!);
      this._cacheCounter++; //need this additional counter (instead of i) not to consider empty molecules
      added = true;
    }
    return added;
  }
}
