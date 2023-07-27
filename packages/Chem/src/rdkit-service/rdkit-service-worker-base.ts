import {RDModule, RDMol, RDSubstructLibrary} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import { LRUCache } from 'typescript-lru-cache';

export const MAX_MOL_CACHE_SIZE = 10000;

export class RdKitServiceWorkerBase {
  _rdKitModule: RDModule;
  _webRoot: string;
  _rdKitMols: (RDMol | null)[] | null = null;
  _molsCache: LRUCache<string, RDMol> | null = null;
  constructor(module: RDModule, webRoot: string) {
    this._rdKitModule = module;
    this._webRoot = webRoot;
    this._molsCache = new LRUCache<string, RDMol>({
      maxSize: MAX_MOL_CACHE_SIZE,
      onEntryEvicted: ({ key, value, isExpired }) => {
        value?.delete();
      },
    });  
  }
}
