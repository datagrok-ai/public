/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {axolabsStyleMap} from '../data-loading-utils/json-loader';

export class ExternalDataManager {
  private static instance: ExternalDataManager;
  private constructor() { }

  static getInstance(): ExternalDataManager {
    if (!ExternalDataManager.instance) {
      ExternalDataManager.instance = new ExternalDataManager();
    }
    return ExternalDataManager.instance;
  }

  fetchNucleotideBases(): string[] {
    const nucleotideBases: string[] = Object.keys(axolabsStyleMap);
    return nucleotideBases;
  }
}
