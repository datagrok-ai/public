/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import {PATTERN_KEY } from './const';
import {PatternConfiguration} from './types';
import {DEFAULT_SEQUENCE_LENGTH, DEFAULT_PHOSPHOROTHIOATE} from './const';
import {ExternalDataManager} from './external-data-manager';

export class PatternConfigurationManager {
  constructor() {
    const defaultPatternConfiguration = this.getDefaultConfiguration();
    this._patternConfiguration = new rxjs.BehaviorSubject<PatternConfiguration>(defaultPatternConfiguration);
  }

  private _patternConfiguration: rxjs.BehaviorSubject<PatternConfiguration>;
  private externalDataManager = new ExternalDataManager();

  getCurrentConfiguration(): PatternConfiguration {
    return this._patternConfiguration.getValue();
  }

  updateConfiguration(newConfig: PatternConfiguration) {
    this._patternConfiguration.next(newConfig);
  }

  private getDefaultConfiguration(): PatternConfiguration {
    const nucleotideBases = this.externalDataManager.fetchNucleotideBases();
    const defaultBase = nucleotideBases[0];
    const defaultBases = new Array<string>(DEFAULT_SEQUENCE_LENGTH).fill(defaultBase);
    const defaultPhosphorothioate = new Array<boolean>(DEFAULT_SEQUENCE_LENGTH).fill(DEFAULT_PHOSPHOROTHIOATE);
    const defaultTerminalModification = '';
    const defaultComment = '';

    const defaultPatternConfiguration: PatternConfiguration = {
      [PATTERN_KEY.SS_BASES]: defaultBases,
      [PATTERN_KEY.AS_BASES]: defaultBases,
      [PATTERN_KEY.SS_PTO]: defaultPhosphorothioate,
      [PATTERN_KEY.AS_PTO]: defaultPhosphorothioate,
      [PATTERN_KEY.SS_3]: defaultTerminalModification,
      [PATTERN_KEY.SS_5]: defaultTerminalModification,
      [PATTERN_KEY.AS_3]: defaultTerminalModification,
      [PATTERN_KEY.AS_5]: defaultTerminalModification,
      [PATTERN_KEY.COMMENT]: defaultComment,
    }

    return defaultPatternConfiguration;
  }
}
