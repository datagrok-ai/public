/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import {PATTERN_KEY, STRANDS, StrandType, TerminalType, TERMINAL_KEYS } from './const';
import {PatternConfiguration} from './types';
import {DEFAULT_SEQUENCE_LENGTH, DEFAULT_PHOSPHOROTHIOATE} from './const';
import {ExternalDataManager} from './external-data-manager';

export class PatternConfigurationManager {
  private _bases = {} as Record<StrandType, string[]>;
  private _phosphorothioate = {} as Record<StrandType, boolean[]>;
  private _terminalModification = {} as Record<StrandType, Record<TerminalType, string>>;
  private _comment: string;

  constructor() {
    const defaultBase = this.fetchDefaultBase();
    STRANDS.forEach((strand) => {
      this._bases[strand] = new Array(DEFAULT_SEQUENCE_LENGTH).fill(defaultBase);
    });
    STRANDS.forEach((strand) => {
      this._phosphorothioate[strand] = new Array(DEFAULT_SEQUENCE_LENGTH).fill(DEFAULT_PHOSPHOROTHIOATE);
    });
    STRANDS.forEach((strand) => {
      this._terminalModification[strand] = {} as Record<TerminalType, string>;
      TERMINAL_KEYS.forEach((terminal) => {
        this._terminalModification[strand][terminal] = '';
      });
    });
    this._comment = '';
  }

  private fetchDefaultBase(): string {
    return new ExternalDataManager().fetchNucleotideBases()[0];
  }

  getBases(strand: StrandType): string[] {
    return [...this._bases[strand]];
  }

  getPhosphorothioate(strand: StrandType): boolean[] {
    return [...this._phosphorothioate[strand]];
  }

  getTerminalModification(strand: StrandType, terminal: TerminalType): string {
    return this._terminalModification[strand][terminal];
  }

  getComment(): string {
    return this._comment;
  }
}
