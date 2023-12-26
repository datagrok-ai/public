/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import {PATTERN_KEY, STRANDS, StrandType, TerminalType, TERMINAL_KEYS } from './const';
import {PatternConfiguration} from './types';
import {DEFAULT_SEQUENCE_LENGTH, DEFAULT_PHOSPHOROTHIOATE} from './const';
import {PatternAppDataManager} from './external-data-manager';
import {EventBus} from './event-bus';

export class PatternConfigurationManager {
  private _bases = {} as Record<StrandType, string[]>;
  private _phosphorothioate = {} as Record<StrandType, boolean[]>;
  private _terminalModification = {} as Record<StrandType, Record<TerminalType, string>>;
  private _comment = '';

  constructor(private eventBus: EventBus, private dataManager: PatternAppDataManager) {
    this.initializeBases();
    this.initializePhosphorothioate();
    this.initializeTerminalModification();
    this.initializeComment();
  }

  private initializeBases(): void {
    const defaultBase = this.fetchDefaultBase();
    STRANDS.forEach((strand) => {
      this._bases[strand] = new Array(DEFAULT_SEQUENCE_LENGTH).fill(defaultBase);
    });
  }

  private fetchDefaultBase(): string {
    return this.dataManager.fetchNucleotideBases()[0];
  }
  
  private initializePhosphorothioate(): void {
    STRANDS.forEach((strand) => {
      this._phosphorothioate[strand] = new Array(DEFAULT_SEQUENCE_LENGTH).fill(DEFAULT_PHOSPHOROTHIOATE);
    });
  }
  
  private initializeTerminalModification(): void {
    STRANDS.forEach((strand) => {
      this._terminalModification[strand] = {} as Record<TerminalType, string>;
      TERMINAL_KEYS.forEach((terminal) => {
        this._terminalModification[strand][terminal] = '';
      });
    });
  }

  private initializeComment(): void {
    this.eventBus.commentChange$.subscribe((comment) => {
      this._comment = comment;
    });
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
