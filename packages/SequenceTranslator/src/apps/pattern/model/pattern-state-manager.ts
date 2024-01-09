/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import {STRANDS, TERMINI, SEQUENCE_LENGTH, DEFAULT_PHOSPHOROTHIOATE} from './const';
import {PatternAppDataManager} from './external-data-manager';
import {EventBus} from './event-bus';
import {StrandType, TerminalType} from './types';

export class PatternConfigurationManager {
  private _patternName = '';
  private _isAntisenseStrandActive = true;
  private _bases = {} as Record<StrandType, string[]>;
  private _phosphorothioateLinkages = {} as Record<StrandType, boolean[]>;
  private _terminalModifications = {} as Record<StrandType, Record<TerminalType, string>>;
  private _comment = '';
  private _modificationsWithNumericLabels = [] as string[];

  constructor(private eventBus: EventBus, private dataManager: PatternAppDataManager) {
    this.initializeBases();
    this.initializePhosphorothioate();
    this.initializeTerminalModification();
    this.initializeComment();
  }

  private initializeBases(): void {
    const defaultBase = this.fetchDefaultBase();
    STRANDS.forEach((strand) => {
      this._bases[strand] = new Array(SEQUENCE_LENGTH.DEFAULT).fill(defaultBase);
    });
  }

  private fetchDefaultBase(): string {
    return this.dataManager.fetchAvailableNucleotideBases()[0];
  }
  
  private initializePhosphorothioate(): void {
    STRANDS.forEach((strand) => {
      this._phosphorothioateLinkages[strand] = new Array(SEQUENCE_LENGTH.DEFAULT).fill(DEFAULT_PHOSPHOROTHIOATE);
    });
  }
  
  private initializeTerminalModification(): void {
    STRANDS.forEach((strand) => {
      this._terminalModifications[strand] = {} as Record<TerminalType, string>;
      TERMINI.forEach((terminus) => {
        this._terminalModifications[strand][terminus] = '';
      });
    });
  }

  private initializeComment(): void {
    // this.eventBus.commentChange$.subscribe((comment) => {
    //   this._comment = comment;
    // });
  }

  getBases(strand: StrandType): string[] {
    return [...this._bases[strand]];
  }

  getPhosphorothioate(strand: StrandType): boolean[] {
    return [...this._phosphorothioateLinkages[strand]];
  }

  getTerminalModification(strand: StrandType, terminal: TerminalType): string {
    return this._terminalModifications[strand][terminal];
  }

  getComment(): string {
    return this._comment;
  }

  setPatternName(patternName: string): void {
    grok.shell.info(`Pattern name changed to ${patternName}`);
    this._patternName = patternName;
  }
}
