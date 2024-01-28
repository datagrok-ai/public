/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {StrandType, TerminalType} from './types';
import {NucleotideSequences, PhosphorothioateLinkageFlags, StrandTerminusModifications} from './types';

import * as rxjs from 'rxjs';
import {DefaultStateConfigurator} from './default-state-configurator';

export class EventBus {
  private _patternName$: rxjs.BehaviorSubject<string>;
  private _isAntisenseStrandVisible$: rxjs.BehaviorSubject<boolean>;
  private _nucleobases$: rxjs.BehaviorSubject<NucleotideSequences>;
  private _phosphorothioateLinkageFlags: rxjs.BehaviorSubject<PhosphorothioateLinkageFlags>;
  private _terminalModifications: rxjs.BehaviorSubject<StrandTerminusModifications>;
  private _comment$: rxjs.BehaviorSubject<string>;
  private _modificationsWithNumericLabels$: rxjs.BehaviorSubject<string[]>;

  private _patternListUpdated$ = new rxjs.Subject<void>();
  private _patternLoadRequested$ = new rxjs.Subject<string>();
  private _patternSaveRequested$ = new rxjs.Subject<string>();

  private _patternDeletionRequested$ = new rxjs.Subject<string>();
  private _tableSelection$ = new rxjs.BehaviorSubject<DG.DataFrame | null>(null);
  private _sequenceBaseChanged$ = new rxjs.Subject<string>();


  constructor(defaults: DefaultStateConfigurator) {
    this.initializeDefaultState(defaults);
  }

  private initializeDefaultState(defaults: DefaultStateConfigurator) {
    this._patternName$ = new rxjs.BehaviorSubject(defaults.getPatternName());
    this._isAntisenseStrandVisible$ = new rxjs.BehaviorSubject(defaults.getAntiSenseStrandVisibilityFlag());
    this._nucleobases$ = new rxjs.BehaviorSubject(defaults.getNucleobases());
    this._phosphorothioateLinkageFlags = new rxjs.BehaviorSubject(defaults.getPhosphorothioateLinkageFlags());
    this._terminalModifications = new rxjs.BehaviorSubject(defaults.getTerminusModifications());
    this._comment$ = new rxjs.BehaviorSubject(defaults.getComment());
    this._modificationsWithNumericLabels$ = new rxjs.BehaviorSubject(defaults.getModificationsWithNumericLabels());
  }

  getPatterName(): string {
    return this._patternName$.getValue();
  }

  setNewPatternName(patternName: string) {
    this._patternName$.next(patternName);
  }

  get antisenseStrandToggled$(): rxjs.Observable<boolean> {
    return this._isAntisenseStrandVisible$.asObservable();
  }

  toggleAntisenseStrand(isActive: boolean) {
    this._isAntisenseStrandVisible$.next(isActive);
  }

  antisenseStrandIsVisible(): boolean {
    return this._isAntisenseStrandVisible$.getValue();
  }

  changeSequenceBase(base: string) {
    this._sequenceBaseChanged$.next(base);
  }

  get patternLoadRequested$(): rxjs.Observable<string> {
    return this._patternLoadRequested$.asObservable();
  }

  requestPatternLoad(patternName: string) {
    this._patternLoadRequested$.next(patternName);
  }

  get patternListUpdated$(): rxjs.Observable<void> {
    return this._patternListUpdated$.asObservable();
  }

  updatePatternList() {
    this._patternListUpdated$.next();
  }

  get tableSelectionChanged$(): rxjs.Observable<DG.DataFrame | null> {
    return this._tableSelection$.asObservable();
  }

  selectTable(table: DG.DataFrame | null) {
    this._tableSelection$.next(table);
  }

  getTableSelection(): DG.DataFrame | null {
    return this._tableSelection$.getValue();
  }

  deletePattern(patternName: string) {
    this._patternDeletionRequested$.next(patternName);
  }

  requestPatternSave(patternName: string) {
    this._patternSaveRequested$.next(patternName);
  }
}
