/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';

export class EventBus {
  private static _instance: EventBus;

  private constructor() { }

  static getInstance(): EventBus {
    if (!EventBus._instance) {
      EventBus._instance = new EventBus();
    }
    return EventBus._instance;
  }

  private _isAntisenseStrandVisible$ = new rxjs.Subject<boolean>();
  private _patternListUpdated$ = new rxjs.Subject<void>();
  private _patternLoadRequested$ = new rxjs.Subject<string>();
  private _patternSaveRequested$ = new rxjs.Subject<string>();
  private _patternDeletionRequested$ = new rxjs.Subject<string>();
  private _tableSelectionChanged$ = new rxjs.BehaviorSubject<DG.DataFrame | null>(null);
  private _nucleotideBaseChanged$ = new rxjs.BehaviorSubject<string>('');


  get isAntisenseStrandActive$(): rxjs.Observable<boolean> {
    return this._isAntisenseStrandVisible$.asObservable();
  }

  get requestLoadPattern$(): rxjs.Observable<string> {
    return this._patternLoadRequested$.asObservable();
  }

  get patternListUpdate$(): rxjs.Observable<void> {
    return this._patternListUpdated$.asObservable();
  }

  get tableSelectionChanged$(): rxjs.Observable<DG.DataFrame | null> {
    return this._tableSelectionChanged$.asObservable();
  }

  get nucleotideBaseChanged$(): rxjs.Observable<string> {
    return this._nucleotideBaseChanged$.asObservable();
  }

  getTableSelection(): DG.DataFrame | null {
    return this._tableSelectionChanged$.getValue();
  }

  toggleAntisenseStrand(isActive: boolean) {
    this._isAntisenseStrandVisible$.next(isActive);
  }

  updatePatternList() {
    this._patternListUpdated$.next();
  }

  requestPatternLoad(patternName: string) {
    this._patternLoadRequested$.next(patternName);
  }

  deletePattern(patternName: string) {
    this._patternDeletionRequested$.next(patternName);
  }

  requestPatternSave(patternName: string) {
    this._patternSaveRequested$.next(patternName);
  }

  selectTable(table: DG.DataFrame | null) {
    this._tableSelectionChanged$.next(table);
  }

  changeNucleotideBase(base: string) {
    this._nucleotideBaseChanged$.next(base);
  }
}
