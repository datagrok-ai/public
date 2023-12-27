/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';

export class EventBus {
  private _isAntisenseStrandVisible$ = new rxjs.Subject<boolean>();
  private _patternListUpdated$ = new rxjs.Subject<void>();
  private _patternLoadRequested$ = new rxjs.Subject<string>();
  private _patternDeletionRequested$ = new rxjs.Subject<string>();

  get isAntisenseStrandActive$(): rxjs.Observable<boolean> {
    return this._isAntisenseStrandVisible$.asObservable();
  }

  get requestLoadPattern$(): rxjs.Observable<string> {
    return this._patternLoadRequested$.asObservable();
  }

  get patternListUpdate$(): rxjs.Observable<void> {
    return this._patternListUpdated$.asObservable();
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
}
