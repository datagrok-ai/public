/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';

export class EventBus {
  private _isAntisenseStrandVisible$ = new rxjs.BehaviorSubject<boolean>(true);
  private _commentChanged$ = new rxjs.BehaviorSubject<string>('');
  private _patternListUpdated$ = new rxjs.Subject<void>();
  private _patternLoadRequested$ = new rxjs.BehaviorSubject<string>('');
  private _userSelectionChanged$ = new rxjs.BehaviorSubject<string>('');
  private _patternDeletionRequested$ = new rxjs.Subject<string>();

  get antisenseStrandActive$(): rxjs.Observable<boolean> {
    return this._isAntisenseStrandVisible$.asObservable();
  }

  get commentChange$(): rxjs.Observable<string> {
    return this._commentChanged$.asObservable();
  }

  get requestLoadPattern$(): rxjs.Observable<string> {
    return this._patternLoadRequested$.asObservable();
  }

  get patternListUpdate$(): rxjs.Observable<void> {
    return this._patternListUpdated$.asObservable();
  }

  setAntisenseStrandVisibility(visible: boolean) {
    this._isAntisenseStrandVisible$.next(visible);
  }

  changeComment(newComment: string) {
    this._commentChanged$.next(newComment);
  }

  updatePatternList() {
    this._patternListUpdated$.next();
  }

  requestPatternLoad(patternName: string) {
    this._patternLoadRequested$.next(patternName);
  }

  getUserNameChoice(): string {
    return this._userSelectionChanged$.getValue();
  }

  deletePattern(patternName: string) {
    this._patternDeletionRequested$.next(patternName);
  }
}
