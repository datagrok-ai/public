/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';

export class EventBus {
  private _antisenseStrandVisible$ = new rxjs.BehaviorSubject<boolean>(true);
  private _commentChange$ = new rxjs.BehaviorSubject<string>('');
  private _updatePatternList$ = new rxjs.Subject<void>();
  private _requestLoadPattern$ = new rxjs.BehaviorSubject<string>('');
  private _selectedUserChanged$ = new rxjs.BehaviorSubject<string>('');
  private _deletePattern$ = new rxjs.Subject<string>();

  get antisenseStrandVisible$(): rxjs.Observable<boolean> {
    return this._antisenseStrandVisible$.asObservable();
  }

  get commentChange$(): rxjs.Observable<string> {
    return this._commentChange$.asObservable();
  }

  get requestLoadPattern$(): rxjs.Observable<string> {
    return this._requestLoadPattern$.asObservable();
  }

  get patternListUpdate$(): rxjs.Observable<void> {
    return this._updatePatternList$.asObservable();
  }

  toggleAntisenseStrand(value: boolean) {
    this._antisenseStrandVisible$.next(value);
  }

  changeComment(value: string) {
    this._commentChange$.next(value);
  }

  updatePatternList() {
    this._updatePatternList$.next();
  }

  loadNewPattern(patternName: string) {
    this._requestLoadPattern$.next(patternName);
  }

  getUserNameChoice(): string {
    return this._selectedUserChanged$.getValue();
  }

  deletePattern(patternName: string) {
    this._deletePattern$.next(patternName);
  }
}
