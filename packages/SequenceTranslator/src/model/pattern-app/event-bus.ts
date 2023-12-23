/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';

export class EventBus {
  private _antisenseStrandVisible$ = new rxjs.BehaviorSubject<boolean>(true);

  get antisenseStrandVisible$(): rxjs.Observable<boolean> {
    return this._antisenseStrandVisible$.asObservable();
  }

  toggleAntisenseStrand(value: boolean) {
    this._antisenseStrandVisible$.next(value);
  }
}
