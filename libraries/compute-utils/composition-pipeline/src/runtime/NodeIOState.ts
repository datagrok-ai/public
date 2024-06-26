import {BehaviorSubject, Observable, of} from 'rxjs';
import {switchMap, takeUntil} from 'rxjs/operators';
import {StateItemConfiguration} from '../config/PipelineConfiguration';
import {ItemName, InputState} from '../config/CommonTypes';
import {debuglog} from '../utils';
import {PipelineDriverState} from './PipelineDriverState';

export class NodeIOState<T = any> {
  private currentSource = new BehaviorSubject<Observable<T>>(of());
  private valueChanges = this.currentSource.pipe(
    switchMap((source) => source),
  );

  public value = new BehaviorSubject<T | undefined>(undefined);

  private setter = (x: T, _inputState?: InputState) => {
    this.currentSource.next(of(x));
  };

  constructor(
    public conf: StateItemConfiguration,
    public pipelineState: PipelineDriverState,
    public parentId: ItemName,
    public notifier?: Observable<true>,
  ) {
    this.valueChanges.pipe(
      takeUntil(this.pipelineState.closed),
    ).subscribe((x) => {
      debuglog(`state updated: ${this.parentId}/${this.conf.id}, new value ${x}`);
      this.value.next(x);
    });
  }

  linkState(source: Observable<T>, setter?: (x: T, inputState?: InputState) => void) {
    this.currentSource.next(source);
    if (setter)
      this.setter = setter;
  }

  getValue() {
    return this.value.value;
  }

  setValue(val: T, inputState?: InputState) {
    this.setter(val, inputState);
  }
}
