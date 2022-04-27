/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {StateService} from '../common/service-interfaces';

export class DefaultStateService implements StateService {
  name = 'Model';

  lastCall: DG.FuncCall | null = null;

  func: DG.Func | null = null;

  prepare() {
    this.lastCall = this.func!.prepare();
  }

  public onComputationError: rxjs.Subject<DG.FuncCall> = new rxjs.Subject();
  public onComputationCompleted: rxjs.Subject<DG.FuncCall> = new rxjs.Subject();

  /** The actual computation function. */
  async compute(call: DG.FuncCall): Promise<void> {await call.call();}

  /** Maps inputs to parameters, computes, and maps output parameters to the UI. */
  async run(): Promise<void> {
    this.lastCall = this.lastCall!.clone();

    try {
      await this.compute(this.lastCall);
      this.onComputationCompleted.next(this.lastCall);
    } catch (e) {
      this.onComputationError.next(this.lastCall);
    }
  }
}
