import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Callback} from './callback-base';
import {CallbackAction} from '../solver-defs';
import {WARNING} from '../../ui-constants';

/** This callback proposes to abort computations if the maximum calculation time is exceeded */
export class DialogCallback extends Callback {
  private maxTime = 0;
  private startingTime = 0;
  private isDlgAllowed = true;
  private toAbortComputations = false;
  private isDlgShown = false;
  private dlg = ui.dialog(WARNING.TITLE);

  constructor(maxTime: number) {
    super();
    //this.dlg.show();
    this.maxTime = maxTime;
    this.startingTime = performance.now();

    this.dlg.getButton('CANCEL').hidden = true;
    //this.dlg.onClose.subscribe(() => this.isDlgShown = false);
  }

  public onIterationStart(): void {
    if (this.toAbortComputations)
      throw new CallbackAction('Computations aborted');

    if (this.isDlgAllowed) {
      if (performance.now() - this.startingTime > this.maxTime) {
        alert('TOO long!');
        console.log('DIALOG!');
        this.isDlgAllowed = false;
        this.isDlgShown = true;
        this.dlg.show();
      }
    }
  }

  public onComputationsCompleted(): void {
    if (this.isDlgShown)
      this.dlg.close();
  }
}
