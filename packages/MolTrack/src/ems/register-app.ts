import {untracked} from '@datagrok-libraries/u2';
import {DomainApp} from '@datagrok-libraries/ems/src/dg';
import type {DomainForm} from '@datagrok-libraries/ems/src/dg';
import type {DomainSource} from '@datagrok-libraries/ems';
import {RegisterDialog} from './register-dialog';

/** Compounds and batches arrive through registration: Import runs the bulk registration, not a plain import. */
export class RegisterApp extends DomainApp {
  async importRows(): Promise<boolean> {
    const target = this.table.address.endsWith('.compound') ? 'Compounds' : 'Batches';
    return await RegisterDialog.open(undefined, target) !== null;
  }

  /** Builds the sections around each entity page's form; a run's scope owns what it builds, so they go with the
   * page they were built for. */
  protected aroundForm(build: (form: DomainForm, source: DomainSource) => void): void {
    this.effect(() => {
      const form = this.form.value;
      const source = this.entitySource.value;
      if (form !== null && source !== null)
        untracked(() => build(form, source));
    });
  }
}
