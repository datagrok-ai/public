import {DomainApp} from '@datagrok-libraries/ems/src/dg';
import type {DomainAppOptions} from '@datagrok-libraries/ems/src/dg';
import {MolTrackEms} from './app';

/** The batch page shows the compound the batch is of above its fields, once one is picked. */
export class BatchesApp extends DomainApp {
  constructor(options: DomainAppOptions) {
    super(options);
    MolTrackEms.aroundForm(this, (form, source) => form.root.before(MolTrackEms.compoundSection('Compound', () => {
      source.rows.items.value;
      const id = source.currentRow.value?.compound_id;
      return typeof id === 'string' && id !== '' ? id : null;
    }).root));
  }
}
