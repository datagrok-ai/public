import {Rows, Section, divV, link, signal, span} from '@datagrok-libraries/u2';
import type {DomainAppOptions} from '@datagrok-libraries/ems/src/dg';
import {MolTrackEms} from './app';
import {RegisterApp} from './register-app';
import {MolTrackTables} from './tables';

/** A new batch registers its compound with it, unless the structure is registered already under the
 * matching rule: then the batch is added to that compound, as the legacy registration did. The page shows the
 * compound the batch is on above the batch's own fields. */
export class BatchesApp extends RegisterApp {
  /** The compound the structure drawn for a new batch turned out to be, which the batch was moved onto. */
  private readonly _matched = signal<{id: string, corporateId: string} | null>(null);
  private _seen = '';

  constructor(options: DomainAppOptions) {
    super(options);
    const compounds = options.inline!.compound_id.table;
    // the compound's derivation has looked its hash up by the time it writes it
    this.effect(() => {
      const compound = this.inline.value.get('compound_id');
      compound?.rows.items.value;
      const draft = compound?.currentRow.value;
      const hash = draft?.hash_mol;
      const batch = this.entitySource.value?.currentRow.peek();
      if (typeof hash !== 'string' || batch === null || batch === undefined || batch.compound_id !== draft!.id ||
          `${draft!.id} ${hash}` === this._seen)
        return;
      this._seen = `${draft!.id} ${hash}`;
      const existing = MolTrackTables.lookedUp([hash]).get(hash);
      if (existing === undefined)
        return;
      this._matched.value = {id: existing.id,
        corporateId: MolTrackTables.corporateId(compounds, 'molregno', existing.molregno)};
      this.entitySource.peek()!.edit.peek()!.setValue(batch.id, 'compound_id', existing.id);
    });
    this.aroundForm((form, source) => {
      this._matched.value = null;
      const compoundOf = () => {
        source.rows.items.value;
        const id = source.currentRow.value?.compound_id;
        return typeof id === 'string' && id !== '' && !Rows.isDraft(id) ? id : null;
      };
      const section = MolTrackEms.compoundSection('Compound', compoundOf);
      const note = span('');
      const matched = divV([note, link('Use a different structure', () => {
        this._matched.value = null;
        this._seen = '';
        const row = source.currentRow.peek()!;
        source.edit.peek()!.setValue(row.id, 'compound_id', null);
      })], 'moltrack-ems-matched');
      section.add(matched);
      section.effect(() => {
        const match = this._matched.value;
        const row = source.currentRow.value;
        matched.hidden = match === null || row === null || !Rows.isDraft(row) || compoundOf() !== match.id;
        note.textContent = match === null ? '' :
          `This structure is registered as ${match.corporateId}: the batch is added to it`;
      });
      const batch = new Section({title: 'Batch', collapsible: false});
      form.root.before(section.root, batch.root);
      batch.add(form.root);
    });
  }
}
