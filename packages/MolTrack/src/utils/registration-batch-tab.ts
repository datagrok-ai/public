import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { EntityBaseView } from './registration-entity-base';
import { MolTrackProp, Scope } from './constants';
import { fetchBatchProperties } from '../package';

export class RegistrationBatchView extends EntityBaseView {
  protected scope: Scope = Scope.BATCHES;
  constructor() {
    super();
    this.view.name = 'Register a batch';
  }

  async getMolTrackDGProperties(): Promise<DG.Property[]> {
    try {
      const propsJson: MolTrackProp[] = (JSON.parse(await fetchBatchProperties())['properties'] as MolTrackProp[])
        .filter((p) => p.name.toLowerCase() !== 'corporate_batch_id');

      return propsJson.map(this.convertToDGProperty);
    } catch (err) {
      grok.shell.error(`Failed to fetch properties: ${err}`);
      return [];
    }
  }
}
