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
}
