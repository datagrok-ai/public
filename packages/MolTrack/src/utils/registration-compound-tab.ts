import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { fetchCompoundProperties } from '../package';
import { MolTrackProp, Scope } from './constants';
import { EntityBaseView } from './registration-entity-base';


export class RegistrationCompoundView extends EntityBaseView {
  protected scope: Scope = Scope.COMPOUNDS;
  constructor() {
    super();
    this.view.name = 'Register a compound';
  }

  async getMolTrackDGProperties(): Promise<DG.Property[]> {
    try {
      const propsJson: MolTrackProp[] = (JSON.parse(await fetchCompoundProperties()) as MolTrackProp[])
        .filter((p) => p.name.toLowerCase() !== 'corporate_compound_id');

      return propsJson.map(this.convertToDGProperty);
    } catch (err) {
      grok.shell.error(`Failed to fetch properties: ${err}`);
      return [];
    }
  }
}
