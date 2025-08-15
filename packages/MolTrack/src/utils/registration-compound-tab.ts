import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { fetchCompoundProperties } from '../package';
import { MolTrackProp, Scope } from './constants';
import { EntityBaseView } from './registration-entity-base';


export class RegistrationCompoundView extends EntityBaseView {
  protected scope: Scope = Scope.COMPOUNDS;
  constructor(buildUI: boolean = true) {
    super(buildUI);
    this.view.name = 'Register a compound';
  }
}
