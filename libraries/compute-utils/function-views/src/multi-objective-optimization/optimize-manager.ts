import * as DG from 'datagrok-api/dg';

import {Validation} from './defs';
import {OptimizationView} from '../optimization-view';

export abstract class OptimizeManager {
  protected parent: OptimizationView;
  constructor(parent: OptimizationView) {
    this.parent = parent;
  };

  public abstract areSettingsValid(): Validation;
  public abstract getInputs(): DG.InputBase[]
}
