import * as DG from 'datagrok-api/dg';

import {Validation} from './defs';

export abstract class OptimizeManager {
  constructor() {};

  public abstract areSettingsValid(): Validation;
  public abstract getInputs(): DG.InputBase[]
}
