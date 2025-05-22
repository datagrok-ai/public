import * as DG from 'datagrok-api/dg';

import {Func, InputOptions, OptResult, Validation} from './defs';
import {OptimizationView} from '../optimization-view';
import {OPT_TYPE} from './defs';

export abstract class OptimizeManager {
  protected parent: OptimizationView;
  constructor(parent: OptimizationView) {
    this.parent = parent;
  };

  public abstract areSettingsValid(): Validation;
  public abstract getInputs(): DG.InputBase[];
  public abstract perform(func: Func, inputOpts: InputOptions, outputDim: number, pi?: DG.ProgressIndicator): OptResult;
  public abstract visualize(view: DG.TableView, inputDim: number, outputDim: number, type: OPT_TYPE): DG.Viewer[];
}
