// The MOEAD manager

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MoeadOptions, DEFAULT_SETTINGS, Validation, OptResult, InputOptions, Func, OPT_TYPE} from './defs';
import {OptimizeManager} from './optimize-manager';
import {OptimizationView} from '../optimization-view';

import {Moead} from './moead';
import {Visualizer} from './visualizer';

enum LIMITS {
  MIN_WEIGHTS = 1,
  MIN_GENS = 1,
  MIN_NEIG = 1,
  MIN_MUT = 0,
  MAX_MUT = 1,
};

export class MoeadManager extends OptimizeManager {
  private methodOpts: MoeadOptions;

  constructor(parent: OptimizationView) {
    super(parent);
    this.methodOpts = DEFAULT_SETTINGS;
  }

  public areSettingsValid(): Validation {
    for (const [key, value] of Object.entries(this.methodOpts)) {
      if ((value === null) || (value === undefined)) {
        return {
          res: false,
          msg: `Invalid value of "${key}" of the MEOA/D method settings`,
        };
      }
    }

    if (this.methodOpts.nWeights < LIMITS.MIN_WEIGHTS) {
      return {
        res: false,
        msg: `Invalid value of "samples" of the MEOA/D method settings`,
      };
    }

    if (this.methodOpts.generations < LIMITS.MIN_GENS) {
      return {
        res: false,
        msg: `Invalid value of "generations" of the MEOA/D method settings`,
      };
    }

    if (this.methodOpts.neighbors < LIMITS.MIN_NEIG) {
      return {
        res: false,
        msg: `Invalid value of "neighbors" of the MEOA/D method settings`,
      };
    }

    if ((this.methodOpts.mutationRate < LIMITS.MIN_MUT) || (this.methodOpts.mutationRate > LIMITS.MAX_MUT)) {
      return {
        res: false,
        msg: `Invalid value of "mutation rate" of the MEOA/D method settings`,
      };
    }

    return {res: true, msg: ''};
  }

  public getInputs(): DG.InputBase[] {
    return [
      ui.input.int('samples', {
        value: this.methodOpts.nWeights + 1,
        tooltipText: 'The number of points to be computed',
        nullable: false,
        min: LIMITS.MIN_WEIGHTS + 1,
        onValueChanged: (val) => {
          this.methodOpts.nWeights = val - 1;
          this.parent.updateApplicabilityState();
        },
      }),
      ui.input.int('generations', {
        value: this.methodOpts.generations,
        tooltipText: 'The total number of iterations',
        nullable: false,
        min: LIMITS.MIN_GENS,
        onValueChanged: (val) => {
          this.methodOpts.generations = val;
          this.parent.updateApplicabilityState();
        },
      }),
      ui.input.int('neighbors', {
        value: this.methodOpts.neighbors,
        tooltipText: 'The number of neighboring weight vectors considered for mating and replacement',
        nullable: false,
        min: LIMITS.MIN_NEIG,
        onValueChanged: (val) => {
          this.methodOpts.neighbors = val;
          this.parent.updateApplicabilityState();
        },
      }),
      ui.input.float('mutation rate', {
        value: this.methodOpts.mutationRate,
        tooltipText: 'The probability of mutating a variable in a solution',
        nullable: false,
        min: LIMITS.MIN_MUT,
        max: LIMITS.MAX_MUT,
        onValueChanged: (val) => {
          this.methodOpts.mutationRate = val;
          this.parent.updateApplicabilityState();
        },
      }),
    ];
  }

  public async perform(func: Func, inputOpts: InputOptions, outputDim: number, pi: DG.ProgressIndicator): OptResult {
    const moead = new Moead(
      func,
      inputOpts,
      outputDim,
      this.methodOpts,
      pi,
    );

    const results = await moead.perform();

    return results;
  };

  public visualize(view: DG.TableView, inputDim: number, outputDim: number, type: OPT_TYPE): DG.Viewer[] {
    const viz = new Visualizer(view, inputDim, outputDim, type);
    return viz.visualize();
  }
};
