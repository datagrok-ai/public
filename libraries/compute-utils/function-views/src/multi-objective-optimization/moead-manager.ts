// The MOEAD manager

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MoeadOptions, DEFAULT_SETTINGS, Validation} from './defs';
import {OptimizeManager} from './optimize-manager';

enum LIMITS {
  MIN_WEIGHTS = 1,
  MIN_GENS = 1,
  MIN_NEIG = 1,
  MIN_MUT = 0,
  MAX_MUT = 1,
};

export class MoeadManager extends OptimizeManager {
  private methodOpts: MoeadOptions;

  constructor() {
    super();
    this.methodOpts = DEFAULT_SETTINGS;
  }

  public areSettingsValid(): Validation {
    for (const [key, value] of Object.entries(this.methodOpts)) {
      if ((value === null) || (value === undefined)) {
        return {
          res: false,
          msg: `Invalid value for "${key}" of the MEOA/D method settings`,
        };
      }
    }

    if (this.methodOpts.nWeights < LIMITS.MIN_WEIGHTS) {
      return {
        res: false,
        msg: `Invalid value for "weights" of the MEOA/D method settings`,
      };
    }

    if (this.methodOpts.generations < LIMITS.MIN_GENS) {
      return {
        res: false,
        msg: `Invalid value for "generations" of the MEOA/D method settings`,
      };
    }

    if (this.methodOpts.neighbors < LIMITS.MIN_NEIG) {
      return {
        res: false,
        msg: `Invalid value for "neighbors" of the MEOA/D method settings`,
      };
    }

    if ((this.methodOpts.mutationRate < LIMITS.MIN_MUT) || (this.methodOpts.mutationRate > LIMITS.MAX_MUT)) {
      return {
        res: false,
        msg: `Invalid value for "mutation rate" of the MEOA/D method settings`,
      };
    }

    return {res: true, msg: ''};
  }

  public getInputs(): DG.InputBase[] {
    return [
      ui.input.int('weights', {
        value: this.methodOpts.nWeights,
        tooltipText: 'The number of weight vectors',
        nullable: false,
        min: LIMITS.MIN_WEIGHTS,
        onValueChanged: (val) => this.methodOpts.nWeights = val,
      }),
      ui.input.int('generations', {
        value: this.methodOpts.generations,
        tooltipText: 'The total number of iterations',
        nullable: false,
        min: LIMITS.MIN_GENS,
        onValueChanged: (val) => this.methodOpts.generations = val,
      }),
      ui.input.int('neighbors', {
        value: this.methodOpts.neighbors,
        tooltipText: 'The number of neighboring weight vectors considered for mating and replacement',
        nullable: false,
        min: LIMITS.MIN_NEIG,
        onValueChanged: (val) => this.methodOpts.neighbors = val,
      }),
      ui.input.float('mutation rate', {
        value: this.methodOpts.mutationRate,
        tooltipText: 'The probability of mutating a variable in a solution',
        nullable: false,
        min: LIMITS.MIN_MUT,
        max: LIMITS.MAX_MUT,
        onValueChanged: (val) => this.methodOpts.mutationRate = val,
      }),
    ];
  }
};
