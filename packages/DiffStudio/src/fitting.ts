// In-worker parameters optimization features for Diff Studio models

import {IVP} from '@datagrok/diff-studio-tools';

/** Parameters optimization view for Diff Studio models */
export class Optimizer {
  /** Check applicability */
  static isApplicable(ivp: IVP): boolean {
    return (ivp.loop === null) && (ivp.updates === null) && (ivp.outputs === null);
  }
}; // Optimizer
