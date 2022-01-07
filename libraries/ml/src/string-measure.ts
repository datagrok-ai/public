import * as fl from 'fastest-levenshtein';
import {jaroWinkler} from 'jaro-winkler-typescript';

import {DistanceMetric} from '@datagrok-libraries/utils/src/type-declarations';

const AvailableMetrics: {[name: string]: DistanceMetric} = {
  'Levenshtein': fl.distance,
  'Jaro-Winkler': jaroWinkler,
};

export type KnownMetrics = keyof typeof AvailableMetrics;

/** Unified class implementing different string measures. */
export class StringMeasure {
  protected method: string;

  /**
   * Creates an instance of StringMeasure with .
   * @param {string} method Method to calculate distance between strings.
   * @memberof Measurer
   */
  constructor(method: KnownMetrics) {
    this.method = method as string;
  }

  /**
   * Returns custom string distance function specified.
   * @return {DistanceMetric} Callback of the measure chosen.
   */
  public getMeasure(): DistanceMetric {
    return AvailableMetrics[this.method];
  }

  /** Returns metric names available. */
  static get availableMeasures(): string[] {
    return Object.keys(AvailableMetrics);
  }
}
