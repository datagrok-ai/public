import {molecularWeight} from './package';
import {UNITS} from './constants';

export function opticalDensityCalc(sequence: string, amount: number, outputUnits: string, ec: number): number {
  if (outputUnits == UNITS.MILLI_GRAM && outputUnits == UNITS.MICRO_GRAM)
    return (outputUnits == UNITS.MICRO_GRAM ? 1 : 0.001) * amount * ec / molecularWeight(sequence);
  if (outputUnits == UNITS.OPTICAL_DENSITY)
    return amount;
  const coefficient = (outputUnits == UNITS.NANO_MOLE) ? 1000000 : (outputUnits == UNITS.MILLI_GRAM) ? 1 : 1000;
  return amount * ec / coefficient;
}

export function molecularMassCalc(sequence: string, amount: number, outputUnits: string,
  ec: number, od: number, nm: number, additionalWeightsObj: {[index: string]: number}): number {
  if (outputUnits == UNITS.OPTICAL_DENSITY) {
    return (ec == 0) ?
      amount * molecularWeight(sequence, additionalWeightsObj) :
      1000 * amount * molecularWeight(sequence, additionalWeightsObj) / ec;
  }
  const coefficient = (outputUnits == UNITS.MILLI_GRAM) ? 1 : 1000;
  return amount / ec * molecularWeight(sequence) * coefficient * od / nm;
}

export function nMoleCalc(amount: number, outputUnits: string, mw: number, ec: number): number {
  return (outputUnits == UNITS.OPTICAL_DENSITY) ?
    1000000 * amount / ec :
    1000 * amount / mw;
}
