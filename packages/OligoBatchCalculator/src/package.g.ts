import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: getUnits
//output: list<string> result
export function getUnits() {
  return PackageFunctions.getUnits();
}

//name: opticalDensity
//input: string sequence 
//input: double amount 
//input: string outputUnits { choices: OligoBatchCalculator:getUnits }
//input: dynamic extCoefsObj 
//output: double result
export async function opticalDensity(sequence: string, amount: number, outputUnits: string, extCoefsObj: any) {
  return PackageFunctions.opticalDensity(sequence, amount, outputUnits, extCoefsObj);
}

//name: nMole
//input: string sequence 
//input: double amount 
//input: string outputUnits { choices: OligoBatchCalculator:getUnits }
//input: dynamic extinctionCoefficientsObj 
//input: dynamic weightsObj 
//output: double result
export async function nMole(sequence: string, amount: number, outputUnits: string, extinctionCoefficientsObj: any, weightsObj: any) {
  return PackageFunctions.nMole(sequence, amount, outputUnits, extinctionCoefficientsObj, weightsObj);
}

//name: molecularMass
//input: string sequence 
//input: double amount 
//input: string outputUnits { choices: OligoBatchCalculator:getUnits }
//output: double result
export async function molecularMass(sequence: string, amount: number, outputUnits: string) {
  return PackageFunctions.molecularMass(sequence, amount, outputUnits);
}

//name: molecularWeight
//input: string sequence 
//input: string additionalWeightsObj 
//output: double result
export function molecularWeight(sequence: string, additionalWeightsObj?: any) {
  return PackageFunctions.molecularWeight(sequence, additionalWeightsObj);
}

//name: Oligo Batch Calculator
//tags: app
//meta.browsePath: Peptides | Oligo Toolkit
export async function OligoBatchCalculatorApp() {
  return PackageFunctions.OligoBatchCalculatorApp();
}
