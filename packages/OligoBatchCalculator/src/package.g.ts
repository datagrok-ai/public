import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//output: list<string> result
export function getUnits() : string[] {
  return PackageFunctions.getUnits();
}

//input: string sequence 
//input: double amount 
//input: string outputUnits { choices: OligoBatchCalculator:getUnits }
//input: dynamic extCoefsObj 
//output: double result
export async function opticalDensity(sequence: string, amount: number, outputUnits: string, extCoefsObj: any) : Promise<number> {
  return await PackageFunctions.opticalDensity(sequence, amount, outputUnits, extCoefsObj);
}

//input: string sequence 
//input: double amount 
//input: string outputUnits { choices: OligoBatchCalculator:getUnits }
//input: dynamic extinctionCoefficientsObj 
//input: dynamic weightsObj 
//output: double result
export async function nMole(sequence: string, amount: number, outputUnits: string, extinctionCoefficientsObj: any, weightsObj: any) : Promise<number> {
  return await PackageFunctions.nMole(sequence, amount, outputUnits, extinctionCoefficientsObj, weightsObj);
}

//input: string sequence 
//input: double amount 
//input: string outputUnits { choices: OligoBatchCalculator:getUnits }
//output: double result
export async function molecularMass(sequence: string, amount: number, outputUnits: string) : Promise<number> {
  return await PackageFunctions.molecularMass(sequence, amount, outputUnits);
}

//input: string sequence 
//input: string additionalWeightsObj 
//output: double result
export function molecularWeight(sequence: string, additionalWeightsObj?: any) : number {
  return PackageFunctions.molecularWeight(sequence, additionalWeightsObj);
}

//name: Oligo Batch Calculator
//meta.browsePath: Peptides | Oligo Toolkit
//meta.role: app
export async function OligoBatchCalculatorApp() : Promise<void> {
  await PackageFunctions.OligoBatchCalculatorApp();
}
