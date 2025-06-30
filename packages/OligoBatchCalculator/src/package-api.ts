import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function getUnits(): Promise<any> {
    return await grok.functions.call('OligoBatchCalculator:GetUnits', {});
  }

  export async function opticalDensity(sequence: string, amount: number, outputUnits: string): Promise<any> {
    return await grok.functions.call('OligoBatchCalculator:OpticalDensity', { sequence, amount, outputUnits });
  }

  export async function nMole(sequence: string, amount: number, outputUnits: string): Promise<any> {
    return await grok.functions.call('OligoBatchCalculator:NMole', { sequence, amount, outputUnits });
  }

  export async function molecularMass(sequence: string, amount: number, outputUnits: string): Promise<any> {
    return await grok.functions.call('OligoBatchCalculator:MolecularMass', { sequence, amount, outputUnits });
  }

  export async function molecularWeight(sequence: string, additionalWeightsObj: string): Promise<any> {
    return await grok.functions.call('OligoBatchCalculator:MolecularWeight', { sequence, additionalWeightsObj });
  }

  export async function oligoBatchCalculatorApp(): Promise<any> {
    return await grok.functions.call('OligoBatchCalculator:OligoBatchCalculatorApp', {});
  }
}
