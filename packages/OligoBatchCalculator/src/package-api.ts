import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function getUnits(): Promise<any> {
    return await grok.functions.call('@datagrok/oligo-batch-calculator:GetUnits', {});
  }

  export async function opticalDensity(sequence: string, amount: number, outputUnits: string): Promise<any> {
    return await grok.functions.call('@datagrok/oligo-batch-calculator:OpticalDensity', { sequence, amount, outputUnits });
  }

  export async function nMole(sequence: string, amount: number, outputUnits: string): Promise<any> {
    return await grok.functions.call('@datagrok/oligo-batch-calculator:NMole', { sequence, amount, outputUnits });
  }

  export async function molecularMass(sequence: string, amount: number, outputUnits: string): Promise<any> {
    return await grok.functions.call('@datagrok/oligo-batch-calculator:MolecularMass', { sequence, amount, outputUnits });
  }

  export async function molecularWeight(sequence: string, additionalWeightsObj: string): Promise<any> {
    return await grok.functions.call('@datagrok/oligo-batch-calculator:MolecularWeight', { sequence, additionalWeightsObj });
  }

  export async function oligoBatchCalculatorApp(): Promise<any> {
    return await grok.functions.call('@datagrok/oligo-batch-calculator:OligoBatchCalculatorApp', {});
  }
}
