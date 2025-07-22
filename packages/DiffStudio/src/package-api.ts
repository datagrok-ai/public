import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('DiffStudio:Info', {});
  }

  export async function init(): Promise<any> {
    return await grok.functions.call('DiffStudio:Init', {});
  }

  export async function solve(problem: any): Promise<any> {
    return await grok.functions.call('DiffStudio:Solve', { problem });
  }

  export async function solveEquations(problem: any, options: any): Promise<any> {
    return await grok.functions.call('DiffStudio:SolveEquations', { problem, options });
  }

  //Solver of ordinary differential equations systems
  export async function runDiffStudio(): Promise<any> {
    return await grok.functions.call('DiffStudio:RunDiffStudio', {});
  }

  //Interactive solver of ordinary differential equations (ODE)
  export async function runDiffStudioDemo(): Promise<any> {
    return await grok.functions.call('DiffStudio:RunDiffStudioDemo', {});
  }

  export async function ivpFileHandler(content: string): Promise<any> {
    return await grok.functions.call('DiffStudio:IvpFileHandler', { content });
  }

  export async function previewIvp(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('DiffStudio:PreviewIvp', { file });
  }

  export async function runDiffStudioTreeBrowser(treeNode: any, browsePanel: any): Promise<any> {
    return await grok.functions.call('DiffStudio:RunDiffStudioTreeBrowser', { treeNode, browsePanel });
  }

  //Ball flight simulation
  export async function ballFlight(dB: number, roB: number, v: number, a: number): Promise<any> {
    return await grok.functions.call('DiffStudio:BallFlight', { dB, roB, v, a });
  }

  //Return serialized initial value problem for ordinary differential equations
  export async function serializeEquations(problem: string): Promise<any> {
    return await grok.functions.call('DiffStudio:SerializeEquations', { problem });
  }

  //Perform ODEs serialization to JS-code
  export async function odesToCode(serialization: any): Promise<any> {
    return await grok.functions.call('DiffStudio:OdesToCode', { serialization });
  }

  //Solve initial value problem for ordinary differential equations
  export async function solveODE(problem: string): Promise<any> {
    return await grok.functions.call('DiffStudio:SolveODE', { problem });
  }

  //In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
  export async function pkPdNew(): Promise<any> {
    return await grok.functions.call('DiffStudio:PkPdNew', {});
  }

  //In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
  export async function demoSimPKPD(): Promise<any> {
    return await grok.functions.call('DiffStudio:DemoSimPKPD', {});
  }

  //Controlled fab-arm exchange mechanism simulation
  export async function bioreactor(): Promise<any> {
    return await grok.functions.call('DiffStudio:Bioreactor', {});
  }

  //In-browser simulation of controlled fab-arm exchange mechanism
  export async function demoBioreactor(): Promise<any> {
    return await grok.functions.call('DiffStudio:DemoBioreactor', {});
  }

  //Run model with Diff Studio UI
  export async function runModel(model: string, inputsTabDockRatio: number, graphsDockRatio: number): Promise<any> {
    return await grok.functions.call('DiffStudio:RunModel', { model, inputsTabDockRatio, graphsDockRatio });
  }
}
