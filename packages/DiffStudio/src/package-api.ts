import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:Info', {});
  }

  export async function init(): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:Init', {});
  }

  export async function solve(problem: any): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:Solve', { problem });
  }

  export async function solveEquations(problem: any, options: any): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:SolveEquations', { problem, options });
  }

  //Solver of ordinary differential equations systems
  export async function runDiffStudio(): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:RunDiffStudio', {});
  }

  //Interactive solver of ordinary differential equations (ODE)
  export async function runDiffStudioDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:RunDiffStudioDemo', {});
  }

  export async function ivpFileHandler(content: string): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:IvpFileHandler', { content });
  }

  export async function previewIvp(file: DG.FileInfo): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:PreviewIvp', { file });
  }

  export async function runDiffStudioTreeBrowser(treeNode: any, browsePanel: any): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:RunDiffStudioTreeBrowser', { treeNode, browsePanel });
  }

  //Ball flight simulation
  export async function ballFlight(dB: number, roB: number, v: number, a: number): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:BallFlight', { dB, roB, v, a });
  }

  //Return serialized initial value problem for ordinary differential equations
  export async function serializeEquations(problem: string): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:SerializeEquations', { problem });
  }

  //Perform ODEs serialization to JS-code
  export async function odesToCode(serialization: any): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:OdesToCode', { serialization });
  }

  //Solve initial value problem for ordinary differential equations
  export async function solveODE(problem: string): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:SolveODE', { problem });
  }

  //In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
  export async function pkPdNew(): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:PkPdNew', {});
  }

  //In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
  export async function demoSimPKPD(): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:DemoSimPKPD', {});
  }

  //Controlled fab-arm exchange mechanism simulation
  export async function bioreactor(): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:Bioreactor', {});
  }

  //In-browser simulation of controlled fab-arm exchange mechanism
  export async function demoBioreactor(): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:DemoBioreactor', {});
  }

  //Run model with Diff Studio UI
  export async function runModel(model: string, inputsTabDockRatio: number, graphsDockRatio: number): Promise<any> {
    return await grok.functions.call('@datagrok/diff-studio:RunModel', { model, inputsTabDockRatio, graphsDockRatio });
  }
}
