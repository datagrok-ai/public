import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Scripts {
  export async function sampleAB(A: number, B: number): Promise<number> {
    return await grok.functions.call('Python:SampleAB', { A, B });
  }

  //Calculates number of cells in the table
  export async function samplePandasInput(table: DG.DataFrame): Promise<number> {
    return await grok.functions.call('Python:SamplePandasInput', { table });
  }

  //Produces dataframe with an integer sequence
  export async function samplePandasOutput(size: number): Promise<DG.DataFrame> {
    return await grok.functions.call('Python:SamplePandasOutput', { size });
  }

  //RDKit-based script.
  export async function sampleGasteigerCharges(mol: string, contours: number): Promise<any> {
    return await grok.functions.call('Python:SampleGasteigerCharges', { mol, contours });
  }

  //Trains sklearn linear regression model. Compatible with predictive modeling toolkit
  export async function sampleSklearnModel(data: DG.DataFrame, predict_column: string): Promise<any> {
    return await grok.functions.call('Python:SampleSklearnModel', { data, predict_column });
  }
}

export namespace Funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('Python:Info', {});
  }
}
