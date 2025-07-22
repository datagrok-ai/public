import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function app(): Promise<any> {
    return await grok.functions.call('EnamineChemSpace:App', {});
  }

  //Chemspace Samples
  export async function samplesPanel(smiles: string): Promise<any> {
    return await grok.functions.call('EnamineChemSpace:SamplesPanel', { smiles });
  }

  //Creates search panel
  export async function createSearchPanel(): Promise<any> {
    return await grok.functions.call('EnamineChemSpace:CreateSearchPanel', {});
  }

  //Chemspace Prices
  export async function pricesPanel(id: string): Promise<any> {
    return await grok.functions.call('EnamineChemSpace:PricesPanel', { id });
  }

  //Perform query with multipart form data
  export async function queryMultipart(path: string, formParamsStr: string, paramsStr: string): Promise<any> {
    return await grok.functions.call('EnamineChemSpace:QueryMultipart', { path, formParamsStr, paramsStr });
  }

  export async function getChemspaceIds(shipToCountry: string): Promise<any> {
    return await grok.functions.call('EnamineChemSpace:GetChemspaceIds', { shipToCountry });
  }

  export async function getChemspacePrices(data: DG.DataFrame, shipToCountry: string): Promise<any> {
    return await grok.functions.call('EnamineChemSpace:GetChemspacePrices', { data, shipToCountry });
  }
}
