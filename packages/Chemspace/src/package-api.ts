import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function app(): Promise<any> {
    return await grok.functions.call('Chemspace:App', {});
  }

  //Chemspace Samples
  export async function samplesPanel(smiles: string): Promise<any> {
    return await grok.functions.call('Chemspace:SamplesPanel', { smiles });
  }

  //Creates search panel
  export async function createSearchPanel(): Promise<any> {
    return await grok.functions.call('Chemspace:CreateSearchPanel', {});
  }

  //Chemspace Prices
  export async function pricesPanel(id: string): Promise<any> {
    return await grok.functions.call('Chemspace:PricesPanel', { id });
  }

  //Perform query with multipart form data
  export async function queryMultipart(path: string, formParamsStr: string, paramsStr: string): Promise<any> {
    return await grok.functions.call('Chemspace:QueryMultipart', { path, formParamsStr, paramsStr });
  }

  export async function getChemspaceIds(shipToCountry: string): Promise<any> {
    return await grok.functions.call('Chemspace:GetChemspaceIds', { shipToCountry });
  }

  export async function getChemspacePrices(data: DG.DataFrame, shipToCountry: string): Promise<any> {
    return await grok.functions.call('Chemspace:GetChemspacePrices', { data, shipToCountry });
  }
}
