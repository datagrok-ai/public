import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function initDrugBank(): Promise<any> {
    return await grok.functions.call('@datagrok/drug-bank:InitDrugBank', {});
  }

  export async function drugBankSubstructureSearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('@datagrok/drug-bank:DrugBankSubstructureSearchPanel', { mol });
  }

  export async function drugBankSimilaritySearchPanel(mol: string): Promise<any> {
    return await grok.functions.call('@datagrok/drug-bank:DrugBankSimilaritySearchPanel', { mol });
  }

  export async function drugNameMolecule(id: string): Promise<any> {
    return await grok.functions.call('@datagrok/drug-bank:DrugNameMolecule', { id });
  }
}
