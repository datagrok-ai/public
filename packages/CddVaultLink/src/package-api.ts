import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function cddVaultApp(path: string, filter: string): Promise<any> {
    return await grok.functions.call('CDDVaultLink:CddVaultApp', { path, filter });
  }

  export async function cddVaultAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('CDDVaultLink:CddVaultAppTreeBrowser', { treeNode, browseView });
  }

  export async function molColumnPropertyPanel(mol: string): Promise<any> {
    return await grok.functions.call('CDDVaultLink:MolColumnPropertyPanel', { mol });
  }

  export async function cddvaultSearchEditor(call: any): Promise<any> {
    return await grok.functions.call('CDDVaultLink:CDDVaultSearchEditor', { call });
  }

  export async function getVaultStats(vaultId: number, vaultName: string): Promise<any> {
    return await grok.functions.call('CDDVaultLink:GetVaultStats', { vaultId, vaultName });
  }

  export async function getVaults(): Promise<any> {
    return await grok.functions.call('CDDVaultLink:GetVaults', {});
  }

  export async function getMolecules(vaultId: number, moleculesIds: string): Promise<any> {
    return await grok.functions.call('CDDVaultLink:GetMolecules', { vaultId, moleculesIds });
  }

  export async function getMoleculesAsync(vaultId: number, moleculesIds: string, timeoutMinutes: number): Promise<any> {
    return await grok.functions.call('CDDVaultLink:GetMoleculesAsync', { vaultId, moleculesIds, timeoutMinutes });
  }

  export async function getProtocolsAsync(vaultId: number, timeoutMinutes: number): Promise<any> {
    return await grok.functions.call('CDDVaultLink:GetProtocolsAsync', { vaultId, timeoutMinutes });
  }

  export async function getCollectionsAsync(vaultId: number, timeoutMinutes: number): Promise<any> {
    return await grok.functions.call('CDDVaultLink:GetCollectionsAsync', { vaultId, timeoutMinutes });
  }

  export async function getSavedSearches(vaultId: number): Promise<any> {
    return await grok.functions.call('CDDVaultLink:GetSavedSearches', { vaultId });
  }

  export async function getSavedSearchResults(vaultId: number, searchId: number, timeoutMinutes: number): Promise<any> {
    return await grok.functions.call('CDDVaultLink:GetSavedSearchResults', { vaultId, searchId, timeoutMinutes });
  }

  export async function cDDVaultSearchAsync(vaultId: number, structure: string, structure_search_type: string, structure_similarity_threshold: number, protocol: number, run: number): Promise<any> {
    return await grok.functions.call('CDDVaultLink:CDDVaultSearchAsync', { vaultId, structure, structure_search_type, structure_similarity_threshold, protocol, run });
  }

  export async function cDDVaultSearch(vaultId: number, molecules: string, names: string, include_original_structures: boolean, only_ids: boolean, only_batch_ids: boolean, created_before: string, created_after: string, modified_before: string, modified_after: string, batch_created_before: string, batch_created_after: string, batch_field_before_name: string, batch_field_before_date: string, batch_field_after_name: string, batch_field_after_date: string, projects: string, data_sets: string, structure: string, structure_search_type: string, structure_similarity_threshold: number, inchikey: string, molecule_fields: any, batch_fields: any, fields_search: any): Promise<any> {
    return await grok.functions.call('CDDVaultLink:CDDVaultSearch', { vaultId, molecules, names, include_original_structures, only_ids, only_batch_ids, created_before, created_after, modified_before, modified_after, batch_created_before, batch_created_after, batch_field_before_name, batch_field_before_date, batch_field_after_name, batch_field_after_date, projects, data_sets, structure, structure_search_type, structure_similarity_threshold, inchikey, molecule_fields, batch_fields, fields_search });
  }

  export async function cDDVaultSearch2(vaultId: number, structure: string, structure_search_type: string, structure_similarity_threshold: number, protocol: number, run: number): Promise<any> {
    return await grok.functions.call('CDDVaultLink:CDDVaultSearch2', { vaultId, structure, structure_search_type, structure_similarity_threshold, protocol, run });
  }
}
