
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from "@datagrok-libraries/utils/src/u2";
import {MoleculeFieldSearch, getVaults, MoleculeQueryParams, queryMolecules, queryReadoutRows, Molecule, Batch, querySavedSearches, SavedSearch, querySavedSearchById, queryExportStatus, queryExportResult, ExportStatus, ApiResponse} from "./cdd-vault-api";
import { CDDVaultSearchType } from './constants';

export const CDD_HOST = 'https://app.collaborativedrug.com/';

export async function getAsyncResults(vaultId: number, exportResponse: ApiResponse<ExportStatus>, timeoutMinutes: number, sdf: boolean): Promise<DG.DataFrame> {

  if (exportResponse.error) {
    grok.shell.error(exportResponse.error);
    return DG.DataFrame.create();
  }

  const exportId = exportResponse.data?.id;
  if (!exportId) {
    grok.shell.error('No export ID received');
    return DG.DataFrame.create();
  }

  const timeoutMs = timeoutMinutes * 60 * 1000;
  const startTime = Date.now();

  while (Date.now() - startTime < timeoutMs) {
    const statusResponse = await queryExportStatus(vaultId, exportId);
    
    if (statusResponse.error) {
      grok.shell.error(statusResponse.error);
      return DG.DataFrame.create();
    }

    const status = statusResponse.data?.status;
    if (status === "finished") {
      const resultResponse = await queryExportResult(vaultId, exportId, sdf);
      if (resultResponse.error) {
        grok.shell.error(resultResponse.error);
        return DG.DataFrame.create();
      }

      let df = DG.DataFrame.create();
      if (sdf) {
        const dfs = await grok.functions.call('Chem:importSdf', {bytes: resultResponse.data});
        if (dfs.length)
            df = dfs[0];
      } else {
          if (resultResponse.data?.objects)
              df = DG.DataFrame.fromObjects(resultResponse.data.objects)!;
      }

      return df;
    }

    // Wait for 2 seconds before next check
    await new Promise(resolve => setTimeout(resolve, 2000));
  }

  grok.shell.error(`Export timed out after ${timeoutMinutes} minutes`);
  return DG.DataFrame.create();
}

export async function createLinksFromIds(vaultId: number, df: DG.DataFrame) {
    const idCol = df.col('id');
    if (idCol) {
        const linkIdsCol = DG.Column.string('id', df.rowCount).init((i) => {
            const id = idCol.get(i);
            return `[${id}](${`${CDD_HOST}vaults/${vaultId}/molecules/${id}/`})`;
        });
        df.columns.replace(idCol, linkIdsCol);
    }
}