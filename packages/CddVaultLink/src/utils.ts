
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {queryExportStatus, queryExportResult, ExportStatus, ApiResponse, Batch, Project} from "./cdd-vault-api";

export const CDD_HOST = 'https://app.collaborativedrug.com/';

export async function getAsyncResultsAsDf(vaultId: number, exportResponse: ApiResponse<ExportStatus>, timeoutMinutes: number, sdf: boolean): Promise<DG.DataFrame> {

    const resultResponse = await getAsyncResults(vaultId, exportResponse, timeoutMinutes, sdf);

    let df = DG.DataFrame.create();
    if (!resultResponse)
        return df;
    if (sdf) {
        const dfs = await grok.functions.call('Chem:importSdf', { bytes: resultResponse.data });
        if (dfs.length)
            df = dfs[0];
    } else {
        if (resultResponse.data?.objects) {
          prepareDataForDf(resultResponse.data.objects as any[]);
          df = DG.DataFrame.fromObjects(resultResponse.data.objects)!;
        }
    }
    return df;

}

const EXCLUDE_FIELDS = ['udfs', 'source_files']; 

function prepareDataForDf(objects: any[]) {
  for (let i = 0; i < objects.length; i++) {
    EXCLUDE_FIELDS.forEach((key) => delete objects[i][key]);
    if (objects[i]['batches']) {
      objects[i]['molecule_batch_identifiers'] = (objects[i]['batches'] as Batch[]).map((it) => it.molecule_batch_identifier).filter((it) => it !== undefined);
      delete objects[i]['batches'];
    }
    if (objects[i]['projects'])
      objects[i]['projects'] = (objects[i]['projects'] as Project[]).map((it) => it.name);
    if (objects[i]['molecule_fields']) {
      Object.keys(objects[i]['molecule_fields']).forEach((key: string) => objects[i][key] = objects[i]['molecule_fields'][key]);
      delete objects[i]['molecule_fields'];
    }
  }
}

export async function getAsyncResults(vaultId: number, exportResponse: ApiResponse<ExportStatus>, timeoutMinutes: number, text: boolean): Promise<ApiResponse<any> | null> {

    if (exportResponse.error) {
      grok.shell.error(exportResponse.error);
      return null;
    }
  
    const exportId = exportResponse.data?.id;
    if (!exportId) {
      grok.shell.error('No export ID received');
      return null;
    }
  
    const timeoutMs = timeoutMinutes * 60 * 1000;
    const startTime = Date.now();
  
    while (Date.now() - startTime < timeoutMs) {
      const statusResponse = await queryExportStatus(vaultId, exportId);
      
      if (statusResponse.error) {
        grok.shell.error(statusResponse.error);
        return null;
      }
  
      const status = statusResponse.data?.status;
      if (status === "finished") {
        const resultResponse = await queryExportResult(vaultId, exportId, text);
        if (resultResponse.error) {
          grok.shell.error(resultResponse.error);
          return null;
        }  
        return resultResponse;
      }
  
      // Wait for 2 seconds before next check
      await new Promise(resolve => setTimeout(resolve, 2000));
    }
  
    grok.shell.error(`Export timed out after ${timeoutMinutes} minutes`);
    return null;
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

export async function reorderColummns(df: DG.DataFrame) {
  const colNames = df.columns.names();
  const firstColumns = ['id', 'name', 'smiles'];
  const newColOrder = [];
  for (const colName of firstColumns) {
    const index = colNames.indexOf(colName);
    if (index > -1) {
      colNames.splice(index, 1);
      newColOrder.push(colName);
    }
  }
  df.columns.setOrder(newColOrder.concat(colNames));
}

export function paramsStringFromObj(params: any): string {
    let paramsStr = '';
    const paramNames = Object.keys(params);
    for (let i = 0; i < paramNames.length; i++) {
        const paramVal = params[paramNames[i]];
        if (paramVal) {
            paramsStr += paramsStr === '' ? `?${paramNames[i]}=${paramVal}` : `&${paramNames[i]}=${paramVal}`;
        }
    }
    return paramsStr;
}

export function createObjectViewer(obj: any, title: string = 'Object Viewer', additionalHeaderEl?: HTMLElement): HTMLElement {
  // Helper function to determine if a value is a dictionary/object
  function isDictionary(value: any): boolean {
    return value !== null && 
           typeof value === 'object' && 
           !Array.isArray(value) && 
           !(value instanceof Date);
  }

  // Helper function to determine if a value is a simple array
  function isSimpleArray(value: any): boolean {
    return Array.isArray(value) && 
           value.every(item => 
             typeof item === 'string' || 
             typeof item === 'number' || 
             typeof item === 'boolean'
           );
  }

  // Helper function to get pane name for array items
  function getPaneName(item: any, index: number, parentName: string): string {
    if (item.name) return item.name;
    if (item.id) return item.id.toString();
    if (parentName === 'protocol_statistics' && item.readout_definition?.name)
      return item.readout_definition.name;
    return `Item ${index}`;
  }

  // Helper function to create a view for a specific value
  function createValueView(value: any, parentName: string): HTMLElement {
    if (value === null || value === undefined) {
      return ui.divText('null');
    }

    if (isSimpleArray(value)) {
      return ui.divText(value.join(', '));
    }

    if (Array.isArray(value)) {
      const accordion = ui.accordion();
      value.forEach((item, index) => {
        if (isDictionary(item) || Array.isArray(item)) {
          accordion.addPane(getPaneName(item, index, parentName), () => createValueView(item, parentName));
        } else {
          accordion.addPane(`Item ${index}`, () => ui.divText(String(item)));
        }
      });
      return accordion.root;
    }

    if (isDictionary(value)) {
      const simpleProperties: { [key: string]: any } = {};
      const complexProperties: { [key: string]: any } = {};

      // Separate simple and complex properties
      for (const [key, val] of Object.entries(value)) {
        if (isDictionary(val) || Array.isArray(val)) {
          complexProperties[key] = val;
        } else {
          simpleProperties[key] = val;
        }
      }

      // Create container for both simple properties table and complex properties accordion
      const container = ui.divV([]);

      // Add simple properties table if any exist
      if (Object.keys(simpleProperties).length > 0) {
        container.appendChild(ui.tableFromMap(simpleProperties));
      }

      // Add complex properties as nested accordions if any exist
      if (Object.keys(complexProperties).length > 0) {
        const accordion = ui.accordion();
        for (const [key, val] of Object.entries(complexProperties)) {
          accordion.addPane(key, () => createValueView(val, parentName));
        }
        container.appendChild(accordion.root);
      }

      return container;
    }

    if (value instanceof Date) {
      return ui.divText(value.toISOString());
    }

    return ui.divText(String(value));
  }

  // Create dictionaries for simple and complex properties
  const simpleProperties: { [key: string]: any } = {};
  const complexProperties: { [key: string]: any } = {};

  // Separate simple and complex properties
  for (const [key, value] of Object.entries(obj)) {
    if (isDictionary(value) || Array.isArray(value)) {
      complexProperties[key] = value;
    } else {
      simpleProperties[key] = value;
    }
  }

  const header = ui.divH([ui.h2(title)]);
  if (additionalHeaderEl)
    header.append(additionalHeaderEl);
  // Create the main container
  const container = ui.divV([header]);

  // Add simple properties table if any exist
  if (Object.keys(simpleProperties).length > 0) {
    container.appendChild(ui.tableFromMap(simpleProperties));
  }

  // Add complex properties as accordion if any exist
  if (Object.keys(complexProperties).length > 0) {
    const accordion = ui.accordion();
    for (const [key, value] of Object.entries(complexProperties)) {
      accordion.addPane(key, () => createValueView(value, key));
    }
    container.appendChild(accordion.root);
  }

  return container;
}