/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {InventoryApp} from './inventory-app';
export * from './package.g';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Inventory
//description: Inventory over entity-mapped domain schemas — batch upsert by SKU, optimistic stock adjustments, per-department column security
//tags: app
//output: view result
export async function inventoryApp(): Promise<DG.ViewBase> {
  return await InventoryApp.run();
}

//name: setupInventoryDemo
//description: Column-security demo setup - verifies the department property schemas are registered and creates the Chemists and Procurement groups (see README for the grant step)
//output: string result
export async function setupInventoryDemo(): Promise<string> {
  const schemas = await grok.dapi.stickyMeta.getSchemas();
  for (const name of ['chemistry', 'procurement', 'quality'])
    if (!schemas.some((s) => s.name === name))
      throw new Error(`Property schema '${name}' not found — is the Inventory package deployed?`);
  const created: string[] = [];
  for (const name of ['Chemists', 'Procurement'])
    if (await grok.dapi.groups.filter(`shortName = "${name}"`).first() == null) {
      await grok.dapi.groups.createNew(name);
      created.push(name);
    }
  const summary = (created.length === 0 ? 'Groups already exist.'
    : `Created groups: ${created.join(', ')}.`) +
    ' Grant the department schemas to them as described in the Inventory README.';
  grok.shell.info(summary);
  return summary;
}
