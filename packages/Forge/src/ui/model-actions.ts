import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {readModelFile} from '../catalog/model-file';
import {MODEL_TYPE} from '../constants';
import {deleteModel, ownBlob} from '../storage/model-store';
import {openApplyDialog} from './apply-model-dialog';
import {openEditModelDialog} from './edit-model-dialog';
import {forgeModelHandler} from './model-handler';
import {reportError} from './report-error';

export interface ModelAction {
  name: string;
  /** Without it, the command applies to every model. */
  isApplicable?: (row: DG.DomainRow) => boolean;
  /** [preferredTable]: the table the catalog's **Applicable to** holds, which the Apply dialog opens on. */
  run: (row: DG.DomainRow, preferredTable: DG.DataFrame | null) => Promise<void> | void;
}

const APPLY: ModelAction = {name: 'Apply...', run: (row, preferredTable) => openModelApply(row.id, preferredTable)};
const DOWNLOAD: ModelAction = {name: 'Download', isApplicable: (row) => ownBlob(row.values.blob) !== null,
  run: downloadModelFile};

/** Forge's commands wherever the platform shows a model (the Domains view, the Browse tree, the panel's title), next
 * to the platform's own Open, Edit..., Clone, Delete, Share..., History and Copy link. */
export const MODEL_ACTIONS: readonly ModelAction[] = [APPLY, DOWNLOAD];

/** The catalog grid's commands: also editing and deleting the model the Forge way (Name, Description and Tags; the row
 * with its file and applications), which the platform's own Edit... and Delete do differently. */
export const CATALOG_ACTIONS: readonly ModelAction[] = [
  APPLY,
  {name: 'Edit model...', run: (row) => openEditModelDialog(row.id)},
  DOWNLOAD,
  {name: 'Delete model', run: (row) => confirmDeleteModel(row.id, row.displayName)},
];

/** Registers {@link MODEL_ACTIONS} as the commands of `forge.model` rows; `_initForge` calls it once per page, with the
 * model handler's registration. */
export function registerModelActions(): void {
  for (const action of MODEL_ACTIONS) {
    grok.functions.registerParamFunc(action.name, MODEL_TYPE, (x: unknown) => {
      const row = forgeModelHandler.modelOf(x);
      if (row !== null)
        void runModelAction(action, row, null);
    }, (x: unknown) => {
      const row = forgeModelHandler.modelOf(x);
      return row !== null && (action.isApplicable?.(row) ?? true);
    });
  }
}

/** The Apply dialog for the model [modelId], from the catalog or a model's menu: on the current table, or outside a
 * table view on the first open one, and the table's view after applying. */
export async function openModelApply(modelId: string, preferredTable: DG.DataFrame | null): Promise<void> {
  await openApplyDialog(grok.shell.currentTable ?? grok.shell.tables[0] ?? null,
    {modelId, switchToTable: true, preferredTable});
}

/** The commands of a catalog row: the catalog's grid shows its own menu, without the platform's model commands. */
export function addModelItems(menu: DG.Menu, row: DG.DomainRow, preferredTable: DG.DataFrame | null): void {
  for (const action of CATALOG_ACTIONS.filter((a) => a.isApplicable?.(row) ?? true))
    menu.item(action.name, () => void runModelAction(action, row, preferredTable));
}

/** The trash icon's and **Delete model**'s confirmation; OK deletes the row, its file and its applications. */
export function confirmDeleteModel(id: string, name: string): void {
  ui.dialog('Delete model')
    .add(ui.divText(`Delete the model "${name}"?`))
    .onOK(async () => {
      try {
        await deleteModel(id);
      } catch (e) {
        reportError(e);
      }
    })
    .show();
}

async function downloadModelFile(row: DG.DomainRow): Promise<void> {
  const file = await readModelFile({name: row.displayName, blob: row.values.blob});
  // A Blob takes only an ArrayBuffer-backed array; the platform's bytes are typed with any buffer.
  DG.Utils.download(file.name, Uint8Array.from(file.bytes));
}

async function runModelAction(action: ModelAction, row: DG.DomainRow,
  preferredTable: DG.DataFrame | null): Promise<void> {
  try {
    await action.run(row, preferredTable);
  } catch (e) {
    reportError(e);
  }
}
