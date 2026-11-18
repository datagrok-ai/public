import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {modelsFor} from '../catalog/applicable-tables';
import {MAX_MODELS} from '../constants';
import {forgeDb, ModelRow} from '../generated/db';
import {openApplyDialog} from './apply-model-dialog';
import {CARD_COLUMNS, ForgeModelHandler, forgeModelHandler} from './model-handler';
import {TrainView} from './train-view';
import '../../css/forge.css';

/** The `model` columns of a line: the card's, which its menu commands and the model's tooltip read. */
type PaneModelRow = Pick<ModelRow, typeof CARD_COLUMNS[number] | 'id' | 'author_id' | 'created_on'>;

/** **Models**, in the look of the platform's own **Models** pane: the models that fit [table], newest first, each a
 * line of the run icon (apply), the name as a link (the model's tooltip, menu and drag) and the menu arrow; then
 * **Train**.
 * Built when the pane opens; it holds no subscription, as the platform never detaches it. */
export function tableModelsPane(table: DG.DataFrame): HTMLElement {
  return ui.wait(async () => {
    const rows: PaneModelRow[] = await forgeDb.models.query().select(...CARD_COLUMNS).orderBy('created_on', true)
      .top(MAX_MODELS);
    const lines = modelsFor(rows, table).map((model) => modelLine(model, table));
    const train = ui.bigButton('Train', () => void TrainView.openFor(table), 'Train a new model on this table.');
    return ui.divV([lines.length === 0 ? ui.divText('No model fits this table.') : ui.div(lines, 'd4-flex-col'),
      ui.div([train], 'forge-model-train')]);
  });
}

function modelLine(model: PaneModelRow, table: DG.DataFrame): HTMLElement {
  const row = ForgeModelHandler.rowOf(model);
  const run = ui.iconSvg('model-run', () => void openApplyDialog(table,
    {modelId: row.id, preferredTable: table, switchToTable: true}), 'Run prediction');
  run.classList.add('d4-entity-markup-icon');
  // The platform's menu of the row: its commands and Forge's Apply... / Download.
  const actions = ui.contextActions(row);
  const name = ui.span([ui.label(row.displayName)], 'd4-link-label forge-model-name');
  // Forced: the platform drops an unforced change for 2 s after an accordion pane is expanded.
  name.addEventListener('click', () => grok.shell.setCurrentObject(row, false, true));
  // A ui.bind of a JS row gets only the platform's generic menu (Properties...): forward to the row's own menu.
  name.addEventListener('contextmenu', (e) => {
    e.preventDefault();
    e.stopPropagation();
    actions.dispatchEvent(new MouseEvent('click', {bubbles: true, clientX: e.clientX, clientY: e.clientY}));
  });
  ui.tooltip.bind(name, () => forgeModelHandler.renderTooltip(row));
  ui.makeDraggable(name, {getDragObject: () => row});
  // The built-in pane's row (d4.css): the 12 px icon hangs left of the name.
  const line = ui.div([run, name, actions], 'forge-model-line');
  line.dataset['wrapper'] = 'true';
  return ui.span([line], 'd4-entity-markup-row');
}
