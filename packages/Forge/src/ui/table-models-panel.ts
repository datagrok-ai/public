import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {modelsFor} from '../catalog/applicable-tables';
import {ModelActivity, modelActivity} from '../catalog/model-activity';
import {MAX_MODELS} from '../constants';
import {forgeDb, ModelRow} from '../generated/db';
import {tagsOf} from '../storage/model-fields';
import {openApplyDialog} from './apply-model-dialog';
import {CARD_COLUMNS, ForgeModelHandler} from './model-handler';
import {modelFacts} from './model-panes';
import {TrainView} from './train-view';
import '../../css/forge.css';

/** The `model` columns of a line: the card's, which its menu commands read, and the tooltip's training table. */
const PANE_COLUMNS = [...CARD_COLUMNS, 'dataset_name'] as const;
type PaneModelRow = Pick<ModelRow, typeof PANE_COLUMNS[number] | 'id' | 'author_id' | 'created_on' | 'updated_on'>;

/** **Forge models**, the built-in **Models** pane's look: the models that fit [table], newest first, each a line of
 * the run icon (apply), the name as a link (the model's tooltip, menu and drag) and the menu arrow; then **Train**.
 * Built when the pane opens; it holds no subscription, as the platform never detaches it. */
export function tableModelsPane(table: DG.DataFrame): HTMLElement {
  return ui.wait(async () => {
    const rows: PaneModelRow[] = await forgeDb.models.query().select(...PANE_COLUMNS).orderBy('created_on', true)
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
  // The last run is read once per line; a failed read is not kept.
  let lastRun: Promise<ModelActivity['lastRun']> | undefined;
  const lastRunOf = () => lastRun ??= modelActivity(model.id).then((activity) => activity.lastRun, (e) => {
    lastRun = undefined;
    throw e;
  });
  // Built inside ui.wait: the platform calls the builder in a timer and drops the whole tooltip on a synchronous
  // error (tooltip.dart), while ui.wait shows it.
  ui.tooltip.bind(name, () => ui.wait(async () => modelSummary(model, lastRunOf)));
  ui.makeDraggable(name, {getDragObject: () => row});
  // The built-in pane's row (d4.css): the 12 px icon hangs left of the name.
  const line = ui.div([run, name, actions], 'forge-model-line');
  line.dataset['wrapper'] = 'true';
  return ui.span([line], 'd4-entity-markup-row');
}

/** The built-in model tooltip: the name, then Created by, Created, Updated, Trained on, Last run, Tags, Inputs,
 * Outputs, Method and Applicable to (when an open table fits), from the pane's row. */
function modelSummary(model: PaneModelRow, lastRun: () => Promise<ModelActivity['lastRun']>): HTMLElement {
  const facts = modelFacts(model, lastRun);
  const details: {[caption: string]: string | HTMLElement} = {
    'Created by': facts.author,
    'Created': facts.created,
    'Updated': facts.updated,
    'Trained on': facts.table,
    'Last run': facts.lastRun,
    'Tags': tagsOf(model.tags).join(', '),
    'Inputs': facts.features.join(', '),
    'Outputs': model.target_name,
    'Method': model.engine_name,
  };
  if (facts.applicableTo !== '')
    details['Applicable to'] = facts.applicableTo;
  return ui.divV([ui.span([ui.label(model.name)], 'd4-link-label'), ui.tableFromMap(details)]);
}
