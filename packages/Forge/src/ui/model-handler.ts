import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {featureFit} from '../catalog/applicable-tables';
import {MODEL_TYPE} from '../constants';
import {forgeDb, ModelRow} from '../generated/db';
import {modelAccordion, modelIcon} from './model-panes';

const HANDLER_NAME = 'Forge model handler';
const CARD_FEATURES = 5;

/** Forge's presentation of `forge.model` rows wherever the platform shows them: the catalog's context panel, the
 * Domains view, the Browse tree, the Predicted by pane. */
export class ForgeModelHandler extends DG.DomainObjectHandler {
  constructor() {
    super(MODEL_TYPE);
  }

  get name(): string {
    return HANDLER_NAME;
  }

  /** Registers the handler once per page, unless [registered] (the handlers' names) has it already (the package-test
   * bundle finds the one its package registered); true when it registered it now. */
  static registerOnce(registered: string[]): boolean {
    if (registered.includes(HANDLER_NAME))
      return false;
    DG.ObjectHandler.register(forgeModelHandler);
    return true;
  }

  /** [values] of one model, such as a catalog row, as a `DomainRow` without a round trip. */
  static rowOf(values: Pick<ModelRow, 'id' | 'name'> & Partial<ModelRow>): DG.DomainRow {
    // The platform parses the system dates from text.
    const {created_on: created, updated_on: updated, ...rest} = values;
    return forgeModelHandler.rowFrom({...rest, ...(created === undefined ? {} : {created_on: created.toISOString()}),
      ...(updated === undefined ? {} : {updated_on: updated.toISOString()})});
  }

  /** The model row behind [x]: a row, a semantic value or the platform's own object, as commands receive it; null
   * for anything else. */
  modelOf(x: unknown): DG.DomainRow | null {
    const isPlatformObject = typeof x === 'object' && x !== null && !('dart' in x);
    return this.rowOf(isPlatformObject ? DG.toJs(x) : x);
  }

  renderIcon(): HTMLElement {
    return modelIcon();
  }

  renderMarkup(x: unknown): HTMLElement {
    return ui.span([this.renderIcon(), ui.label(this.rowOrThrow(x).displayName)]);
  }

  renderTooltip(x: unknown): HTMLElement {
    const row = this.rowOrThrow(x);
    const v = row.values;
    return ui.divV([ui.label(row.displayName), ui.tableFromMap({'Method': v.engine_name, 'Task': v.task,
      'Target': v.target_name, 'Training rows': v.row_count ?? ''})]);
  }

  /** The built-in tool's card: name, the open tables it fits, target, features, method and creation date. */
  renderCard(x: unknown): HTMLElement {
    const row = this.rowOrThrow(x);
    const v = row.values;
    const {names: features, tables} = featureFit({features: v.features, options: v.options}, grok.shell.tables);
    const lines: HTMLElement[] = [ui.label(row.displayName, 'grok-gallery-grid-item-title')];
    if (tables.length > 0)
      lines.push(ui.label(`Applicable to ${tables.map((t) => t.name).join(', ')}`));
    lines.push(ui.label(`Predict ${v.target_name}`));
    if (features.length > 0) {
      const more = features.length > CARD_FEATURES ? '...' : '';
      lines.push(ui.label(`by ${features.slice(0, CARD_FEATURES).join(', ')}${more}`));
    }
    lines.push(ui.label(`using ${v.engine_name}`));
    const created = row.createdOn;
    if (created !== null)
      lines.push(ui.label(`Created on ${created.format('YYYY-MM-DD')}`, 'grok-gallery-grid-item-date'));
    return ui.bind(row, ui.divV(lines, 'd4-gallery-item'));
  }

  /** The model's accordion, from its full row read once. */
  renderProperties(x: unknown): HTMLElement {
    const row = this.rowOrThrow(x);
    return ui.wait(async () => {
      const model: ModelRow | null = await forgeDb.models.get(row.id);
      return model === null ? ui.divText('The model is no longer available.') : modelAccordion(row, model).root;
    });
  }
}

/** The one instance: registered with the platform and used by Forge's own code. */
export const forgeModelHandler = new ForgeModelHandler();
