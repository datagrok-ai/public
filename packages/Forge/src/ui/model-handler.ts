import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {featureSchemasOf} from '../apply/apply-model';
import {featureFit} from '../catalog/applicable-tables';
import {MODEL_TYPE} from '../constants';
import {forgeDb, ModelRow, ModelStorageMode} from '../generated/db';
import {tagsOf} from '../storage/model-fields';
import {modelAccordion, modelIcon} from './model-panes';
import {richTooltip, shortList} from './rich-tooltip';
import {STORAGE_CAPTIONS} from './save-model-dialog';
import '../../css/forge.css';

const HANDLER_NAME = 'Forge model handler';
const CARD_FEATURES = 5;
const DATE_FORMAT = 'YYYY-MM-DD';
/** The authors' logins by user id, read once per session for the tooltips; a failed read is tried again. */
const LOGINS = new Map<string, Promise<string>>();
type FeatureRow = Pick<ModelRow, 'features' | 'options'>;
/** The feature lists and options of the models whose rows came without them, by id, read once per session; a failed
 * read is tried again. */
const FEATURE_ROWS = new Map<string, Promise<FeatureRow | null>>();

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

  /** The name; what it predicts from how many features with which method, the task and the training rows, then the
   * features; in grey the creation date and author, the data storage and the tags. */
  renderTooltip(x: unknown): HTMLElement {
    const row = this.rowOrThrow(x);
    const v = row.values;
    const facts = ui.divText('');
    const list = ui.divText('');
    const fill = (features: string[]) => {
      const count = typeof v.feature_count === 'number' ? v.feature_count : features.length;
      const from = count > 0 ? ` from ${ForgeModelHandler.featureCount(count)}` : '';
      facts.textContent = [`Predicts ${v.target_name}${from} with ${v.engine_name}`, v.task,
        ...(typeof v.row_count === 'number' ? [`${v.row_count} rows`] : [])].join(' · ');
      list.textContent = shortList(features);
      ui.setDisplay(list, features.length > 0);
    };
    fill([]);
    ForgeModelHandler.withFeatures(row, (fit) => fill((featureSchemasOf(fit.features) ?? []).map((f) => f.name)));
    return richTooltip({title: row.displayName, main: [facts, list], minor: [ForgeModelHandler.minorLine(row)]});
  }

  /** The built-in tool's card in five lines of one height: the name (with an icon when an open table fits), target,
   * features, method and creation date; the rest is in the tooltip. */
  renderCard(x: unknown): HTMLElement {
    const row = this.rowOrThrow(x);
    const v = row.values;
    const title = ui.divH([ui.label(row.displayName, 'grok-gallery-grid-item-title')], 'forge-model-card-title');
    const count = typeof v.feature_count === 'number' && v.feature_count > 0 ? v.feature_count : 0;
    const by = ui.label(count > 0 ? `by ${ForgeModelHandler.featureCount(count)}` : '');
    ForgeModelHandler.withFeatures(row, (fit) => {
      const {names, tables} = featureFit(fit, grok.shell.tables);
      if (names.length > 0)
        by.textContent = `by ${names.slice(0, CARD_FEATURES).join(', ')}${names.length > CARD_FEATURES ? '...' : ''}`;
      if (tables.length > 0)
        title.append(ui.iconFA('table', null, `Applicable to ${tables.map((t) => t.name).join(', ')}`));
    });
    const created = row.createdOn;
    return ui.bind(row, ui.divV([title, ui.label(`Predict ${v.target_name}`), by, ui.label(`using ${v.engine_name}`),
      ui.label(created === null ? '' : `Created on ${created.format(DATE_FORMAT)}`, 'grok-gallery-grid-item-date'),
    ], 'd4-gallery-item forge-model-card'));
  }

  /** Calls [use] with the model's feature list and options: at once when [row] has them, otherwise when they are read
   * (the Domains view's rows come without the json columns), once per model and session. */
  private static withFeatures(row: DG.DomainRow, use: (fit: FeatureRow) => void): void {
    const v = row.values;
    if (v.features != null && v.options != null) {
      use({features: v.features, options: v.options});
      return;
    }
    let read = FEATURE_ROWS.get(row.id);
    if (read === undefined) {
      // A model deleted or not shared with the user is not found: the lines stay as they are.
      read = forgeDb.models.query().where('id', '=', row.id).select('features', 'options').first()
        .then((model): FeatureRow | null => model ?? null, () => {
          FEATURE_ROWS.delete(row.id);
          return null;
        });
      FEATURE_ROWS.set(row.id, read);
    }
    read.then((fit) => {
      if (fit !== null)
        use(fit);
    });
  }

  /** `Created <date> by <author> · Data storage: <mode> · Tags: <tags>`, the parts the row has; the author's login is
   * read once and filled in when it comes. */
  private static minorLine(row: DG.DomainRow): HTMLElement {
    const v = row.values;
    const created = row.createdOn;
    const storage: string | undefined = STORAGE_CAPTIONS[v.storage_mode as ModelStorageMode];
    const tags = tagsOf(typeof v.tags === 'string' ? v.tags : '');
    const text = (author: string | undefined) => [
      ...(created === null ? [] : [`Created ${created.format(DATE_FORMAT)}${author ? ` by ${author}` : ''}`]),
      ...(storage === undefined ? [] : [`Data storage: ${storage}`]),
      ...(tags.length === 0 ? [] : [`Tags: ${tags.join(', ')}`]),
    ].join(' · ');
    const authorId = row.authorId;
    const line = ui.divText(text(undefined), 'forge-tooltip-minor');
    if (authorId) {
      let login = LOGINS.get(authorId);
      if (login === undefined) {
        login = grok.dapi.users.find(authorId).then((user: DG.User) => user.login);
        LOGINS.set(authorId, login);
      }
      // A deleted author is not found: the line stays without a login.
      login.then((name) => {
        line.textContent = text(name);
      }, () => LOGINS.delete(authorId));
    }
    return line;
  }

  private static featureCount(count: number): string {
    return `${count} ${count === 1 ? 'feature' : 'features'}`;
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
