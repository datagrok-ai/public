import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {featureSchemasOf} from '../apply/apply-model';
import {DATE_FORMAT, MODEL_TYPE} from '../constants';
import {forgeDb, ModelRow} from '../generated/db';
import {modelsChanged} from '../storage/model-store';
import {isShownVersion, modelAccordion, modelIcon, userOf} from './model-panes';
import {reportError} from './report-error';
import {richTooltip} from './rich-tooltip';
import '../../css/forge.css';

/** The `model` columns the markup, the card, their rich tooltip and the commands (Download reads blob) read. */
export const CARD_COLUMNS = ['name', 'features', 'options', 'target_name', 'engine_name', 'task', 'row_count',
  'storage_mode', 'tags', 'blob'] as const;
const HANDLER_NAME = 'Forge model handler';
const CARD_FEATURES = 5;
/** The tooltip's feature list is cut to this many characters. */
export const MAX_FEATURES_TEXT = 50;
const FOLLOW_DELAY_MS = 300;
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

  /** Keeps the context panel on the model it shows as that model is now, whoever changed it (the catalog, **Edit
   * model**, the API); `_initForge` subscribes once per page. */
  static followChanges(): rxjs.Subscription {
    // A burst of writes (tag chips) reads the model once.
    return DG.debounce(modelsChanged, FOLLOW_DELAY_MS).subscribe(() => void ForgeModelHandler.showChanged());
  }

  /** The model the context panel shows, read again: shown anew when it changed since its accordion was built, and the
   * panel emptied, as after the catalog's delete, when it is gone. */
  private static async showChanged(): Promise<void> {
    const shown = forgeModelHandler.modelOf(grok.shell.o);
    if (shown === null)
      return;
    try {
      const model = await forgeDb.models.query().where('id', '=', shown.id).select(...CARD_COLUMNS).first();
      if (forgeModelHandler.modelOf(grok.shell.o)?.id !== shown.id)
        return;
      // Forced: the platform takes a row of the same id for the object it already shows.
      if (model === null)
        grok.shell.setCurrentObject(null, true, true);
      else if (!isShownVersion(model))
        grok.shell.setCurrentObject(ForgeModelHandler.rowOf(model), true, true);
    } catch (e) {
      reportError(e);
    }
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

  /** The name; `Predicts <target> with <method> from <features>`; in grey `Created <date> by <author>`. */
  renderTooltip(x: unknown): HTMLElement {
    const row = this.rowOrThrow(x);
    const v = row.values;
    // A row made by hand may lack any of these.
    const [target, method] = [v.target_name, v.engine_name].map((s) => typeof s === 'string' ? s : '');
    const predicts = ui.divText(ForgeModelHandler.predictsText(target, method, []));
    ForgeModelHandler.withFeatures(row, (fit) => {
      predicts.textContent = ForgeModelHandler.predictsText(target, method,
        (featureSchemasOf(fit.features) ?? []).map((f) => f.name));
    });
    return richTooltip({title: row.displayName, main: [predicts], minor: [ForgeModelHandler.createdLine(row)]});
  }

  /** `Predicts <target> with <method> from <features>`, leaving out the parts that are empty; the features joined by
   * commas and cut to {@link MAX_FEATURES_TEXT} characters with a trailing `…`; empty when all three are. */
  static predictsText(target: string, method: string, features: string[]): string {
    if (target === '' && method === '' && features.length === 0)
      return '';
    const joined = features.join(', ');
    const from = joined.length > MAX_FEATURES_TEXT ?
      `${joined.slice(0, MAX_FEATURES_TEXT - 1).replace(/[\s,]+$/, '')}…` : joined;
    return ['Predicts', target, ...(method === '' ? [] : ['with', method]), ...(from === '' ? [] : ['from', from])]
      .filter((s) => s !== '').join(' ');
  }

  /** The built-in tool's card in five lines of one height: the name, target,
   * features, method and creation date; the rest is in the tooltip. */
  renderCard(x: unknown): HTMLElement {
    const row = this.rowOrThrow(x);
    return ui.bind(row, ForgeModelHandler.cardContent(row));
  }

  /** The card's five lines, not bound to the row: the catalog's gallery handles the clicks itself. */
  static cardContent(row: DG.DomainRow): HTMLElement {
    const v = row.values;
    const title = ui.divH([ui.label(row.displayName, 'grok-gallery-grid-item-title')], 'forge-model-card-title');
    const count = typeof v.feature_count === 'number' && v.feature_count > 0 ? v.feature_count : 0;
    const by = ui.label(count > 0 ? `by ${ForgeModelHandler.featureCount(count)}` : '');
    ForgeModelHandler.withFeatures(row, (fit) => {
      const names = (featureSchemasOf(fit.features) ?? []).map((f) => f.name);
      if (names.length > 0)
        by.textContent = `by ${names.slice(0, CARD_FEATURES).join(', ')}${names.length > CARD_FEATURES ? '...' : ''}`;
    });
    const created = row.createdOn;
    return ui.divV([title, ui.label(`Predict ${v.target_name}`), by, ui.label(`using ${v.engine_name}`),
      ui.label(created === null ? '' : `Created on ${created.format(DATE_FORMAT)}`, 'grok-gallery-grid-item-date'),
    ], 'd4-gallery-item forge-model-card');
  }

  /** The catalog's list row: the name as a link, as the platform's brief lists show entities (a label: it goes nowhere,
   * the gallery handles the clicks). */
  static listLine(row: DG.DomainRow): HTMLElement {
    return ui.label(row.displayName, 'd4-link-label');
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

  /** `Created <date> by <author>`, empty without a date; the author's login is read once and filled in when it
   * comes. */
  private static createdLine(row: DG.DomainRow): HTMLElement {
    const created = row.createdOn;
    const text = (author: string | undefined) =>
      created === null ? '' : `Created ${created.format(DATE_FORMAT)}${author ? ` by ${author}` : ''}`;
    const authorId = row.authorId;
    const line = ui.divText(text(undefined), 'forge-tooltip-minor');
    // A deleted author is not found: the line stays without a login.
    if (authorId) {
      userOf(authorId).then((user) => {
        line.textContent = text(user.login);
      }, () => {});
    }
    return line;
  }

  /** `1 feature`, `<n> features`. */
  static featureCount(count: number): string {
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
