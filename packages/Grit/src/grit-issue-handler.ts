import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {gritDb, IssueRow} from './generated/db';

/** semType (= entity type) of `grit.issue` rows. */
const ISSUE_TYPE = 'grit.issue';

/** A Grit issue handle: `<PROJECTKEY>-<number>`, e.g. `GRITEST-1` (per-project numbering). */
const ISSUE_HANDLE_RE = /^([A-Za-z][A-Za-z0-9_]*)-(\d+)$/;

const statusColors: {[key: string]: string} = {
  'open': '#1f8fff', 'in progress': '#e8912d', 'resolved': '#2a9d3a', 'closed': '#8b95a1',
};
const priorityColors: {[key: string]: string} = {
  'low': '#8b95a1', 'medium': '#1f8fff', 'high': '#e8912d', 'critical': '#d9463f',
};

/** Injects the static badge presentation once; the per-badge color stays inline. */
function ensureBadgeStyles(): void {
  if (document.getElementById('grit-badge-styles') != null)
    return;
  const style = document.createElement('style');
  style.id = 'grit-badge-styles';
  style.textContent =
    '.grit-badge{color:#fff;border-radius:4px;padding:1px 6px;margin-left:4px;font-size:11px;}';
  document.head.appendChild(style);
}

function badge(text: string, color: string): HTMLElement {
  ensureBadgeStyles();
  const b = ui.span([text], 'grit-badge');
  b.style.backgroundColor = color;
  return b;
}

/** Custom presentation of `grit.issue` rows — badges on cards/tooltips/markup and a
 * grid tweak; everything NOT overridden here (the Entity View, the context panel,
 * commands, pickers) is the platform default via {@link DG.DomainObjectHandler}'s
 * delegation. Also claims `<KEY>-<number>` issue handles typed into global search
 * and resolves them through the typed `gritDb` business-key lookups. */
export class GritIssueHandler extends DG.DomainObjectHandler {
  constructor() { super(ISSUE_TYPE); }

  get name(): string { return 'Grit issue handler'; }

  isApplicable(x: any): boolean {
    return super.isApplicable(x) || this.handleString(x) != null;
  }

  getCaption(x: any): string {
    return this.handleString(x) ?? this.rowOrThrow(x).semValue;
  }

  /** Suggestion pattern for the global-search bar (DG.ObjectHandler.regexpExample).
   * Issue handles are `<PROJECTKEY>-<number>` (per-project numbering), so there is no
   * single fixed prefix; `_initGrit` registers a `<KEY>-\d+` detector per existing
   * project so `SemanticValue.parse('GRITEST-1')` tags it as a `grit.issue`. */
  get regexpExample(): {regexpMarkup: string, example: string, nonVariablePart: string} {
    return {regexpMarkup: '<project>-[0-9]+', example: 'GRITEST-1', nonVariablePart: '-'};
  }

  /** Raw handle string when [x] is a parsed SemanticValue (e.g. `GRIT-123` typed
   * into global search) rather than a resolved DomainRow. Claims ONLY the
   * `<KEY>-<number>` shape; the `grit.issue:<uuid>` colon form falls through to
   * the generic DomainHandleMeta (which navigates directly). */
  private handleString(x: any): string | null {
    return x instanceof DG.SemanticValue && typeof x.value === 'string' &&
      x.semType === ISSUE_TYPE && ISSUE_HANDLE_RE.test(x.value) ? x.value : null;
  }

  /** Resolves a `GRIT-123` handle (project key + per-project number) to its issue row
   * via the business-key lookups (`project.key`, then `issue.(project_id, number)`). */
  private async resolveHandle(handle: string): Promise<IssueRow | null> {
    const m = ISSUE_HANDLE_RE.exec(handle);
    if (m == null)
      return null;
    const project = await gritDb.projects.getByKey({key: m[1]});
    if (project == null)
      return null;
    return await gritDb.issues.getByKey({project_id: project.id, number: Number(m[2])});
  }

  /** Card for a handle typed into search: resolves asynchronously, shows the issue,
   * and opens its Entity View on click (the platform then renders it through this
   * same handler). */
  private renderHandleCard(handle: string): HTMLElement {
    const card = ui.divV([ui.divText(handle, {style: {fontWeight: 'bold'}})], 'd4-gallery-item');
    this.resolveHandle(handle).then((row) => {
      ui.empty(card);
      if (row == null) {
        card.appendChild(ui.divText(`${handle} — not found`));
        return;
      }
      card.appendChild(ui.divH([ui.divText(handle, {style: {fontWeight: 'bold'}}),
        ...this.badgeEls(row.status, row.priority)]));
      card.appendChild(ui.divText(row.title ?? ''));
      card.style.cursor = 'pointer';
      card.onclick = () => grok.shell.route(`/domains/grit/issue/${row.id}`);
    }).catch(() => {});
    return card;
  }

  private badgeEls(status?: string, priority?: string): HTMLElement[] {
    const res: HTMLElement[] = [];
    if (status != null)
      res.push(badge(status, statusColors[status] ?? '#8b95a1'));
    if (priority != null)
      res.push(badge(priority, priorityColors[priority] ?? '#8b95a1'));
    return res;
  }

  private badges(row: DG.DomainRow): HTMLElement[] {
    return this.badgeEls(row.values.status, row.values.priority);
  }

  renderIcon(x: any): HTMLElement { return ui.iconFA('bug'); }

  renderMarkup(x: any): HTMLElement {
    const h = this.handleString(x);
    if (h != null)
      return ui.span([this.renderIcon(x), ui.label(h)]);
    const row = this.rowOrThrow(x);
    return ui.span([this.renderIcon(row), ui.label(row.values.title ?? row.semValue), ...this.badges(row)]);
  }

  renderTooltip(x: any): HTMLElement {
    const h = this.handleString(x);
    if (h != null)
      return ui.divText(h);
    const row = this.rowOrThrow(x);
    return ui.divV([
      ui.divH([ui.label(row.values.title ?? row.semValue), ...this.badges(row)]),
      ui.divText(row.values.description ?? ''),
    ]);
  }

  renderCard(x: any): HTMLElement {
    const h = this.handleString(x);
    if (h != null)
      return this.renderHandleCard(h);
    const row = this.rowOrThrow(x);
    return ui.bind(x, ui.divV([
      ui.divH([
        ui.divText(row.values.title ?? row.semValue, {style: {fontWeight: 'bold'}}),
        ...this.badges(row),
      ]),
      ui.divText(row.values.description ?? ''),
    ], 'd4-gallery-item'));
  }

  renderGrid(grid: DG.Grid, options?: {items?: DG.DataFrame}): void {
    super.renderGrid(grid, options);
    const title = grid.col('title');
    if (title != null)
      title.width = 300;
  }
}
