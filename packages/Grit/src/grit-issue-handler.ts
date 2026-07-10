import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {gritDb, IssueRow} from './generated/db';

/** semType (= entity type) of `grit.issue` rows. */
const ISSUE_TYPE = 'grit.issue';

/** A Grit issue handle: `<PROJECTKEY>-<number>`, e.g. `GRITEST-1` (per-project numbering). */
const ISSUE_HANDLE_RE = /^([A-Za-z][A-Za-z0-9_]*)-(\d+)$/;

/** Badge presentation lives in css/grit.css (`.grit-badge`, `.grit-status-*`, `.grit-priority-*`). */
function badge(text: string, kind: 'status' | 'priority'): HTMLElement {
  return ui.span([text], `grit-badge grit-${kind}-${text.toLowerCase().replace(/\s+/g, '-')}`);
}

/** Renders a `grit.issue` domain row throughout the platform — Domain View cards,
 * the context panel, and the Entity View — overriding the generic Dart
 * `DomainRowMeta` for this table only (other grit/plates tables keep the defaults).
 * Registered at startup via {@link DG.ObjectHandler.register} (wins `forEntity`
 * dispatch) and discoverable by semType through the `gritIssueHandler` package
 * function (wins `forSemType('grit.issue')`). */
export class GritIssueHandler extends DG.ObjectHandler<DG.DomainRow> {
  get type(): string { return ISSUE_TYPE; }
  get name(): string { return 'Grit issue handler'; }

  isApplicable(x: any): boolean {
    if (x instanceof DG.DomainRow)
      return x.typeName === ISSUE_TYPE;
    if (x instanceof DG.SemanticValue) {
      if (x.value instanceof DG.DomainRow)
        return x.value.typeName === ISSUE_TYPE;
      // Claim ONLY the `<KEY>-<number>` handle shape; the `grit.issue:<uuid>` colon
      // form falls through to the generic DomainHandleMeta (which navigates directly).
      if (typeof x.value === 'string')
        return x.semType === ISSUE_TYPE && ISSUE_HANDLE_RE.test(x.value);
    }
    return false;
  }

  private issueOf(x: DG.DomainRow | DG.SemanticValue): DG.DomainRow {
    return x instanceof DG.SemanticValue ? x.value : x;
  }

  getCaption(x: DG.DomainRow): string {
    const h = this.handleString(x);
    return h != null ? h : this.issueOf(x).semValue;
  }

  /** Suggestion pattern for the global-search bar (DG.ObjectHandler.regexpExample).
   * Issue handles are `<PROJECTKEY>-<number>` (per-project numbering), so there is no
   * single fixed prefix; `_initGrit` registers a `<KEY>-\d+` detector per existing
   * project so `SemanticValue.parse('GRITEST-1')` tags it as a `grit.issue`. */
  get regexpExample(): {regexpMarkup: string, example: string, nonVariablePart: string} {
    return {regexpMarkup: '<project>-[0-9]+', example: 'GRITEST-1', nonVariablePart: '-'};
  }

  /** Raw handle string when [x] is a parsed SemanticValue (e.g. `GRIT-123` typed
   * into global search) rather than a resolved DomainRow. */
  private handleString(x: any): string | null {
    return x instanceof DG.SemanticValue && typeof x.value === 'string' ? x.value : null;
  }

  /** Resolves a `GRIT-123` handle (project key + per-project number) to its issue row. */
  private async resolveHandle(handle: string): Promise<IssueRow | null> {
    const m = ISSUE_HANDLE_RE.exec(handle);
    if (m == null)
      return null;
    const projects = await gritDb.project.query({filter: `key = "${m[1]}"`, limit: 1});
    if (projects.length === 0)
      return null;
    const issues = await gritDb.issue.query(
      {filter: `project_id = "${projects[0].id}" and number = ${m[2]}`, limit: 1});
    return issues.length === 1 ? issues[0] : null;
  }

  /** Card for a handle typed into search: resolves asynchronously, shows the issue,
   * and opens its Entity View on click (the platform then renders it through this
   * same handler's renderView). */
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
      res.push(badge(status, 'status'));
    if (priority != null)
      res.push(badge(priority, 'priority'));
    return res;
  }

  private badges(row: DG.DomainRow): HTMLElement[] {
    return this.badgeEls(row.values.status, row.values.priority);
  }

  renderIcon(x: DG.DomainRow): HTMLElement { return ui.iconFA('bug'); }

  renderMarkup(x: DG.DomainRow): HTMLElement {
    const h = this.handleString(x);
    if (h != null)
      return ui.span([this.renderIcon(x), ui.label(h)]);
    const row = this.issueOf(x);
    return ui.span([this.renderIcon(row), ui.label(row.values.title ?? row.semValue), ...this.badges(row)]);
  }

  renderTooltip(x: DG.DomainRow): HTMLElement {
    const h = this.handleString(x);
    if (h != null)
      return ui.divText(h);
    const row = this.issueOf(x);
    return ui.divV([
      ui.divH([ui.label(row.values.title ?? row.semValue), ...this.badges(row)]),
      ui.divText(row.values.description ?? ''),
    ]);
  }

  renderCard(x: DG.DomainRow): HTMLElement {
    const h = this.handleString(x);
    if (h != null)
      return this.renderHandleCard(h);
    const row = this.issueOf(x);
    return ui.bind(x, ui.divV([
      ui.divH([
        ui.divText(row.semValue, {style: {fontWeight: 'bold'}}),
        ...this.badges(row),
      ]),
      ui.divText(row.values.title ?? ''),
    ], 'd4-gallery-item'));
  }

  renderView(x: DG.DomainRow): HTMLElement {
    const row = this.issueOf(x);
    const v = row.values;
    const timeline = ui.divV([ui.loader()]);
    this.buildTimeline(row).then((el) => ui.empty(timeline).appendChild(el)).catch(() => {});
    return ui.divV([
      ui.divH([ui.h1(row.semValue), ...this.badges(row)]),
      ui.h2(v.title ?? ''),
      ui.divText(v.description ?? ''),
      ui.h3('Timeline'),
      timeline,
    ], 'grit-issue-view');
  }

  renderProperties(x: DG.DomainRow): HTMLElement {
    const row = this.issueOf(x);
    const v = row.values;
    const timeline = ui.divV([ui.loader()]);
    this.buildTimeline(row).then((el) => ui.empty(timeline).appendChild(el)).catch(() => {});
    return ui.divV([
      ui.divH([ui.h2(row.semValue), ...this.badges(row)]),
      ui.tableFromMap({
        'Title': v.title ?? '',
        'Description': v.description ?? '',
        'Status': v.status ?? '',
        'Priority': v.priority ?? '',
        'Project': v.project ?? '',
      }),
      ui.h3('Timeline'),
      timeline,
    ], 'grit-issue-properties');
  }

  private async buildTimeline(row: DG.DomainRow): Promise<HTMLElement> {
    const users = new Map<string, string>();
    for (const u of await grok.dapi.users.list())
      users.set(u.id, u.friendlyName);
    const audit = await gritDb.issue.audit(row.id);
    if (audit.length === 0)
      return ui.divText('No changes recorded.');
    return ui.divV(audit.map((a: any) => {
      const who = users.get(a.actor_id) ?? 'unknown';
      const when = a.ts?.substring(0, 16).replace('T', ' ') ?? '';
      const changes = a.op === 'update' && a.after != null ?
        ': ' + Object.keys(a.after).map((k) => `${k} → ${a.after[k]}`).join(', ') : '';
      return ui.divText(`${when}  ${who} — ${a.op}${changes}`);
    }));
  }
}
