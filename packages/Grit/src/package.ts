/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {badge, div, divV, span} from '@datagrok-libraries/u2';
import type {BadgeVariant, RowView} from '@datagrok-libraries/u2';
import {DomainApp} from '@datagrok-libraries/u2/src/dg/index.js';
import type {DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';
import {GritIssueHandler, lookupName, warmLookups} from './grit-issue-handler';
import {gritDb, IssueRow} from './generated/db';
import {getGritDb, GritDb} from './generated/db-ui';
export * from './package.g';

// the skins of everything the u2 domain stack renders
import '@datagrok-libraries/u2/src/dg/domain/styles.js';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: autostart
//meta.autostartImmediate: true
export function _initGrit(): void {
  DG.ObjectHandler.register(new GritIssueHandler());
  registerIssueHandleDetectors();
}

/** Registers a `<KEY>-\d+` detector per existing Grit project so that a real issue
 * handle (e.g. `GRITEST-1`) typed into global search is tagged as a `grit.issue`
 * and resolved by the handler (WO-31, ARCHITECTURE §8.6 discovery). Issue numbers
 * are per-project, so the fixed part is the project's own key, not a literal. */
async function registerIssueHandleDetectors(): Promise<void> {
  try {
    const projects = await gritDb.projects.query({});
    for (const p of projects) {
      const key = (p.key ?? '').replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
      if (key.length > 0)
        DG.SemanticValue.registerRegExpDetector('grit.issue', `${key}-\\d+`, `Grit issue (${p.key})`);
    }
  } catch (e) {
    console.error('Grit issue-handle detectors not registered:', e);
  }
}

//name: gritIssueHandler
//description: Custom rendering for grit.issue domain rows (status/priority badges)
//tags: objectHandler, objectHandler-grit.issue
//output: object handler
export function gritIssueHandler(): DG.ObjectHandler<DG.DomainRow> {
  return new GritIssueHandler();
}

const priorityVariants: {[name: string]: BadgeVariant} = {critical: 'error', high: 'warning', medium: 'accent'};

/** The id of the `closed` status; null where the lookup table has no such row, so nothing closes. */
let closedId: string | null = null;

let _grit: Promise<GritDb> | undefined;

/** The typed schema handles (one await, cached per page) with what Grit declares on `issue`:
 * the two actions, the closing rule and the card — made once, shared by every view and test. */
export function openGrit(): Promise<GritDb> {
  return _grit ??= Promise.all([getGritDb(), gritDb.statuses.getByKey({name: 'closed'}), warmLookups()])
    .then(([db, closed]) => {
      closedId = closed?.id ?? null;
      declareIssues(db.tables.issues);
      return db;
    });
}

function declareIssues(issues: DomainTable<IssueRow>): void {
  const me = () => grok.shell.user.id;
  issues.actions.add({name: 'Assign to me', icon: 'user', requires: 'edit',
    when: (r) => r.assignee !== me(), run: (r) => { r.assignee = me(); }});
  issues.actions.add({name: 'Close', icon: 'check', requires: 'edit',
    when: (r) => closedId !== null && r.status_id !== closedId, run: (r) => { r.status_id = closedId!; }});
  issues.validators.add('status_id', (v, r) =>
    closedId !== null && v === closedId && !r.assignee ? 'Assign before closing' : null);
  issues.renderer = {...issues.renderer, card: (r) => divV([
    div([span(`#${r.number} ${r.title}`, 'u2-domain-card-title'), ...priorityBadge(r)]),
    ...(r.description ? [span(r.description, 'u2-domain-card-description')] : []),
  ], 'u2-domain-card')};
}

function priorityBadge(r: RowView<IssueRow>): HTMLElement[] {
  const name = lookupName(r.priority_id);
  return name === null ? [] : [badge(name, {variant: priorityVariants[name] ?? 'default'})];
}

/** The Issues app: the platform ribbon plus the Mine / Open presets, and a key for each action. */
export class IssuesApp extends DomainApp {
  shortcuts = {'m': 'Assign to me', 'c': 'Close'};
  private _full: ReturnType<DomainApp['ribbon']> | undefined;

  ribbon(): ReturnType<DomainApp['ribbon']> {
    return this._full ??= [...super.ribbon(),
      [this.presets(['Mine', 'assignee = $me'], ['Open', 'status_id.name != "closed"'])]];
  }
}

//name: Issues
//description: Issue tracker over entity-mapped domain schemas — the u2 app over grit.issue with Grit's actions, presets and shortcuts
//tags: app
//meta.icon: images/bug.svg
//input: string path {meta.url: true; optional: true}
//output: view result
export async function issuesApp(path?: string): Promise<DG.ViewBase> {
  const {issues} = (await openGrit()).tables;
  // live: the list follows the server, so an issue another session files shows up here
  const view = issues.app({name: 'Issues', path: '/apps/Grit/Issues', app: IssuesApp, mode: 'cards',
    children: {tables: ['comment']}, live: true});
  void DomainApp.of(view)!.open(path || undefined);
  return view;
}
