/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {badge, div, divV, span} from '@datagrok-libraries/u2';
import type {BadgeVariant, RowView} from '@datagrok-libraries/u2';
import {DomainApp, domains} from '@datagrok-libraries/ems/src/dg';
import type {DomainRibbon} from '@datagrok-libraries/ems/src/dg';
import {GritIssueHandler, lookupName, warmLookups} from './grit-issue-handler';
import {gritDb, IssueRow} from './generated/db';
import {getGritDb} from './generated/db-ui';
export * from './package.g';

// the skins of everything the u2 domain stack renders
import '@datagrok-libraries/ems/src/dg/styles.js';

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

let _closed: Promise<string | null> | undefined;

/** The id of the `closed` status, looked up once per page; null where the lookup table has no such row, so nothing
 * closes. */
function closedStatus(): Promise<string | null> {
  if (_closed === undefined) {
    _closed = gritDb.statuses.getByKey({name: 'closed'}).then((closed) => closed?.id ?? null);
    _closed.catch(() => _closed = undefined);
  }
  return _closed;
}

/** What Grit declares on `issue` — the two actions, the closing rule and the card — on every handle over it: the
 * Issues app's, a child tab's, the `/domains` route's. */
domains.configure<IssueRow>('grit.issue', async (issues) => {
  const [closedId] = await Promise.all([closedStatus(), warmLookups()]);
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
}, {id: 'grit'});

function priorityBadge(r: RowView<IssueRow>): HTMLElement[] {
  const name = lookupName(r.priority_id);
  return name === null ? [] : [badge(name, {variant: priorityVariants[name] ?? 'default'})];
}

/** The Issues app: the platform ribbon plus the Mine / Open presets, and a key for each action. */
export class IssuesApp extends DomainApp {
  shortcuts = {'m': 'Assign to me', 'c': 'Close'};
  private _full: DomainRibbon | undefined;

  ribbon(): DomainRibbon {
    return this._full ??= {...super.ribbon(),
      presets: [this.presets(['Mine', 'assignee = $me'], ['Open', 'status_id.name != "closed"'])]};
  }
}

/** A view per Grit table — the issues where they were, at `/apps/Grit/Issues` — over the handles {@link getGritDb}
 * opens, so the actions, the closing rule and the card hold in each. Grit's one app answers everything under
 * `/apps/Grit`. */
const views = domains.views({
  path: '/apps/Grit',
  handle: async (address) => Object.values((await getGritDb()).tables).find((t) => t.address === address)!,
  views: [
    {table: 'grit.issue', name: 'Issues', options: {app: IssuesApp, mode: 'cards',
      children: {tables: ['comment']}}},
    {table: 'grit.project', name: 'Projects'},
    {table: 'grit.label', name: 'Labels'},
    {table: 'grit.status', name: 'Statuses'},
    {table: 'grit.priority', name: 'Priorities'},
    {table: 'grit.issue_type', name: 'Issue Types'},
  ],
});

//name: Issues
//description: Issue tracker over entity-mapped domain schemas — the u2 app over grit.issue with Grit's actions, presets and shortcuts
//tags: app
//meta.icon: images/bug.svg
//input: string path {meta.url: true; optional: true}
//output: view result
export async function issuesApp(path?: string): Promise<DG.ViewBase> {
  return views.open(path ?? '');
}

//name: issuesTreeBrowser
//input: dynamic treeNode
//meta.role: appTreeBrowser
//meta.app: Issues
export function issuesTreeBrowser(treeNode: DG.TreeViewGroup): void {
  views.tree(treeNode);
}
