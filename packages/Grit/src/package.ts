/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Splitter} from '@datagrok-libraries/u2';
import {appView, domains, domainForm, domainList, newButton, saveButton, discardButton}
  from '@datagrok-libraries/u2/src/dg/index.js';
import {GritIssueHandler} from './grit-issue-handler';
import {gritDb} from './generated/db';
export * from './package.g';

// the skins of everything the two u2 views render (tokens first: every other sheet reads them)
import '@datagrok-libraries/u2/css/tokens.css';
import '@datagrok-libraries/u2/css/elements.css';
import '@datagrok-libraries/u2/css/buttons.css';
import '@datagrok-libraries/u2/css/inputs.css';
import '@datagrok-libraries/u2/css/number.css';
import '@datagrok-libraries/u2/css/date.css';
import '@datagrok-libraries/u2/css/choice.css';
import '@datagrok-libraries/u2/css/combobox.css';
import '@datagrok-libraries/u2/css/tags.css';
import '@datagrok-libraries/u2/css/typeahead.css';
import '@datagrok-libraries/u2/css/entity.css';
import '@datagrok-libraries/u2/css/form.css';
import '@datagrok-libraries/u2/css/list.css';
import '@datagrok-libraries/u2/css/menu.css';
import '@datagrok-libraries/u2/css/splitter.css';
import '@datagrok-libraries/u2/css/async.css';
import '@datagrok-libraries/u2/css/dialog.css';
import '@datagrok-libraries/u2/css/notify.css';
import '@datagrok-libraries/u2/css/tooltip.css';
import '@datagrok-libraries/u2/css/badge.css';
import '@datagrok-libraries/u2/css/domain.css';

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

//name: Issues
//description: Issue tracker over entity-mapped domain schemas — the platform's Domain View with the filter panel
//tags: app
//meta.icon: images/bug.svg
//input: string path {meta.url: true; optional: true}
//output: view result
export async function issuesApp(path?: string): Promise<DG.ViewBase> {
  const route = (path ?? '').split('?')[0].replace(/^\/+|\/+$/g, '').toLowerCase();
  if (route === 'create-issue')
    return await createIssueView();
  if (route === 'projects')
    return await projectsView();
  const view = DG.DomainView.create({schema: 'grit', table: 'issue'});
  whenDocked(view, () => view.showFilters());
  return view;
}

/** Runs [action] once [view] is in the DOM: `DomainView.showFilters` and the ribbon
 * name slot need the view docked, and an app function returns its view BEFORE the
 * shell docks it. */
function whenDocked(view: DG.ViewBase, action: () => void, tries: number = 100): void {
  view.root.isConnected ? action() : tries > 0 ? void setTimeout(() => whenDocked(view, action, tries - 1), 100) : null;
}

//name: issuesTreeBrowser
//input: dynamic treeNode
//meta.role: appTreeBrowser
//meta.app: Issues
export function issuesTreeBrowser(treeNode: DG.TreeViewGroup): void {
  issuesNode = treeNode;
  treeNode.item('Create Issue').onSelected.subscribe(async () =>
    grok.shell.preview = await createIssueView());
  treeNode.item('Projects').onSelected.subscribe(async () =>
    grok.shell.preview = await projectsView());
}

/** The 'Issues' app group node, once the browse tree has expanded it. */
let issuesNode: DG.TreeViewGroup | null = null;

async function createIssueView(): Promise<DG.ViewBase> {
  const [issues, me] = await Promise.all([domains.table('grit.issue'), grok.dapi.users.current()]);
  const src = issues.draft({reporter: me.id});
  return childView(appView({name: 'Create Issue', content: domainForm(src), own: [src],
    ribbon: [[newButton(src, (last) => ({reporter: me.id, project_id: last?.project_id})),
      saveButton(src), discardButton(src)]],
    status: src.summary}), 'create-issue');
}

async function projectsView(): Promise<DG.ViewBase> {
  const projects = await domains.table('grit.project');
  const src = projects.source();
  const content = new Splitter([domainList(src, {mode: 'cards'}), domainForm(src)],
    {direction: 'horizontal', sizes: [40, 60]});
  return childView(appView({name: 'Projects', content, own: [src],
    ribbon: [[newButton(src), saveButton(src), discardButton(src)]], status: src.summary}), 'projects');
}

/** Decorates a child view of the Issues app: the `/apps/Grit/<route>` address (the
 * platform strips the duplicated app prefix when it is already implied) and
 * `Issues / <name>` breadcrumbs; clicking 'Issues' goes back to the main view. */
function childView(view: DG.ViewBase, route: string): DG.ViewBase {
  view.basePath = `/apps/Grit/${route}`;
  whenDocked(view, () => {
    const crumbs = ui.breadcrumbs(['Issues', view.name]);
    crumbs.onPathClick.subscribe((path) => {
      if (path[path.length - 1] === 'Issues')
        issuesNode != null ? issuesNode.currentItem = issuesNode : issuesApp().then((v) => grok.shell.addPreview(v));
    });
    const nameRoot = view.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
    if (nameRoot != null) {
      nameRoot.textContent = '';
      nameRoot.appendChild(crumbs.root);
    }
  });
  return view;
}
