/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {button, computed, Control, Section, Splitter, SuggestInput} from '@datagrok-libraries/u2';
import type {RowView, Text} from '@datagrok-libraries/u2';
import {SharedSession} from '@datagrok-libraries/ems';
import {appView, Editors} from '@datagrok-libraries/u2/src/dg/index.js';
import {domains, DomainApp} from '@datagrok-libraries/ems/src/dg';
import type {DomainTable} from '@datagrok-libraries/ems/src/dg';
import {stockroomDb} from './generated/db';
export * from './package.g';

// the skins of everything the u2 domain stack renders
import '@datagrok-libraries/ems/src/dg/styles.js';
import '../css/stockroom.css';

export const _package = new DG.Package();

// a site is typed freely, offered from the sites the locations already name, so a container's
// location picker is not gated on a name typed from memory
Editors.register({
  match: (prop) => prop.semType === 'Stockroom site',
  create: (_, options) => new SuggestInput({...options, minChars: 0, openOnFocus: true, source: async (text) => {
    const sites = await stockroomDb.locations.aggregate({groupBy: ['site'], measures: [{fn: 'count'}], sort: 'site'});
    return sites.map((s) => String(s.site)).filter((s) => s.toLowerCase().includes(text.toLowerCase()));
  }}),
});

/** A view per stockroom table under the Stockroom app, the substances at its root. */
const views = domains.views({
  path: '/apps/Stockroom/Stockroom',
  views: [
    {table: 'stockroom.substance', name: 'Substances', segment: ''},
    {table: 'stockroom.container', name: 'Containers', options: {sort: 'label'}},
    {table: 'stockroom.vendor', name: 'Vendors'},
    {table: 'stockroom.purchase_order', name: 'Purchase Orders'},
    {table: 'stockroom.sds_document', name: 'Safety Data Sheets'},
    {table: 'stockroom.disposal', name: 'Disposals'},
    {table: 'stockroom.hazard_class', name: 'Hazard Classes'},
    {table: 'stockroom.h_statement', name: 'Hazard Statements'},
    {table: 'stockroom.p_statement', name: 'Precautionary Statements'},
  ],
});

//name: Stockroom
//description: Chemical stockroom on the GHS classification — the zero-code app over databases/stockroom/schema.json
//tags: app
//meta.icon: images/flask.svg
//input: string path {meta.url: true; optional: true}
//output: view result
export async function stockroomApp(path?: string): Promise<DG.ViewBase> {
  return views.open(path ?? '');
}

//name: stockroomTreeBrowser
//input: dynamic treeNode
//meta.role: appTreeBrowser
//meta.app: Stockroom
export function stockroomTreeBrowser(treeNode: DG.TreeViewGroup): void {
  views.tree(treeNode);
}

/** The Containers app over everything stored under `location`, opened straight on a draft that
 * starts there — the Locations view's New container, and what its test drives. The view is named
 * after the node, because the query it carries reads as a uuid. */
export function newContainerAt(containers: DomainTable, location: string,
  options: {site?: string, caption?: string} = {}): DG.ViewBase {
  const view = containers.app({
    name: options.caption === undefined ? 'Containers' : `Containers in ${options.caption}`,
    path: '/apps/Stockroom/Containers', query: `location_id under "${location}"`, sort: 'label',
    // `site` is required on a container and is what narrows its location picker
    // (`location_id.filter: site = $site`); the picked node carries it, so the draft starts with both
    defaults: {location_id: location, ...(options.site === undefined ? {} : {site: options.site})}});
  grok.shell.addView(view);
  void DomainApp.of(view)!.create();
  return view;
}

//name: Locations
//description: The stockroom location tree, and the containers stored anywhere under the selected node
//tags: app
//meta.icon: images/flask.svg
//output: view result
export async function stockroomLocations(): Promise<DG.ViewBase> {
  const [locations, containers] = await Promise.all([
    domains.table('stockroom.location'), domains.table('stockroom.container')]);
  const session = new SharedSession();
  // live: a bare source is not polled by default (R-b flips the app's option, not this one), so a
  // container another session moves in shows up here only because this asks for it
  const source = SharedSession.runWith(session, () =>
    containers.source({session, live: true, pageSize: 100, sort: 'label'}));
  const addAt = (node: RowView) => newContainerAt(containers, node.id,
    {site: typeof node.site === 'string' ? node.site : undefined, caption: locations.renderer.caption(node)});
  // the "all" row is the way back to every container once a node has been picked
  const tree = domains.tree(locations, {allNode: 'All locations',
    actions: [{name: 'New container', icon: 'plus', run: addAt}]});
  const table = domains.dataTable(source,
    {columns: ['label', 'substance_id', 'location_id', 'lot', 'quantity', 'unit', 'expires'],
      onActivate: (row) => containers.open(row)});
  // short values get short tracks, so the names keep their room
  const widths: Record<string, string> = {label: '56px', quantity: '72px', unit: '48px'};
  table.effect(() => {
    const grid = table.table.value;
    if (grid !== null)
      grid.columns = grid.columns.map((c) => widths[c.name] === undefined ? c : {...c, width: widths[c.name]});
  });
  // one node, one subtree: `under` matches the node and everything below it, through the ref
  // column into the tree — nothing selected is the whole table
  tree.effect(() => {
    const node = tree.selected.value;
    source.query.value = node === null ? '' : `location_id under "${node.id}"`;
  });
  // the two panes say what they are, and the containers' header says what the tree did to the collection
  const where = computed(() => {
    const node = tree.selected.value;
    return node === null ? source.summary.value :
      `${source.summary.value} in ${locations.renderer.caption(node)}, including sublocations`;
  });
  const addButton = button('New container', () => {
    const node = tree.selected.value;
    if (node !== null)
      addAt(node);
  });
  const add = new Control(addButton);
  add.root.dataset.u2 = 'new-container-button';
  add.effect(() => {
    const none = tree.selected.value === null;
    addButton.disabled = none;
    // a disabled button takes no pointer events, so the hover that shows the hint lands on the host
    add.root.title = addButton.title = none ? 'Select a location first' : '';
  });
  const section = (title: string, content: Control, hint?: Text): Section =>
    new Section({title, collapsible: false, hint}).add(content);
  const panes = new Splitter([section('Locations', tree), section('Containers', table, where)],
    {direction: 'horizontal', sizes: [0.22, 0.78]});
  panes.root.classList.add('stockroom-locations');
  return appView({
    name: 'Locations',
    content: panes,
    own: [source],
    ribbon: [[add]],
    path: '/apps/Stockroom/Locations',
  });
}
