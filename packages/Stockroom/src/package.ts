/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {button, computed, Control, Section, SharedSession, Splitter} from '@datagrok-libraries/u2';
import {appView, domains, DomainApp} from '@datagrok-libraries/u2/src/dg/index.js';
import type {DomainTable} from '@datagrok-libraries/u2/src/dg/index.js';
export * from './package.g';

// the skins of everything the u2 domain stack renders
import '@datagrok-libraries/u2/src/dg/domain/styles.js';

export const _package = new DG.Package();

//name: Stockroom
//description: Chemical stockroom on the GHS classification — the zero-code app over databases/stockroom/schema.json
//tags: app
//meta.icon: images/flask.svg
//output: view result
export async function stockroomApp(): Promise<DG.ViewBase> {
  return (await domains.table('stockroom.substance')).app();
}

/** The Containers app over everything stored under `location`, opened straight on a draft that
 * starts there — the Locations view's New container, and what its test drives. The view is named
 * after the node, because the query it carries reads as a uuid. */
export function newContainerAt(containers: DomainTable, location: string,
  options: {site?: string, caption?: string} = {}): DG.ViewBase {
  const view = containers.app({
    name: options.caption === undefined ? 'Containers' : `Containers in ${options.caption}`,
    path: '/apps/Stockroom/Containers', query: `location_id under "${location}"`,
    // `site` is required on a container and is what narrows its location picker
    // (`location_id.filter: site = $site`); the picked node carries it, so the draft starts with both
    defaults: {location_id: location, ...(options.site === undefined ? {} : {site: options.site})}});
  grok.shell.addView(view);
  void DomainApp.of(view)!.goTo('entity', DomainApp.NEW);
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
    containers.source({session, live: true, pageSize: 100}));
  // the "all" row is the way back to every container once a node has been picked
  const tree = domains.tree(locations, {allNode: 'All locations'});
  const table = domains.dataTable(source,
    {columns: ['label', 'substance_id', 'location_id', 'lot', 'quantity', 'unit', 'expires'],
      onActivate: (row) => containers.open(row)});
  // one node, one subtree: `under` matches the node and everything below it, through the ref
  // column into the tree — nothing selected is the whole table
  tree.effect(() => {
    const node = tree.selected.value;
    source.query.value = node === null ? '' : `location_id under "${node.id}"`;
  });
  // the two panes say what they are, and the status says what the tree did to the collection
  const where = computed(() => {
    const node = tree.selected.value;
    return node === null ? source.summary.value :
      `${source.summary.value} in ${locations.renderer.caption(node)}, including sublocations`;
  });
  const addButton = button('New container', () => {
    const node = tree.selected.value;
    if (node !== null) {
      newContainerAt(containers, node.id, {site: typeof node.site === 'string' ? node.site : undefined,
        caption: locations.renderer.caption(node)});
    }
  });
  const add = new Control(addButton);
  add.root.dataset.u2 = 'new-container-button';
  add.effect(() => {
    const none = tree.selected.value === null;
    addButton.disabled = none;
    // a disabled button takes no pointer events, so the hover that shows the hint lands on the host
    add.root.title = addButton.title = none ? 'Select a location first' : '';
  });
  const section = (title: string, content: Control): Section => new Section({title, collapsible: false})
    .add(content);
  const panes = new Splitter([section('Locations', tree), section('Containers', table)],
    {direction: 'horizontal', sizes: [0.35, 0.65]});
  return appView({
    name: 'Locations',
    content: panes,
    own: [source],
    ribbon: [[add]],
    status: where,
    path: '/apps/Stockroom/Locations',
  });
}
