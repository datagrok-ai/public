/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import {computed, Control, Section, SharedSession, Splitter} from '@datagrok-libraries/u2';
import {appView, domains} from '@datagrok-libraries/u2/src/dg/index.js';
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
  // live: the list follows the server, so another session's insert shows up here within a poll
  return (await domains.table('stockroom.substance')).app({live: true});
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
  // live: the list follows the server, so a container another session moves in shows up here
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
  const section = (title: string, content: Control): Section => new Section({title, collapsible: false})
    .add(content);
  const panes = new Splitter([section('Locations', tree), section('Containers', table)],
    {direction: 'horizontal', sizes: [0.35, 0.65]});
  return appView({
    name: 'Locations',
    content: panes,
    own: [source],
    status: where,
    path: '/apps/Stockroom/Locations',
  });
}
