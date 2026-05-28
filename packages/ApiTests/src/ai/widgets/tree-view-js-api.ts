// DG.TreeViewNode / DG.TreeViewGroup — core/client/d4/lib/src/widgets/tree_view/tree_view.dart
// (scenario: tree-view-js-api)
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectBoolGetSet, expectFiresWithin, subscribeAll, expectNoThrow, wait} from '../helpers';

category('AI: Widgets: TreeView JS API', () => {
  test('build structure: children / items / item / group / addNode / getOrCreateGroup / rootNode / ' +
    'parent', async () => {
    const root = ui.tree();
    expect(root instanceof DG.TreeViewGroup, true);
    const g1 = root.group('Group 1');
    expect(g1 instanceof DG.TreeViewGroup, true);
    const leaf = g1.item('Leaf', {id: 1});
    expect(leaf instanceof DG.TreeViewNode, true);
    expect(leaf.text, 'Leaf');
    // captionLabel + text round-trip
    expect(leaf.captionLabel instanceof HTMLElement, true);
    leaf.text = 'Leaf-renamed';
    expect(leaf.text, 'Leaf-renamed');
    // rootNode / parent navigation
    expect(leaf.parent instanceof DG.TreeViewGroup, true);
    expect(leaf.rootNode instanceof DG.TreeViewGroup, true);
    expect(leaf.rootNode.dart, root.dart);
    // children / items
    expect(root.children.length, 1);
    expect(g1.items.length, 1);
    // getOrCreateGroup is idempotent for the same text
    const a = root.getOrCreateGroup('Shared');
    const b = root.getOrCreateGroup('Shared');
    expect(a.dart, b.dart);
    expect(root.children.length, 2);
    // addNode at an explicit index reorders, remove prunes
    const moved = root.group('Moved');
    moved.remove();
    expect(root.children.length, 2);
    root.addNode(moved, 0);
    expect(root.children.length, 3);
    expect(root.children[0].text, 'Moved');
    root.root.remove();
  });

  test('addItems bulk + boundary empty tree', async () => {
    const root = ui.tree();
    expect(root.items.length, 0);
    // boundary: empty bulk add is a no-op
    expectNoThrow(() => root.addItems([]));
    expect(root.items.length, 0);
    root.addItems(['alpha', 'beta', 'gamma']);
    expect(root.items.length, 3);
    expect(root.items[0].text, 'alpha');
    expect(root.items[2].text, 'gamma');
    root.root.remove();
  });

  test('checkbox state: enableCheckBox / checkBox / checked + onNodeCheckBoxToggled fires on parent', async () => {
    const root = ui.tree();
    const item = root.item('checkable');
    item.enableCheckBox(false);
    expect(item.checkBox instanceof HTMLInputElement, true);
    expect(item.checked, false);
    // NODE_CHECKBOX_TOGGLED fires on the parent group (rootNode here), via the checked= setter.
    await expectFiresWithin(root.onNodeCheckBoxToggled, () => {item.checked = true;});
    expect(item.checked, true);
    root.root.remove();
  });

  test('autoCheckChildren get/set round-trip', async () => {
    const root = ui.tree();
    expectBoolGetSet(root, 'autoCheckChildren');
    root.root.remove();
  });

  test('currentItem get/set + onSelectedNodeChanged (root) + onSelected (item) fire on selection', async () => {
    const root = ui.tree();
    const g = root.group('Sel group');
    const item = g.item('selectable');
    // onSelectedNodeChanged fires on the root; trigger via currentItem= on root.
    await expectFiresWithin(root.onSelectedNodeChanged, () => {root.currentItem = item;});
    expect(root.currentItem.dart, item.dart);
    // onSelected fires on the item itself when it becomes selected again.
    const item2 = g.item('selectable-2');
    await expectFiresWithin(item2.onSelected, () => {root.currentItem = item2;});
    expect(root.currentItem.dart, item2.dart);
    root.root.remove();
  });

  test('icon + tag round-trip on a node', async () => {
    const root = ui.tree();
    const item = root.item('tagged');
    const iconEl = ui.iconFA('star');
    item.icon = iconEl;
    expect(item.icon != null, true);
    item.tag = {kind: 'demo', n: 7};
    expect(item.tag.kind, 'demo');
    expect(item.tag.n, 7);
    root.root.remove();
  });

  test('removeChildrenWhere prunes recursively; boundary no-op on empty tree', async () => {
    // boundary: empty tree
    const empty = ui.tree();
    expectNoThrow(() => empty.removeChildrenWhere(() => true));
    expect(empty.children.length, 0);
    empty.root.remove();

    const root = ui.tree();
    const keep = root.group('keep');
    keep.item('keep-leaf');
    const drop = root.group('drop-me');
    drop.item('drop-leaf');
    const nested = keep.group('nested');
    nested.item('drop-me');
    expect(root.children.length, 2);
    // Prune every node whose text contains 'drop' — recurses into groups.
    root.removeChildrenWhere((n) => n.text.indexOf('drop') >= 0);
    expect(root.children.length, 1);
    expect(root.children[0].text, 'keep');
    expect(nested.items.length, 0);
    root.root.remove();
  });

  test('onNodeExpanding / onChildNodeExpanding / onChildNodeExpandedChanged fire on first expand; ' +
    'DOM-only events subscribe cleanly', async () => {
    const root = ui.tree();
    // Create the group COLLAPSED before subscribing so the first expand fires the events.
    const g = root.group('expandable', null, false);
    expect(g.expanded, false);
    // onNodeExpanding fires on the group's own bus; the child-* events fire on the root.
    await expectFiresWithin(g.onNodeExpanding, () => {g.expanded = true;});
    expect(g.expanded, true);

    const g2 = root.group('expandable-2', null, false);
    await expectFiresWithin(root.onChildNodeExpanding, () => {g2.expanded = true;});

    const g3 = root.group('expandable-3', null, false);
    await expectFiresWithin(root.onChildNodeExpandedChanged, () => {g3.expanded = true;});

    // Remaining DOM-only events: smoke-assert they are healthy rxjs Observables.
    subscribeAll([
      root.onNodeAdded, root.onNodeContextMenu, root.onNodeMouseEnter,
      root.onNodeMouseLeave, root.onNodeEnter,
    ])();

    // loadSources is referenced (symbol touched) but not invoked — it needs an HttpDataSource.
    expect(typeof root.loadSources, 'function');
    await wait(10);
    root.root.remove();
  });
}, {owner: 'agolovko@datagrok.ai'});
