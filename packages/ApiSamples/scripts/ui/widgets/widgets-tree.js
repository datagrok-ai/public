// Builds a visual tree of currently active widgets,
// along with the functions and properties that these widgets expose

function addWidgetToTree(widget, parentGroup) {
  const childCount = widget.children?.length ?? 0;
  const label = widget.type + (childCount === 0 ? '' : '(' + childCount + ')');
  const group = parentGroup.group(label, widget, false);

  // Only widgets that opt in expose these (DomainGrid / DomainForm do); the base
  // Widget has neither, so a plain widget contributes just its label.
  const functions = typeof widget.getFunctions === 'function' ? widget.getFunctions() : null;
  if (functions && functions.length > 0) {
    const funcNode = group.group('Functions (' + functions.length + ')', null, false);
    for (const f of functions)
      funcNode.item(f.name);
  }

  const properties = typeof widget.getProperties === 'function' ? widget.getProperties() : null;
  if (properties && properties.length > 0) {
    const propNode = group.group('Properties(' + properties.length + ')', null, false);
    for (const p of properties)
      propNode.item(p.name);
  }

  if (widget.children) {
    for (const child of widget.children)
      addWidgetToTree(child, group);
  }
}

const tree = DG.TreeViewGroup.tree();
for (const widget of DG.Widget.getAll()) {
  if (widget.parent == null)
    addWidgetToTree(widget, tree);
}

tree.root;