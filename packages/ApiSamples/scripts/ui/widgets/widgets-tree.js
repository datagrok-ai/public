// Builds a visual tree of currently active widgets,
// along with the functions and properties that these widgets expose

function addWidgetToTree(widget, parentGroup) {
  const label = widget.type + (widget.children.length === 0 ? '' : '(' + widget.children.length + ')');
  const group = parentGroup.group(label, widget, false);

  const functions = widget.getFunctions();
  if (functions && functions.length > 0) {
    const funcNode = group.group('Functions (' + functions.length + ')', null, false);
    for (const f of functions)
      funcNode.item(f.name);
  }

  const properties = widget.getProperties();
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