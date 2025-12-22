// viewer descriptors let you inspect widget properties without creating an instance of the widget
// here, we construct a tree with all viewers and all registered properties.

const tree = ui.tree();

for (const d of DG.WidgetDescriptor.getDescriptors()) {
  const viewerNode = tree.group(d.name);
  viewerNode.icon = d.createIcon();
  viewerNode.expanded = false;
  for (const prop of d.properties)
    viewerNode.item(prop.name, prop);
}

tree.root