// Demonstrates drag-and-drop APIs: ui.makeDraggable and ui.makeDroppable.

let simple = ui.input.search('', {value: 'Drop a column or the Heinlein image here...'});

// Simple drop zone: accept anything, inspect via args.dragObject.
ui.makeDroppable(simple.input, {
  acceptDrop: (_) => true,
  doDrop: (args) => {
    const o = args.dragObject;
    if (o instanceof DG.Column)
      simple.value = `${o.name} dropped (copying=${args.copying}, link=${args.link})`;
    else if (o instanceof HTMLDivElement)
      simple.value = o.style.backgroundImage;
  }
});

let advanced = ui.divText('Drop a column here', {style: {
  padding: '20px', border: '2px dashed var(--grey-3)', borderRadius: '4px'
}});

// Rich drop zone: veto via acceptDrag, react to begin/end/hover, customize suggestion text.
ui.makeDroppable(advanced, {
  acceptDrop: (o) => o instanceof DG.Column,
  acceptDrag: (args) => {
    if (args.dragObject instanceof DG.Column && args.dragObject.type === DG.TYPE.STRING) {
      args.handled = true; // Veto: no string columns allowed.
      return false;
    }
    return true;
  },
  dropSuggestion: 'Drop a numeric column',
  onBeginDrag: (args) => console.log(`begin drag: ${args.dragObjectType}`),
  onEndDrag: (_) => console.log('end drag'),
  onMouseEnter: (_, args) => advanced.style.borderColor = 'var(--blue-1)',
  onMouseLeave: (_, args) => advanced.style.borderColor = 'var(--grey-3)',
  doDrop: (args) => {
    const col = args.dragObject;
    advanced.textContent = `Dropped ${col.name} (type=${col.type}, copying=${args.copying})`;
    advanced.style.borderColor = 'var(--grey-3)';
  }
});

let image = ui.image(
  'https://datagrok.ai/img/heinlein.jpg',
  200, 200,
  {target: 'https://datagrok.ai', tooltipMsg: 'Drag Heinlein to the first field'}
);

ui.makeDraggable(image, {
  getDragObject: () => image,
  getDragCaption: () => 'You are dragging Heinlein!'
});

ui.dialog('Drag and Drop')
  .add(simple)
  .add(DG.Viewer.grid(grok.data.demo.demog()))
  .add(advanced)
  .add(image)
  .show();
