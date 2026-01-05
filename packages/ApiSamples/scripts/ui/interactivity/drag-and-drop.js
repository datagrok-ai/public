// Demonstrates how to create and manipulate drag-n-dropped objects.

let control = ui.input.search('', {value: 'Drag column header or image here...'});

// Allow dragging and dropping objects into the search field:
ui.makeDroppable(control.input, {
  acceptDrop: (_) => true, // Allow to drag anything into the field
  doDrop: (o, _) => { // Let's check what fell into the input field
    if (o instanceof DG.Column)
      control.value = o.name.toUpperCase() + ' was dropped';
    else if (o instanceof HTMLDivElement)
      control.value = o.style.backgroundImage;
  }
});

let image = ui.image(
  'https://datagrok.ai/img/heinlein.jpg',
  200, 200,
  {target: 'https://datagrok.ai', tooltipMsg: 'Drag Heinlein to Input field'}
);

// Let's grab and drag the object (image):
ui.makeDraggable(image, {
  getDragObject: () => image,
  getDragCaption: () => 'You are dragging Heinlein!'
});

ui.dialog('Drag and Drop')
  .add(control)
  .add(DG.Viewer.grid(grok.data.demo.demog()))
  .add(image)
  .show();
