// Take over file drag-and-drop: prevent the platform's drop overlay and handle drops yourself

let sub = grok.events.onFileDragEnter.subscribe((ev) => {
  ev.preventDefault(); // the platform overlay stays away while this view owns the drop
  grok.shell.info('File drag detected — attach your own drop handlers');
});

// With the overlay suppressed, plain DOM handlers on your view receive the files:
// myView.root.addEventListener('dragover', (e) => e.preventDefault());
// myView.root.addEventListener('drop', (e) => console.log(e.dataTransfer.files));
