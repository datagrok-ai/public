grok.events.onBrowseNodeCreated
  .pipe(rxjs.operators.filter(node => node.text === 'Spaces'))
  .subscribe(node => node.item('Custom node'));