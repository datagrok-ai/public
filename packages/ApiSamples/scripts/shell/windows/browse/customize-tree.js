grok.events.onBrowseNodeCreated
  .pipe(rxjs.operators.filter(node => node.text === 'Namespaces'))
  .subscribe(node => node.item('Custom node'));