const typeAhead = ui.typeAhead('Countries', {
  source: {
    local: ['USA', 'Ukraine', 'Antigua', 'United Kingdom', 'United Arab Emirates']
  },
  minLength: 2, 
  limit: 3, 
  hint: true, 
  autoSelect: true, 
  highlight: true, 
  diacritics: true,
  onSubmit: (event, value) => grok.shell.info(value.label),
  debounceRemote: 100
});
  
typeAhead.onInput.subscribe(() => grok.shell.info('input event'));
typeAhead.onChanged.subscribe(() => grok.shell.info('change event'));

grok.shell.newView('TypeAhead', [typeAhead.root]);