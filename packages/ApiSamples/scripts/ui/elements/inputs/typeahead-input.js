const typeAhead = ui.typeAhead('Countries', {source: {local: ['USA', 'Ukraine', 'Antigua', 'United Kingdom', 'United Arab Emirates']},
    minLength: 2, limit: 3, hint: true, autoSelect: true, highlight: true, diacritics: true,
    onSubmit: (event, value) => console.log(`The value is ${value?.label}`), debounceRemote: 100});

typeAhead.onInput(() => console.log('input event'));
typeAhead.onChanged(() => console.log('change event'));

grok.shell.newView('TypeAhead', [ui.ribbonPanel([typeAhead.root])]);
