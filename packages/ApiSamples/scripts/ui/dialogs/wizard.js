// Creating custom dialogs, and add the help url

let name = ui.input.string('Title');
let type = ui.input.choice('Type', {items: ['Bio', 'Chem'], value: 'Bio'});

new DG.Wizard({title: 'Onboard a model'})
  .page({
    caption: 'Details',
    root: ui.inputs([name, type]),
    validate: () => name.value !== '' ? null : 'Please provide a name'
  })
  .page({caption: 'Model Source', root: ui.divText('Second page')})
  .page({caption: 'Finalize', root: ui.divText('Finish')})
  .onOK(() => {})
  .show({width: 500, height: 300});
