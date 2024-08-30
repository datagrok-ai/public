// Manual handling of input-originated events

const medication = {name: 'Aspirin', quantity: '20 mg'};

// Core value input implemented in Dart
const nameProp = DG.Property.js('name', DG.TYPE.STRING);
const nameEditor = DG.InputBase.forProperty(nameProp, medication);
nameEditor.onChanged.subscribe((value) => grok.shell.info('name changed: ' + value));
nameEditor.onInput.subscribe((_) => grok.shell.info('name input: ' + nameEditor.value));

// 'foo' editor is defined in the Widgets package.
const fooProp = DG.Property.fromOptions({name: 'quantity', type: DG.TYPE.STRING, semType: 'foo'});
const quantityEditor = DG.InputBase.forProperty(fooProp, medication);
quantityEditor.onChanged.subscribe((value) => grok.shell.info('quantity changed: ' + value));
quantityEditor.onInput.subscribe((_) => grok.shell.info('quantity input: ' + quantityEditor.value));

quantityEditor.onInput.subscribe((_) => grok.shell.o = medication);

ui.dialog()
  .add(ui.form([nameEditor, quantityEditor]))
  .show();

