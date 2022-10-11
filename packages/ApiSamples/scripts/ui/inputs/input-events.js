// Manual handling of input-originated events

var medication = { name: 'Aspirin', quantity: '20 mg'};

// Core value input implemented in Dart
var nameProp = DG.Property.js('name', DG.TYPE.STRING);
var nameEditor = DG.InputBase.forProperty(nameProp, medication);
nameEditor.onChanged((_) => grok.shell.info('name changed: ' + nameEditor.value));
nameEditor.onInput((_) => grok.shell.info('name input: ' + nameEditor.value));

// 'foo' editor is defined in the Widgets package.
var fooProp = DG.Property.fromOptions({name: 'quantity', type: DG.TYPE.STRING, semType: 'foo'});
var quantityEditor = DG.InputBase.forProperty(fooProp, medication);
quantityEditor.onChanged((_) => grok.shell.info('quantity changed: ' + quantityEditor.value));
quantityEditor.onInput((_) => grok.shell.info('quantity input: ' + quantityEditor.value));

quantityEditor.onInput((_) => grok.shell.o = medication);

ui.dialog()
  .add(ui.form([nameEditor, quantityEditor]))
  .show();

