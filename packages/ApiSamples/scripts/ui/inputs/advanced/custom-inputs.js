// Using custom inputs implemented in JS
// 'foo' editor is defined in the Widgets package.

var medication = { name: 'Aspirin', quantity: '20 mg'};
var nameProp = DG.Property.js('name', DG.TYPE.STRING);
var fooProp = DG.Property.fromOptions({name: 'quantity', type: DG.TYPE.STRING, semType: 'foo'});

ui.dialog()
  .add(ui.input.form(medication, [nameProp, fooProp]))
  .show();

