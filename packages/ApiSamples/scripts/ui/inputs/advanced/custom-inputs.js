// Using custom inputs implemented in JS
// 'foo' editor is defined in the Widgets package.

let medication = {name: 'Aspirin', quantity: '20 mg'};
let nameProp = DG.Property.js('name', DG.TYPE.STRING);
let fooProp = DG.Property.fromOptions({name: 'quantity', type: DG.TYPE.STRING, semType: 'foo'});

ui.dialog()
  .add(ui.input.form(medication, [nameProp, fooProp]))
  .show();

