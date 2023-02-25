// user-editable form

let t = grok.data.demo.demog();
let form = DG.Form.forDataFrame(t, {columns: ['sex', 'race', 'age']});
form.editable = true;
form.designMode = true;
form.row = new DG.Row(t, 4);

ui.dialog().add(form).show();