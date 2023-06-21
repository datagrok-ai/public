// https://datagrok.ai/help/visualize/viewers/form

let form = DG.FormViewer.createDefault(grok.data.demo.demog(), {columns: ['sex', 'race', 'age']});
form.editable = true;
form.designMode = true;

ui.dialog().add(form).show();
