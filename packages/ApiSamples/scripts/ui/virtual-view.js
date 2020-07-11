// View creates DOM elements only when necessary.

let data = grok.data.testData('demog', 100000);

grok.shell.newView('virtual view', [
    ui.virtualView(data.rowCount, (i) => ui.card(ui.tableFromMap({
        id: data.row(i).subj,
        age: data.row(i).age,
        sex: data.row(i).sex,
        height: data.row(i).height,
        weight: data.row(i).weight
    }))).root
]);
