// View creates DOM elements only when necessary.

let data = gr.testData('demog', 100000);

gr.newView('virtual view', [
    ui.virtualView(data.rowCount, (i) => ui.card(ui.tableFromMap({
        id: data.row(i).subj,
        age: data.row(i).age,
        sex: data.row(i).sex,
        height: data.row(i).height,
        weight: data.row(i).weight
    })))
]);