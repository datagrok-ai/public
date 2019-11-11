// https://datagrok.ai/help/viewers/histogram

let view = gr.addTableView(gr.testData('demog', 5000));

view.histogram({
    value: 'age'
});