// https://datagrok.ai/help/visualize/viewers/markup

let html = await fetch('/demo/markup_viewer/eeg_21_10-20.svg')
let t = await grok.data.getDemoTable("sensors/eeg.csv");

let view = grok.shell.addTableView(t);
view.markup({content: html});