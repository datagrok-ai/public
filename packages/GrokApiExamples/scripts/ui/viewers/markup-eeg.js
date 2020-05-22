// https://datagrok.ai/help/viewers/markup-viewer

fetch('/demo/markup_viewer/eeg_21_10-20.svg').then((r) => r.text().then((html) =>
    grok.data.getDemoTable("sensors/eeg.csv").then(function (t) {
        let view = grok.shell.addTableView(t);
        view.markup({content: html});
    })));
