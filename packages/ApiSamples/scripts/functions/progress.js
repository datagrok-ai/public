let p = DG.ProgressIndicator.create();
p.onProgressUpdated.subscribe((p) => grok.shell.info(p.percent));

grok.functions.call("ProgressIndicator:countToTen",
    {s: "a"}, true, p)
    .then((res) => {grok.shell.info(res); });