let p = DG.ProgressIndicator.create();
p.onProgressUpdated.subscribe((p) => grok.shell.info(p.percent));
p.onLogUpdated.subscribe((line) => console.log(line));
p.update(50, "Loading foo.csv");

grok.functions.call("Samples:RDup",
  {s: "hi"}, true, p)
  .then((res) => {
    grok.shell.info(res);
  });