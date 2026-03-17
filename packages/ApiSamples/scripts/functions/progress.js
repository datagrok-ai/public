let p = DG.ProgressIndicator.create();
p.onProgressUpdated.subscribe((p) => grok.shell.info(p.percent));
p.onLogUpdated.subscribe((line) => console.log(line));
p.update(50, 'Loading foo.csv');

let res = await grok.functions.call('ApiSamples:RDup', {s: 'hi'}, true, p);
grok.shell.info(res); 