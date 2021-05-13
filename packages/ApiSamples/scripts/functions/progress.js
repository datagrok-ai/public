//tags: Func
// See Discovery:countToTen function:
// https://github.com/datagrok-ai/public/blob/master/packages/Discovery/src/package.js

let p = DG.ProgressIndicator.create();
p.onProgressUpdated.subscribe((p) => grok.shell.info(p.percent));
p.onLogUpdated.subscribe((line) => console.log(line));

grok.functions.call("Discovery:countToTen",
  {s: "a"}, true, p)
  .then((res) => {
    grok.shell.info(res);
  });