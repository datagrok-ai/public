// See ProgressIndicator:countToTen function:
// https://github.com/datagrok-ai/public/blob/a53c5606186eb88eb9085559bb2568aabd79c729/packages/ProgressIndicator/src/package.js#L12

let p = DG.ProgressIndicator.create();
p.onProgressUpdated.subscribe((p) => grok.shell.info(p.percent));

grok.functions.call("ProgressIndicator:countToTen",
    {s: "a"}, true, p)
    .then((res) => {grok.shell.info(res); });