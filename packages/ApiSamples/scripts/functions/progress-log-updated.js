// A function call in debug mode reports its timing events to the progress indicator's log;
// each event is a plain {level, message, flag, params, time} object
let p = DG.ProgressIndicator.create();
let sub = p.onLogUpdated.subscribe((m) => console.log(`${m.time} ${m.flag}: ${m.message}`));

let call = DG.Script.create('//language: javascript\n//output: int x\nx = 1 + 1;').prepare({});
call.options['debug'] = true;
await call.call(false, p, {processed: false, report: false});
sub.unsubscribe();
grok.shell.info(`x = ${call.outputs.get('x')}`);
