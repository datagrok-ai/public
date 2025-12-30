// Synchronous programmatic ObjectHandler resolution.
// See dynamic-loading.js for dynamic resolution based on the semType

class AppleHandler extends DG.ObjectHandler {
  get type() {return 'fruit / apple';}
}

class OrangeHandler extends DG.ObjectHandler {
  get type() {return 'fruit / orange';}
}

DG.ObjectHandler.onResolve((args) => {
  if (args.value === 'apple')
    args.handler = new AppleHandler();
  if (args.value === 'orange')
    args.handler = new OrangeHandler();
});

let meta = DG.ObjectHandler.forEntity('orange');
grok.shell.info(meta.type);