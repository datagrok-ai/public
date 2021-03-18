// Programmatically resolving an ObjectHandler

class AppleHandler extends DG.ObjectHandler {
  get type() { return 'fruit / apple' }
}

class OrangeHandler extends DG.ObjectHandler {
  get type() { return 'fruit / orange' }
}

DG.ObjectHandler.onResolve((args) => {
  if (args.object === 'apple')
    args.handler = new AppleHandler();
  if (args.object === 'orange')
    args.handler = new OrangeHandler();
})

let meta = DG.ObjectHandler.forEntity('orange');
grok.shell.info(meta.type);