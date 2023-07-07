let handler = DG.ObjectHandler.list().find(h => h.type === 'brief handler');

DG.ObjectHandler.onResolve((args) => {
  if (args.semType === 'plate')
    args.handler = handler;
});

let h = DG.ObjectHandler.forEntity(DG.SemanticValue.fromValueType('2134', 'plate'));
console.log(h);