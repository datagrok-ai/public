// creating a form, and accessing its individual inputs

const sin = DG.Func.find({name: 'Sin'});

let fc = sin[0].prepare({dosage: 77});
DG.InputForm.forFuncCall(fc).then((inputForm) => {
  inputForm.getInput('x').addPostfix('foo');
  inputForm.onInputChanged.subscribe((_) => grok.shell.info('changed'));
  grok.shell.o = inputForm.root;
});
