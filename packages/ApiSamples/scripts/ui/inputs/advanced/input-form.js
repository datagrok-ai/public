// creating a form, and accessing its individual inputs

DG.Func.findAll({name: 'Pkpd'}).then((fs) => {
  let fc = fs[0].prepare({dosage: 77});
  DG.InputForm.forFuncCall(fc).then((inputForm) => {
    inputForm.getInput('Clearance').addPostfix('foo');
    inputForm.onInputChanged.subscribe((_) => grok.shell.info('changed'));
    grok.shell.o = inputForm.root;
  });
});