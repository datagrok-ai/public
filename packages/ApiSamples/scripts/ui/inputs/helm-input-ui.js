const view = grok.shell.newView('demo: inputs');

const helmValue = 'PEPTIDE1{[D-Pyr].Y.C.A.T.P.F.G.T.S.R.T}$PEPTIDE1,PEPTIDE1,3:R3-12:R2$$$V2.0';
  view.append(ui.inputs([await ui.input.helmAsync('Macromolecule', {value: helmValue})]));
 
  view.append(ui.inputs([await grok.functions.call(
    'Helm:helmInput', {name: 'Macromolecule', options: {value: helmValue}})]));