// Test rendering a single molecule to a canvas with an optional SMARTS-scaffold

(async () => {
  let root = ui.div();
  const width = 300; const height = 200;
  let canvas1 = document.createElement('canvas');
  canvas1.width = width; canvas1.height = height;
  canvas1.id = 'canvas1';
  let canvas2 = document.createElement('canvas');
  canvas2.width = width; canvas2.height = height;
  canvas2.id = 'canvas2';
  root.appendChild(canvas1);
  root.appendChild(canvas2);
  // TODO: should be syncronous
  await grok.chem.canvasMol(
    0, 0, width, height, canvas1,
    'COc1ccc2cc(ccc2c1)C(C)C(=O)OCCCc3cccnc3',
    'c1ccccc1');
  await grok.chem.canvasMol(
    0, 0, width, height, canvas2,
    'CN1CCC(O)(CC1)c2ccccc2');
  ui
    .dialog({title: 'Test molecule'})
    .add(root)
    .show();
})();