import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib';

export function molfileWidget(smiles: string) {
  const molfileStr = OCL.Molecule.fromSmiles(smiles).toMolfile();
  const copyButton = ui.button('Copy', () => {
    navigator.clipboard.writeText(molfileStr);
    grok.shell.info('Molfile copied to your clipboard');
  });
  const molfileInput = ui.textInput('', molfileStr);
  (molfileInput.input as HTMLElement).style.height = '300px';
  //TODO: show copy and reset icons on hover over host
  return new DG.Widget(ui.divV([molfileInput.root, copyButton]));
}
