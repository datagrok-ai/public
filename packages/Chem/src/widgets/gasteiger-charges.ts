import * as grok from 'datagrok-api/grok'
import * as ui from  'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';

export async function gasteigerChargesWidget(molString: string): Promise<DG.Widget> {
  const rdKitModule = getRdKitModule();
  try {
    molString = _convertMolNotation(molString, 'unknown', 'molblock', rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  const parameters = {mol: molString, contours: 10};
  let element: string;
  try {
    element = await grok.functions.call('Chem:GasteigerCharges', parameters);
  } catch (e) {
    return new DG.Widget(ui.divText('Couldn\'t visualize molecule'));
  }
  return new DG.Widget(ui.image(`data:image/png;base64,${element}`, 300, 150));
}
