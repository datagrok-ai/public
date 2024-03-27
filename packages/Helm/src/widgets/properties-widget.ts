import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

import {removeGapsFromHelm} from '../utils';

export function getPropertiesDict(seq: string, uh: UnitsHandler): {} {
  const helmString = uh.isHelm() ? seq : uh.getConverter(NOTATION.HELM)(seq);
  // Remove gap symbols for compatibility with WebEditor
  const helmString2 = removeGapsFromHelm(helmString);

  const host = ui.div([], {style: {width: '0', height: '0'}});
  document.documentElement.appendChild(host);
  try {
    const editor = new JSDraw2.Editor(host, {viewonly: true});
    editor.setHelm(helmString2);

    const formula = editor.getFormula(true);
    const molWeight = Math.round(editor.getMolWeight() * 100) / 100;
    const coef = Math.round(editor.getExtinctionCoefficient(true) * 100) / 100;
    return {
      'formula': formula.replace(/<sub>/g, '').replace(/<\/sub>/g, ''),
      'molecular weight': molWeight,
      'extinction coefficient': coef,
    };
  } finally {
    $(host).empty();
    host.remove();
  }
}


export function getPropertiesWidget(value: DG.SemanticValue): DG.Widget {
  const uh = UnitsHandler.getOrCreate(value.cell.column);
  const propDict = getPropertiesDict(value.value, uh);
  return new DG.Widget(
    ui.tableFromMap(propDict)
  );
}
