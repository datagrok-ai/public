import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {ISeqHandler} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {removeGapsFromHelm} from '../utils';

export class SeqPropertiesError extends Error {
  constructor(message?: string, options?: ErrorOptions) {
    super(message, options);
  }
}

// hwe migration (Phase 7): formula / MW / extinction coefficient now come from
// the standalone `@datagrok-libraries/hwe` editor (via the adapter's
// `LegacyEditorWrapper`), replacing the legacy `new JSDraw2.Editor(...)`. The
// helper is obtained through the async `getHelmHelper()` factory (which ensures
// the Helm package is initialized) rather than `_package.helmHelper` — the
// latter is not initialized in the standalone test bundle. This makes the
// property calc async (it previously relied on the eagerly-loaded JSDraw2
// global). The adapter emits the formula in the legacy element order
// (C,H,N,O,S,P) with no `<sub>` markup, so the strip below is a defensive no-op.
export async function getPropertiesDict(seq: string, sh: ISeqHandler): Promise<{}> {
  const host = ui.div([], {style: {width: '0', height: '0'}});
  document.documentElement.appendChild(host);
  try {
    if (seq.length > 1000)
      throw new SeqPropertiesError('Too long sequence to calculate molecular properties.');

    const helmString = sh.isHelm() ? seq : sh.getConverter(NOTATION.HELM)(seq);
    // Remove gap symbols for compatibility with the editor
    const helmString2 = removeGapsFromHelm(helmString);

    const helmHelper = await getHelmHelper();
    const editor = helmHelper.createHelmWebEditor(host).editor;
    editor.setHelm(helmString2);

    const formula = editor.getFormula(true);
    const molWeight = Math.round(editor.getMolWeight() * 100) / 100;
    const coef = Math.round(editor.getExtinctionCoefficient() * 100) / 100;
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


export async function getPropertiesWidget(value: DG.SemanticValue<string>): Promise<DG.Widget> {
  const helmHelper = await getHelmHelper();
  const sh = helmHelper.seqHelper.getSeqHandler(value.cell.column as DG.Column<string>);
  try {
    const propDict = await getPropertiesDict(value.value, sh);
    return new DG.Widget(
      ui.tableFromMap(propDict, true)
    );
  } catch (err: any) {
    if (err instanceof SeqPropertiesError) {
      // Display warning from err.message
      return new DG.Widget(ui.divText(err.message));
    }
    throw err;
  }
}
