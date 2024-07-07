/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';

import {ALPHABET, ALIGNMENT, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

export function _setPeptideColumn(col: DG.Column): void {
  addCommonTags(col);
  col.meta.units = NOTATION.SEPARATOR;
  col.setTag('separator', '-');
  // col.setTag('cell.renderer', 'sequence');
}

function addCommonTags(col: DG.Column<any>) {
  col.setTag('quality', DG.SEMTYPE.MACROMOLECULE);
  col.setTag('aligned', ALIGNMENT.SEQ);
  col.setTag('alphabet', ALPHABET.PT);
}

export function handleError(err: any): void {
  const errMsg: string = err instanceof Error ? err.message : err.toString();
  const stack: string | undefined = err instanceof Error ? err.stack : undefined;
  grok.shell.error(errMsg);
  _package.logger.error(err.message, undefined, stack);
}
