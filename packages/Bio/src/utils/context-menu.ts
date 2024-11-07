import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {_package} from '../package';


export function addCopyMenuUI(cell: DG.Cell, menu: DG.Menu, seqHelper: ISeqHelper): void {
  const sh = seqHelper.getSeqHandler(cell.column);
  const tgtNotationList: string[] = Object.values(NOTATION).filter((v) => v !== sh.units);

  menu.group('Copy')
    .items(tgtNotationList, (tgtNotation) => {
      const srcCol = cell.column;
      const srcRowIdx = cell.rowIndex;
      const srcSh = seqHelper.getSeqHandler(srcCol);
      const separator = tgtNotation === NOTATION.SEPARATOR ? _package.properties.defaultSeparator : undefined;
      const joiner = srcSh.getJoiner({notation: tgtNotation as NOTATION, separator});
      const srcSS = srcSh.getSplitted(srcRowIdx);
      const tgtSeq = joiner(srcSS);

      if (!navigator.clipboard) {
        grok.shell.warning('The clipboard functionality requires a secure origin — either HTTPS or localhost');
      } else {
        navigator.clipboard.writeText(tgtSeq);
        grok.shell.info(`Value of notation '${tgtNotation}' copied to clipboard`);
      }
    });
}
