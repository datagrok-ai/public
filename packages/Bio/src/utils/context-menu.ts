import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

import {_package} from '../package';


export function addCopyMenuUI(cell: DG.Cell, menu: DG.Menu): void {
  const uh = UnitsHandler.getOrCreate(cell.column);
  const tgtNotationList: string[] = Object.values(NOTATION).filter((v) => v !== uh.units);

  menu.group('Copy')
    .items(tgtNotationList, (tgtNotation) => {
      const srcCol = cell.column;
      const srcRowIdx = cell.rowIndex;
      const srcUh = UnitsHandler.getOrCreate(srcCol);
      const separator = tgtNotation === NOTATION.SEPARATOR ? _package.properties.DefaultSeparator : undefined;
      const joiner = srcUh.getJoiner({notation: tgtNotation as NOTATION, separator});
      const srcSS = srcUh.splitted[srcRowIdx];
      const tgtSeq = joiner(srcSS);


      if (!navigator.clipboard) {
        grok.shell.warning('The clipboard functionality requires a secure origin — either HTTPS or localhost');
      } else {
        navigator.clipboard.writeText(tgtSeq);
        grok.shell.info(`Value of notation '${tgtNotation}' copied to clipboard`);
      }
    });
}
