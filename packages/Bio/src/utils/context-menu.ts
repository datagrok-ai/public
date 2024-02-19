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
      const ncUH = UnitsHandler.getOrCreate(cell.column);
      const separator = tgtNotation === NOTATION.SEPARATOR ? _package.properties.DefaultSeparator : undefined;
      const converter = ncUH.getConverter(tgtNotation as NOTATION, separator);
      const tgtSeq = converter(cell.value);

      if (!navigator.clipboard) {
        grok.shell.warning('The clipboard functionality requires a secure origin — either HTTPS or localhost');
      } else {
        navigator.clipboard.writeText(tgtSeq);
        grok.shell.info(`Value of notation '${tgtNotation}' copied to clipboard`);
      }
    });
}
