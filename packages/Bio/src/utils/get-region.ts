import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {getRegion} from '../package';
import {TaskBarProgressIndicator} from 'datagrok-api/dg';

export function getRegionUI(col: DG.Column<string>): void {
  const uh = UnitsHandler.getOrCreate(col);

  const nameInput = ui.stringInput('Name', '');
  const startPositionInput = ui.choiceInput('Start Position', uh.posList[0], uh.posList,
    () => { /* TODO: update name placeholder with getDefaultName() */ });
  const endPositionInput = ui.choiceInput('End Position', uh.posList[uh.posList.length], uh.posList,
    () => { /* TODO: update name placeholder with getDefaultName() */ });

  const getDefaultName = (): string => {
    return `${col.name}:${startPositionInput.value}-${endPositionInput.value}`;
  };

  ui.dialog({title: 'Get Region'}).add(ui.inputs([
    nameInput,
    startPositionInput,
    endPositionInput,
  ])).onOK(() => {
    const pi = TaskBarProgressIndicator.create('Getting region...');
    try {
      const name: string = nameInput.value ?? getDefaultName();
      const regCol = getRegionDo(col, name, startPositionInput.value, endPositionInput.value);
      col.dataFrame.columns.add(regCol);
      regCol.setTag(DG.TAGS.CELL_RENDERER, 'sequence');
    } catch (err: any) {
      grok.shell.error(err.toString());
    } finally { pi.close(); }
  });
}

/** {@link startPosName} and {@link endPosName} are according positionNames tag (or default ['1', '2',...]) */
export function getRegionDo(
  col: DG.Column<string>, startPosName: string | null, endPosName: string | null, name: string | null
): DG.Column<string> {
  const uh = UnitsHandler.getOrCreate(col);

  let startPosIdx: number | null = null;
  let endPosIdx: number | null = null;

  for (let posJ: number = 0; posJ < uh.posList.length; ++posJ) {
    if (uh.posList[posJ] == startPosName) startPosIdx = posJ;
    if (uh.posList[posJ] == endPosName) endPosIdx = posJ;
  }
  if (startPosIdx === null && startPosName !== null)
    throw new Error(`Start position ${startPosName} not found.`);
  if (endPosIdx === null && endPosName !== null)
    throw new Error(`End position ${endPosName} not found.`);

  if (uh.posList.length < endPosIdx!)
    throw new Error(`End position ${endPosIdx} exceeds positions length`);

  const regColName: string = !!name ? name : `${col.name}: (${startPosName ?? ''}-${endPosName ?? ''})`;

  const regCol = uh.getRegion(startPosIdx, endPosIdx, regColName);
  return regCol;
}
