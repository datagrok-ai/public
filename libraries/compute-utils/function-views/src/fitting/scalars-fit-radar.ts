/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {INDICES, MIN_RADAR_COLS_COUNT, NAME} from './constants';

export class ScalarsFitRadar {
  private table: DG.DataFrame;

  constructor(table: DG.DataFrame) {
    this.validate(table);
    this.table = table.clone();
    this.prepareTable();
  }

  public async getRoot(): Promise<HTMLElement> {
    return await this.getViewerRoot();
  }

  private async getViewerRoot(): Promise<HTMLElement> {
    const radar = await this.table.plot.fromType('Radar', {
      colorColumnName: NAME.CATEGORY,
      currentRowColor: DG.Color.categoricalPalette[INDICES.LIGHT_GREY],
      showValue: true,
      resizeScheduled: true,
    });

    return radar.root;
  }

  private validate(table: DG.DataFrame): void {
    const cols = table.columns;
    const colsCount = cols.length;

    if (colsCount < MIN_RADAR_COLS_COUNT)
      throw new Error(`Incorrect number of column for radar viewer: ${colsCount}. Expected: > ${MIN_RADAR_COLS_COUNT - 1}`);

    const catColType = cols.byIndex(colsCount - 1).type;
    if (catColType !== DG.COLUMN_TYPE.STRING)
      throw new Error(`Incorrect category column type of a table for radar viewer: ${catColType}. Expected: ${DG.COLUMN_TYPE.STRING}`);

    for (let i = 0; i < colsCount - 2; ++i) {
      const col = cols.byIndex(i);

      if (!col.isNumerical)
        throw new Error(`Incorrect values column '${col.name}' of a table for radar viewer: ${col.type}. Expected: numerical`);
    }
  }

  private prepareTable(): void {
    const catCol = this.table.col(NAME.CATEGORY);

    if (catCol !== null) {
      catCol.colors.setCategorical({
        'Simulation': DG.Color.categoricalPalette[INDICES.BLUE],
        'Target': DG.Color.categoricalPalette[INDICES.ORANGE],
      });
    }
  }
}
