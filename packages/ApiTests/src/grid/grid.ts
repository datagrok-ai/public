import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';


export function hasTag(colTags: string[], colTagValue: string): boolean {
  for (let i = 0; i < colTags.length; i++) {
    for (let j = 0; j < colTags[i].length; j++) {
      // console.log('TAG Value: ' + colTags[i][j]);
      if (colTags[i][j] == colTagValue)
        return true;
    }
  }
  return false;
}

category('Grid', () => {
  let grid: DG.Grid;
  const demog = grok.data.demo.demog(100);
  demog.columns.byName('study').name = '~study';

  before(async () => {
    grid = demog.plot.grid();
  });

  test('setOrder', async () => {
    grid.columns.setOrder(['race', 'age']);
    const firstCol = grid.columns.byIndex(4);
    const secondCol = grid.columns.byIndex(6);
    expect(firstCol?.name, 'race');
    expect(secondCol?.name, 'age');
  });

  test('resizeColumn', async () => {
    grid.columns.byName('age')!.width = 200;
    expect(grid.columns.byName('age')!.width, 200);
  });

  test('filter', async () => {
    demog.rows.match('sex = M').filter();
    expect(demog.filter.trueCount, 73);
  });

  test('colorCoding', async () => {
    grid.col('race')!.categoryColors = {
      'Asian': 0xFF0000FF,
      'Black': 0xFF800080,
      'Caucasian': 0xFF05754A,
      'Other': 0XFFE4DD47,
    };

    demog.col('height')?.meta.colors.setConditional({'20-170': '#00FF00', '170-190': '#220505'});
    demog.col('age')?.meta.colors.setLinear([DG.Color.orange, DG.Color.green]);

    //categorical RACE column check
    const raceCol = demog.col('race')!;
    if (raceCol.getTag(DG.TAGS.COLOR_CODING_CATEGORICAL) !== '{"Asian":4278190335,"Black":4286578816,"Caucasian":4278547786,"Other":4293188935}')
      throw new Error('Categorical Color Coding error');

    //numerical HEIGHT column check for Conditional ColorCoding
    const heightCol = demog.col('height')!;
    if (heightCol.meta.colors.getType() !== DG.COLOR_CODING_TYPE.CONDITIONAL ||
      heightCol.getTag(DG.TAGS.COLOR_CODING_CONDITIONAL) !== '{"20-170":"#00FF00","170-190":"#220505"}')
      throw new Error('Conditional Color Coding error');

    //numerical AGE column check for Linear ColorCoding
    const ageCol = demog.col('age')!;
    if (ageCol.meta.colors.getType() !== DG.COLOR_CODING_TYPE.LINEAR ||
      ageCol.getTag(DG.TAGS.COLOR_CODING_LINEAR) !== '[4294944000,4278255360]')
      throw new Error('Linear Color Coding error');
  });

  test('columnVisibility', async () => {
    const studyColVisible = grid.columns.byName('~study')!.visible;
    grid.columns.setVisible(['age', 'sex', 'race', 'height', 'weight', 'site', 'subj', 'started']);
    const diseaseColVisible = grid.columns.byName('disease')!.visible;
    expect(studyColVisible, false, 'Hiding a column by adding ~ to the name doesn\'t work');
    expect(diseaseColVisible, false, 'Hiding a column by using columns.setVisible doesn\'t work');
  });

  test('columnControlledValues', async () => {
    const col = demog.col('site');
    col!.meta.choices = ['New York', 'Buffalo'];
    col!.meta.autoChoices = false;

    if (JSON.stringify(col?.meta.choices) !== '["New York","Buffalo"]' || col?.meta.autoChoices !== false)
      throw new Error('Column Controlled Values (Choices) error');
  });

  test('GridColumn.renderer', async () => {
    for (const col of demog.columns.numerical)
      expect(grid.col(col.name)?.renderer.cellType, 'number');

    for (const col of demog.columns.categorical)
      expect(grid.col(col.name)?.renderer.cellType, DG.TYPE.STRING);
  });

  test('getOptions', async () => {
    expect(Object.keys(grid.getOptions().look).length, 2);
  });

  test('setOptions', async () => {
    grid.setOptions({allowEdit: false, showColumnLabels: false, colHeaderHeight: 100});
    expect(Object.keys(grid.getOptions().look).length, 5);
  });
});
