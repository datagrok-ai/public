import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';


export function hasTag(colTags: string[], colTagValue: string): boolean {
  for (let i = 0; i < colTags.length; i++) {
    for (let j = 0; j < colTags[i].length; j++) {
      console.log('TAG Value: ' + colTags[i][j]);
      if (colTags[i][j] == colTagValue)
        return true;
    }
  }
  return false;
}

category('Grid', () => {
  let v: DG.TableView;
  let grid: DG.Grid;
  const demog = grok.data.demo.demog(1000);
  demog.columns.byName('study').name = '~study';

  before(async () => {
    v = grok.shell.addTableView(demog);
    grid = v.grid;
  });

  test('setOrder', async () => {
    grid.columns.setOrder(['race', 'age']);
    const firstCol = grid.columns.byIndex(4);
    const secondCol = grid.columns.byIndex(6);
    if (firstCol?.name != 'race' || secondCol?.name != 'age')
      
      throw new Error('grid.setOrder does not work');
  });

  test('resizeColumn', async () => {
    grid.columns.byName('age')!.width = 200;
    expect(grid.columns.byName('age')!.width, 200);
  });

  test('filter', async () => {
    demog.rows.match('sex = M').filter();
    if (demog.filter.trueCount != 605)
      
      throw new Error('Filtering error');
  });

  test('colorCoding', async () => {
    grid.col('race')!.categoryColors = {
      'Asian': 0xFF0000FF,
      'Black': 0xFF800080,
      'Caucasian': 0xFF05754A,
      'Other': 0XFFE4DD47,
    };

    demog.col('height')!.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Conditional';
    demog.col('height')!.tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = `{"20-170":"#00FF00","170-190":"#220505"}`;

    demog.col('age')!.tags[DG.TAGS.COLOR_CODING_TYPE] = 'Linear';
    demog.col('age')!.tags[DG.TAGS.COLOR_CODING_LINEAR] = `[${DG.Color.orange}, ${DG.Color.green}]`;

    //categorical RACE column check
    const raceTags: string[] = Array.from(demog.col('race')!.tags);
    if (!hasTag(raceTags, '.color-coding-categorical') ||
      !hasTag(raceTags, '{"Asian":4278190335,"Black":4286578816,"Caucasian":4278547786,"Other":4293188935}'))
      
      throw new Error('Categorical Color Coding error');

    //numerical HEIGHT column check for Conditional ColorCoding
    const heightTags: string[] = Array.from(demog.col('height')!.tags);
    if (!hasTag(heightTags, '.color-coding-type') ||
      !hasTag(heightTags, 'Conditional') ||
      !hasTag(heightTags, '.color-coding-conditional') ||
      !hasTag(heightTags, '{"20-170":"#00FF00","170-190":"#220505"}'))
      
      throw new Error('Conditional Color Coding error');

    //numerical AGE column check for Linear ColorCoding
    const ageTags: string[] = Array.from(demog.col('age')!.tags);
    if (!hasTag(ageTags, '.color-coding-type') ||
      !hasTag(ageTags, 'Linear') ||
      !hasTag(ageTags, '.color-coding-linear') ||
      !hasTag(ageTags, '[4294944000, 4278255360]'))
      
      throw new Error('Linear Color Coding error');
  });

  test('columnVisibility', async () => {
    const studyColVisible = grid.columns.byName('~study')!.visible;

    grid.columns.setVisible(['age', 'sex', 'race', 'height', 'weight', 'site', 'subj', 'started']);
    const diseaseColVisible = grid.columns.byName('disease')!.visible;

    if (studyColVisible)
      
      throw new Error('Hiding a column by adding ~ to the name doesn\'t work');

    if (diseaseColVisible)
      
      throw new Error('Hiding a column by using columns.setVisible doesn\'t work');
  });

  test('columnControlledValues', async () => {
    demog.col('site')!.tags[DG.TAGS.CHOICES] = '["New York", "Buffalo"]';
    demog.col('site')!.tags[DG.TAGS.AUTO_CHOICES] = 'New York';

    const siteTags: string[] = Array.from(demog.col('site')!.tags);

    if (!hasTag(siteTags, '.choices') ||
      !hasTag(siteTags, '["New York", "Buffalo"]') ||
      !hasTag(siteTags, '.auto-choices') ||
      !hasTag(siteTags, 'New York'))
      
      throw new Error('Column Controlled Values (Choices) error');
  });

  test('GridColumn.renderer', async () => {
    for (const col of demog.columns.numerical)
      expect(grid.col(col.name)?.renderer.cellType, 'number');

    for (const col of demog.columns.categorical)
      expect(grid.col(col.name)?.renderer.cellType, DG.TYPE.STRING);
  });
});
