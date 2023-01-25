import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, after, category, test} from '@datagrok-libraries/utils/src/test';
import {hasTag} from './grid';


category('Grid: Color Coding', () => {
  let v: DG.TableView;
  let grid: DG.Grid;
  let demog = grok.data.demo.demog(1000);
  
  before(async () => {
    v = grok.shell.addTableView(demog);
    grid = v.grid;
  });

  test('colorCoding.api', async () => {
    demog.col('age')!.colors.setLinear();
    demog.col('age')!.colors.setConditional({'<30': DG.Color.green, '30-70': '#ff0000'});
    demog.col('sex')!.colors.setCategorical({'M': 0xFF0000FF, 'F': 0xFF800080});
    demog.col('started')!.colors.setLinear([DG.Color.white, DG.Color.red]);
    grid.setOptions({colorCoding: 'All'});
    grid.setOptions({colorCoding: 'None'});
    grid.setOptions({colorCoding: 'Auto'});
    testTags();

    const layout = v.saveLayout();
    v.close();
    grok.shell.closeTable(demog);
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    //grid = v.grid;
    v.loadLayout(layout);
    testTags();
  });

  after(async () => {
    v.close();
    grok.shell.closeTable(demog);
  });

  function testTags() {
    const ageTags: any[] = Array.from(demog.col('age')!.tags);
    if (!hasTag(ageTags, '.color-coding-type') ||
      !hasTag(ageTags, 'Conditional') ||
      !hasTag(ageTags, '.color-coding-conditional') ||
      !hasTag(ageTags, '{"<30":"#00ff00","30-70":"#ff0000"}'))
      throw new Error('Conditional Color Coding error on Age column');

    const sexTags: any[] = Array.from(demog.col('sex')!.tags);
    if (!hasTag(sexTags, '.color-coding-type') ||
      !hasTag(sexTags, 'Categorical') ||
      !hasTag(sexTags, '.color-coding-categorical') ||
      !hasTag(sexTags, '{"M":4278190335,"F":4286578816}'))
      throw new Error('Categorical Color Coding error on Sex column');
   
    const startedTags: any[] = Array.from(demog.col('started')!.tags);
    if (!hasTag(startedTags, '.color-coding-type') ||
      !hasTag(startedTags, 'Linear') ||
      !hasTag(startedTags, '.color-coding-linear') ||
      !hasTag(startedTags, '[4294967295,4294901760]'))
      throw new Error('Linear Color Coding error on Started column');
  }
});
