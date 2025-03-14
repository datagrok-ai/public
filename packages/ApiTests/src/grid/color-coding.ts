import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/utils/src/test';
import {hasTag} from './grid';


category('Grid: Color Coding', () => {
  let v: DG.TableView;
  let grid: DG.Grid;
  let demog = grok.data.demo.demog(1000);

  test('colorCoding.api', async () => {
    v = grok.shell.addTableView(demog);
    grid = v.grid;
    demog.col('age')!.meta.colors.setLinear();
    demog.col('age')!.meta.colors.setConditional({'<30': DG.Color.green, '30-70': '#ff0000'});
    demog.col('sex')!.meta.colors.setCategorical({'M': 0xFF0000FF, 'F': 0xFF800080});
    demog.col('started')!.meta.colors.setLinear([DG.Color.white, DG.Color.red]);
    demog.col('height')!.meta.colors.setLinearAbsolute({58.31: '#73aff5', 97.98: '#ff5140', 137.65: '#ffa500', 177.32: '#50af28', 197.16: '#d9d9d9', 217: '#9467bd'});
    grid.setOptions({colorCoding: 'None'});
    grid.setOptions({colorCoding: 'Auto'});
    testTags();
    grid.setOptions({colorCoding: 'All'});
    grid.setOptions({colorScheme: [4294922560, 4294944000, 4283477800]});
    if (grid.cell('weight', 0).color !== 4288522774)
      throw new Error('Grid Color Coding > Color Scheme was not applied');
    const layout = v.saveLayout();
    v.close();
    grok.shell.closeTable(demog);
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    v.loadLayout(layout);
    testTags(' after layout');
  }, {owner: 'dkovalyov@datagrok.ai', skipReason: 'GROK-17728'});

  function testTags(after: string = '') {
    const ageCol = demog.col('age')!;
    if (ageCol.meta.colors.getType() !== DG.COLOR_CODING_TYPE.CONDITIONAL ||
      ageCol.getTag(DG.TAGS.COLOR_CODING_CONDITIONAL) !== '{"<30":"#00ff00","30-70":"#ff0000"}')
      throw new Error('Conditional Color Coding error on Age column' + after);

    const sexCol = demog.col('sex')!;
    if (sexCol.meta.colors.getType() !== DG.COLOR_CODING_TYPE.CATEGORICAL ||
      sexCol.getTag(DG.TAGS.COLOR_CODING_CATEGORICAL) !== '{"M":4278190335,"F":4286578816}')
      throw new Error('Categorical Color Coding error on Sex column' + after);

    const startedCol = demog.col('started')!;
    if (startedCol.meta.colors.getType() !== DG.COLOR_CODING_TYPE.LINEAR ||
      startedCol.getTag(DG.TAGS.COLOR_CODING_LINEAR) !== '[4294967295,4294901760]')
      throw new Error('Linear Color Coding error on Started column' + after);

    const heightCol = demog.col('height')!;
    if (heightCol.meta.colors.getType() !== DG.COLOR_CODING_TYPE.LINEAR ||
      heightCol.getTag(DG.TAGS.COLOR_CODING_LINEAR) !== '[4285771765,4294922560,4294944000,4283477800,4292467161,4287915965]' ||
      heightCol.getTag(DG.TAGS.COLOR_CODING_LINEAR_IS_ABSOLUTE) !== 'true' ||
      heightCol.getTag(DG.TAGS.COLOR_CODING_LINEAR_ABSOLUTE) !== '{"58.31":"#73aff5","97.98":"#ff5140","137.65":"#ffa500","177.32":"#50af28","197.16":"#d9d9d9","217":"#9467bd"}')
      throw new Error('Linear Color Coding with absolute values error on Height column' + after);
  }
});
