import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {after, awaitCheck, before, category, expect, expectArray, expectFloat,
  expectObject, isDialogPresent, test} from '@datagrok-libraries/utils/src/test';
import {AddNewColumnDialog} from './add-new-column';


category('Dialogs', () => {
  let tv: DG.TableView;
  let df: DG.DataFrame;
  const dialogs: DG.Dialog[] = [];

  before(async () => {
    df = grok.data.demo.demog(1000);
    tv = grok.shell.addTableView(df);
  });

  test('AddNewColumn', async () => {
    const dlg = new AddNewColumnDialog();
    await awaitCheck(() => isDialogPresent(dlg.addColumnTitle));
    dialogs.push(dlg.uiDialog!);
    const funcs = DG.Func.find();

    for (const f of funcs) {
      if (!(f.name in FUNC_TESTS))
        continue;
      for (const [expression, result] of Object.entries(FUNC_TESTS[f.name])) {
        dlg.inputName!.value = expression;
        dlg.inputExpression!.value = expression;
        dlg.uiDialog!.getButton('OK').click();
        await awaitCheck(() => df.col(expression) != null);
        expectTyped(df.col(expression)!.get(0), result);
      }
    }
  });

  test('FormulaLines', async () => {});

  after(async () => {
    dialogs.forEach((d) => d.close());
    tv.close();
    grok.shell.closeTable(df);
  });
});

// See ApiTests/src/functions/
const FUNC_TESTS: {[f: string]: {[test: string]: any}} = {
  Boolean: {
    'Boolean(true)': true,
    'Boolean("true")': true,
    'Boolean("y")': true,
    'Boolean(1)': true,
    'Boolean(10)': true,
    'Boolean(DateParse("20200131T132700"))': true,
    'Boolean(false)': false,
    'Boolean("false")': false,
    'Boolean("n")': false,
    'Boolean("abc")': false,
    'Boolean("")': false,
    'Boolean(null)': false,
    'Boolean(0)': false,
  },
  ParseFloat: {
    'ParseFloat("2025")': 2025,
    'ParseFloat("12.78")': 12.78,
    'ParseFloat("-012.150")': -12.15,
  },
  ParseInt: {
    'ParseInt("2025")': 2025,
    'ParseInt("-012")': -12,
    'ParseInt(" 0101 ")': 101,
  },
  ParseQnum: {
    'ParseQnum("100")': 100,
    'ParseQnum("<100")': 100,
    'ParseQnum(">100")': 100,
    'ParseQnum("   100 ")': 100,
    'ParseQnum(" < 100 ")': 100,
    'ParseQnum(" > 100 ")': 100,
    'Qualifier(ParseQnum("100"))': '=',
    'Qualifier(ParseQnum("<100"))': '<',
    'Qualifier(ParseQnum(">100"))': '>',
    'QnumToString(ParseQnum(" < 100 "))': '<100',
  },
  ToString: {
    'ToString(1)': '1',
    'ToString(3.14)': '3.14',
    'ToString(true)': 'true',
  },
};

function expectTyped(actual: any, expected: any) {
  if (Array.isArray(expected))
    expectArray(actual, expected);
  else if (typeof expected === 'object' && expected instanceof dayjs)
    expect((<dayjs.Dayjs>expected).isSame(actual), true);
  else if (typeof expected === 'object' && expected !== null)
    expectObject(actual, expected);
  else if (typeof expected === 'number' && Number.isFinite(expected) && !Number.isInteger(expected))
    expectFloat(actual, expected);
  else
    expect(actual == expected, true);
}
