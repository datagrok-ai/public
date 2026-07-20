/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export async function processChoiceInput(input: DG.ChoiceInput<any>, func: DG.Func, inputProperty: DG.Property) {
  ui.setUpdateIndicator(input.input, true, 'loadeing');

  try {
    if (!inputProperty.choices || inputProperty.choices.length !== 1)
      return;

    const prevInputValue = input.value;
    let choiceString = inputProperty.choices[0].toLowerCase().replaceAll('\\"', '"');
    if (func instanceof DG.DataQuery && func.connection) {
      // property tends to strip the ends...
      if (choiceString.startsWith('uery("') && choiceString.endsWith('"'))
        choiceString = `q${choiceString})`;
      if (choiceString.startsWith('query("') && choiceString.endsWith('")')) {
        const query =
            '--name: choicesQuery\n--output: dataframe out\n' + choiceString.substring(7, choiceString.length - 2);
        const chq = func.connection.query('choicesQuery', query);
        const res = await chq.apply({});
        if (res instanceof DG.DataFrame && res.columns.length > 0 && res.rowCount > 0) {
          input.items = res.columns.byIndex(0).categories;
          if (prevInputValue && input.items.includes(prevInputValue) && input.value !== prevInputValue)
            input.value = prevInputValue;
        }
        return;
      }
    }
  } catch (e) {
    console.error('error loading options for input ' + input.caption, e);
  } finally {
    ui.setUpdateIndicator(input.input, false);
  }
}
