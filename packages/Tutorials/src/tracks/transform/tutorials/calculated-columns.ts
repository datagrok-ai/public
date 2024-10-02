import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial, TutorialPrerequisites } from '@datagrok-libraries/tutorials/src/tutorial';
import { Observable } from 'rxjs';


export class CalculatedColumnsTutorial extends Tutorial {
  get name(): string {
    return 'Calculated Columns';
  }
  get description(): string {
    return 'Learn about calculated columns, how to add them to a dataframe, and how to edit predefined formulas.';
  }
  get steps(): number {
    return 13;
  }

  helpUrl: string = 'https://datagrok.ai/help/transform/add-new-column';
  prerequisites: TutorialPrerequisites = { packages: ['PowerPack'] };
  // To do: adjust the tutorial based on presence/absence of the package

  protected async _run(): Promise<void> {
    grok.shell.windows.showContextPanel = true;
    this.header.textContent = this.name;

    this.describe('Columns based on an evaluated expression are called <b>calculated</b>. Such expressions, or ' +
      'formulas, can include mathematical functions, constants, platform objects properties and functions. Moreover, ' +
      'you can use data from existing columns in your formula. Calculated columns are a powerful way to transform data.');

    const addNCDlg = await this.openAddNCDialog();

    this.describe('The dialog contains column name and type inputs and an expression field. Below it, there is a preview ' +
      'that shows the columns used in a formula and the result. Next to it, there is a column and function search.');

    const columnName = 'Height, m';
    await this.dlgInputAction(addNCDlg, `Name a column "${columnName}"`, '', columnName);

    const simpleFormula = 'Div(170, 100)';
    await this.dlgInputAction(addNCDlg, `Enter the expression "${simpleFormula}"`, '', simpleFormula, '', false, 2);

    await this.action(`Enter the expression "${simpleFormula}"`, new Observable((subscriber: any) => {
      const observer = new MutationObserver((mutationsList, observer) => {
        mutationsList.forEach((m) => {
          //@ts-ignore
          if (m.target.innerText === simpleFormula) {
            subscriber.next(true);
            observer.disconnect();
          }
        });
      });
      observer.observe(addNCDlg!.root.querySelector('.cm-line')!, {childList: true, attributes: true});
    }));

    await this.action('Click "OK"', this.t!.onColumnsAdded.pipe(filter((data) =>
      data.args.columns.some((col: DG.Column) => col.name === columnName && col.meta.formula !== null &&
      col.meta.formula === simpleFormula))), addNCDlg.getButton('OK'), 'Once a valid expression is entered, ' +
      'a new column will appear in the preview. Note that the column type is set automatically to "double". The type is ' +
      'determined based on the function output parameter type. You can change the column type manually, if necessary. For ' +
      'convenience, we\'ll automatically change the number formatting to match the format of the original column. You ' +
      'will see the formatted results once the column is added to the grid.');

    const heightCol = this.t!.getCol('height');
    const heightInMetersCol = this.t!.getCol(columnName);
    heightInMetersCol.meta.format = heightCol.meta.format;

    this.describe('Now you should see a new column added to the grid. Let\'s learn how to edit a formula to replace ' +
      'the constant value with the actual patient height from the corresponding column.');

    let accordion: DG.Accordion;
    await this.action(`Click on the "${columnName}" column header`, grok.events.onAccordionConstructed.pipe(
      filter((acc) => {
        if (acc.context instanceof DG.Column && acc.context?.name === columnName) {
          accordion = acc;
          return true;
        }
        return false;
      })));

    accordion!.getPane('Formula').expanded = true;
    const editDlg = await this.openDialog('Click the "Edit" button under the formula field in the context panel',
      'Edit Column Formula', $(accordion!.root).find('div.d4-pane-formula button.ui-btn').filter((idx, el) =>
      el.textContent?.toLowerCase() === 'edit')[0], 'The <b>Formula</b> pane contains the expression the column ' +
      'is calculated on. You can edit it in the field and apply the changes directly from the context panel, ' +
      'or re-open the dialog by pressing "Edit".');

    const formulaWithColInfo = 'To apply a function to column values, drag the column into the dialog formula ' +
      'field either from the grid or from the column list in the dialog (use search input to find a column in ' +
      'large datasets). These actions will create a column reference. The notation is <b>${COLUMN_NAME}</b>, ' +
      'e.g., <b>Div(${HEIGHT}, 100)</b>.<br>You can also press "$" while editing the formula to open up a ' +
      'column list popup and use arrow keys + "Enter" to select a column. If you are editing the formula from ' +
      'the context panel, or in the column properties dialog, type the column name in this notation manually.<br>' +
      'Note that the column type is not updated automatically during editing.<br>Some mathematical functions, ' +
      'such as <i>Div, Mul</i>, and <i>Pow</i>, have equivalent operators. Check out our wiki to learn more about ' +
      ui.link('operators', 'https://datagrok.ai/help/transform/functions/operators').outerHTML;
    const tolerance = 1e-3;

    await this.action('Edit the formula to use the "HEIGHT" column values and click "OK"',
      grok.functions.onAfterRunAction.pipe(filter((call) => {
        const column = call.outputs.get('result');
        return call.func.name === 'AddNewColumn' && column.name === columnName &&
          Math.abs(column.min - 1.275) < tolerance && Math.abs(column.max - 2.033) < tolerance;
      })), editDlg.inputs.filter((input) => input.caption == '')[2]?.root, formulaWithColInfo);

    await this.action('Change the "HEIGHT" value in the first row to "170"', this.t!.onValuesChanged.pipe(filter(() =>
      this.t!.cell(0, 'HEIGHT').value === 170)), null, 'Now we will examine in which circumstances the values of a ' +
      'calculated column get re-calculated. There is a distinction between column data and metadata changes. As ' +
      'you can see, the value of height in meters in the first row doesn\'t change along with our value change. ' +
      'However, if you want to refresh computations after a value change, you can click the <b>Apply</b> button ' +
      'in the <b>Formula</b> pane of the context panel.');

    const addNCDlgBMI = await this.openAddNCDialog('Add a new column that calculates BMI');
    const columnNameBMI = 'BMI';
    await this.dlgInputAction(addNCDlgBMI, `Name a column "${columnNameBMI}"`, '', columnNameBMI);

    await this.action('Enter the BMI formula and click "OK"', grok.functions.onAfterRunAction.pipe(filter((call) => {
      const column = call.outputs.get('result');
      return call.func.name === 'AddNewColumn' && column.name === columnNameBMI &&
        Math.abs(column.min - 12.891) < tolerance && Math.abs(column.max - 62.932) < tolerance;
      })), addNCDlgBMI.inputs.filter((input) => input.caption == '')[2]?.root, 'The body mass index (BMI) is ' +
      `calculated as mass (kg) divided by height (m) raised to power 2:<br>BMI = weight / height^2<br>Use the "WEIGHT" and "${columnName}" ` +
      'columns and functions "Div" and "Pow" (or the corresponding operators). We will use this new column to ' +
      'check what happens when we change the column metadata.');

    await this.action(`Update the formula for "${columnName}" to round the values to 2 decimal places`,
      grok.functions.onAfterRunAction.pipe(filter((call) => {
        const column = call.outputs.get('result');
        return call.func.name === 'AddNewColumn' && column.name === columnName &&
          Math.abs(column.min - 1.279) < tolerance && Math.abs(column.max - 2.029) < tolerance;
      })), null, 'You can apply the new formula from the <b>Formula</b> pane of the context panel. Use the ' +
      '"RoundFloat" function with two arguments (the previous expression column and the number of decimal places).' + 
      `Enter the new formula and click \'APPLY\' button. Pay attention to the "${columnNameBMI}" column. ` +
      `When we change the formula of the underlying column (that is, its metadata), re-calculation is triggered automatically.`);
    
    this.describe('Calculated columns can be based on various functions: core functions (shown in the function search), ' +
      'platform commands, scripts, and package functions. Aside from core functions, you need to specify a fully-' +
      'qualified function name.');
  }
}
