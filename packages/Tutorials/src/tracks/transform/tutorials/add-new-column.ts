import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { interval } from 'rxjs';
import { filter } from 'rxjs/operators';
import { Tutorial, TutorialPrerequisites } from '@datagrok-libraries/tutorials/src/tutorial';


export class CalculatedColumnsTutorial extends Tutorial {
  get name(): string {
    return 'Calculated Columns';
  }
  get description(): string {
    return 'Learn about calculated columns, how to add them to a dataframe, and how to edit predefined formulas.';
  }
  get steps(): number {
    return 10;
  }

  helpUrl: string = 'https://datagrok.ai/help/transform/add-new-column';
  prerequisites: TutorialPrerequisites = { packages: ['PowerPack'] };
  // To do: adjust the tutorial based on presence/absence of the package

  protected async _run(): Promise<void> {
    grok.shell.windows.showContextPanel = true;
    this.header.textContent = this.name;
    // To do: remove once https://github.com/datagrok-ai/public/issues/1840 is fixed
    const heightCol = this.t!.getCol('height');
    this.t!.rows.removeWhereIdx((i) => heightCol.isNone(i));

    this.describe('Columns based on an evaluated expression are called <b>calculated</b>. Such expressions, or ' +
      'formulas, can include mathematical functions, constants, platform objects properties and functions. Moreover, ' +
      'you can use data from existing columns in your formula. Calculated columns are a powerful way to transform data.');

    const addNCDlg = await this.openAddNCDialog();

    this.describe('The dialog contains column name and type inputs and an expression field. Below it, there is a preview ' +
      'that shows the columns used in a formula and the result. Next to it, there is a column and function search.');

    const columnName = 'Rounded height';
    await this.dlgInputAction(addNCDlg, `Name a column "${columnName}"`, '', columnName);

    const simpleFormula = 'Round(170.5)';
    await this.dlgInputAction(addNCDlg, `Enter the expression "${simpleFormula}"`, '', simpleFormula, '', false, 2);

    await this.action('Click "OK"', this.t!.onColumnsAdded.pipe(filter((data) =>
      data.args.columns.some((col: DG.Column) => col.name === columnName && col.tags.has(DG.TAGS.FORMULA) &&
      col.tags[DG.TAGS.FORMULA] === simpleFormula))), addNCDlg.getButton('OK'), 'Once a valid expression is entered, ' +
      'a new column will appear in the preview. Note that the column type is set automatically to "int". The type is ' +
      'determined based on the function output parameter type. You can change the column type manually, if necessary.');

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

    const formulaWithCol = 'Round(${HEIGHT})';
    await this.dlgInputAction(editDlg, 'Edit the formula to use the "HEIGHT" column values', '', formulaWithCol,
      '<p>To apply a function to column values, drag the column into the dialog formula field either from the grid ' +
      'or from the column list in the dialog (use search input to find a column in large datasets). These ' +
      'actions will create a column reference. The notation is <b>${COLUMN_NAME}</b>.</p><p>You can also press "$" ' +
      'while editing the formula to open up a column list popup and use arrow keys + "Enter" to select a column. ' +
      'If you are editing the formula from the context panel, or in the column properties dialog, type the column ' +
      'name in this notation manually.</p><p>Note that the column type is not updated automatically during editing.</p>',
      false, 2);

    await this.action('Click "OK"', this.t!.onMetadataChanged.pipe(filter((data) =>
      (data.args.key as unknown as string) === 'formula' && (data.args.change as unknown as string) === 'set' &&
      (data.args.value as unknown as string).toLowerCase() === formulaWithCol.toLowerCase())), editDlg.getButton('OK'));
  }
}
