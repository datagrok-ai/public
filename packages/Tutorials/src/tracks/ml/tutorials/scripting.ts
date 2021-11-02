import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';


export class ScriptingTutorial extends Tutorial {
  get name() {
    return 'Scripting';
  }
  get description() {
    return 'Scripting is an integration mechanism with languages for statistical computing';
  }
  get steps() { return 1; }

  helpUrl: string = 'https://datagrok.ai/help/develop/scripting';

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Scripting is an integration mechanism with languages for statistical computing');

    this.describe(ui.link('More about ' + this.name, this.helpUrl).outerHTML);

    this.title('Create a script');

    const funcPaneHints = this.getSidebarHints('Functions', DG.View.SCRIPTS);
    const editorIntro = 'This is a script editor. Here, you write code and bind the parameters to the ' +
      'sample dataset (press <b>F1</b> to get help on parameter format). Also, the editor lets you load ' +
      'previously saved scripts, including the samples designed to help better understand the platform.';
    const sv = await this.openViewByType(
      'Click on "Functions | Scripts | New Python Script" to open a script editor',
      'ScriptView', funcPaneHints, editorIntro);

    // UI generation delay
    await new Promise((resolve) => setTimeout(resolve, 1000));

    const sampleDfName = 'cars';
    const sampleScriptName = 'Template';
    const sampleScriptOutput = 510;

    await this.action('Open a sample table for the script', grok.events.onViewAdded.pipe(
      filter((v) => v.type === DG.VIEW_TYPE.TABLE_VIEW && (<DG.TableView>v).table?.name === sampleDfName)),
      $('div.d4-ribbon-item').has('i.grok-icon.fa-asterisk')[0],
      'This simple script calculates the number of cells in a dataframe. The <i class="grok-icon fal fa-asterisk"></i> ' +
      'icon opens a demo table for you. It appears only for scripts annotated with a special <i>sample</i> parameter.');

    const callEditorDlg = await this.openDialog('Run a script', sampleScriptName,
      $('div.d4-ribbon-item').has('i.grok-icon.fa-play')[0],
      'Before a script gets executed, all its input parameters should be set. If there are any, a dialog like this one ' +
      'will appear. Here we need provide only one input named "Table".');

    await this.dlgInputAction(callEditorDlg, `Set "Table" to ${sampleDfName}`, 'Table', sampleDfName,
      'As we have added a sample table in the previous step, the input choices include this table. ' +
      'There you will also see the names of other currently opened tables.');

    await this.action('Click "OK"', grok.functions.onAfterRunAction.pipe(
      filter((c) => c.func.name.startsWith(sampleScriptName) && c.getOutputParamValue() === sampleScriptOutput)),
      $(callEditorDlg.root).find('button.ui-btn').filter((idx, btn) => btn.textContent === 'OK')[0]);

  }
}
