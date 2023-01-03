import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';
import { interval } from 'rxjs';


export class ScriptingTutorial extends Tutorial {
  get name() {
    return 'Scripting';
  }
  get description() {
    return 'Scripting is an integration mechanism with languages for statistical computing';
  }
  get steps() { return 11; }

  demoTable: string = '';
  helpUrl: string = 'https://datagrok.ai/help/compute/scripting';

  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Scripting is an integration mechanism with languages for statistical computing.');

    this.describe(ui.link('More about ' + this.name, this.helpUrl).outerHTML);

    this.title('Run a script');

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

    await this.action('Open a sample table for the script', grok.events.onTableAdded.pipe(
      filter((data) => data.args.dataFrame.name === sampleDfName)),
      $('div.d4-ribbon-item').has('i.grok-icon.fa-asterisk')[0],
      'In front of you is a valid script. The commented out section on top defines script parameters. ' +
      'This simple script calculates the number of cells in a dataframe. The <i class="grok-icon fal fa-asterisk"></i> ' +
      'icon opens a demo table for you. It appears only for scripts annotated with a special <i>sample</i> parameter.');

    const playBtn = $('div.d4-ribbon-item').has('i.grok-icon.fa-play')[0];
    let callEditorDlg = await this.openDialog('Run the script', sampleScriptName, playBtn,
      'Before a script gets executed, all its input parameters should be set. If there are any, a dialog like this one ' +
      'will appear. Here we need provide only one input named "Table".');

    await this.dlgInputAction(callEditorDlg, `Set "Table" to ${sampleDfName}`, 'Table', sampleDfName,
      'As we have added a sample table in the previous step, the input choices include this table. ' +
      'There you will also see the names of other currently opened tables.');

    await this.action('Click "OK"', grok.functions.onAfterRunAction.pipe(
      filter((c) => c.func.name.startsWith(sampleScriptName) && c.getOutputParamValue() === sampleScriptOutput)),
      $(callEditorDlg.root).find('button.ui-btn').filter((idx, btn) => btn.textContent === 'OK')[0]);

    const scriptOutputInfo = 'If a script returns a scalar value, it gets printed to the console. ' +
      'For dataframe outputs, the platform additionally opens a table view. In our case, the result ' +
      'is a numeric value, so let\'s open the console to see it. Go to <b>Windows | Console</b> and ' +
      'turn on the switch. You can use the <b>~</b> key to control the console visibility.';
    grok.shell.windows.showConsole = false;

    await this.action('Find the results in the console',
      interval(1000).pipe(filter(() => grok.shell.windows.showConsole)),
      this.getSidebarHints('Windows', 'Console'),
      scriptOutputInfo);

    // @ts-ignore
    const editor = sv.root.lastChild.lastChild.CodeMirror;
    const doc = editor.getDoc();
    const scriptBodyIndex = doc.getValue().split('\n').findIndex((line: string) => !line.startsWith('#'));
    doc.replaceRange('\n', { line: scriptBodyIndex - 1 });
    const lastLineIndex = doc.lineCount() - 1;
    const newOutputParam = '#output: dataframe clone';
    const newOutputParamDef = 'clone = table';

    const dfCloneTip = 'Scripts are not limited to one output parameter. Let\'s check this by appending ' +
    `<b>${newOutputParam}</b> to the script annotation (line ${scriptBodyIndex + 1}). The cloning part is ` +
    `quite straightforward: put the initial dataframe into a new variable (paste <b>${newOutputParamDef} ` +
    `</b> to line ${lastLineIndex + 1}).`;

    await this.action('Add the second output value to the script', interval(1000).pipe(
      filter(() => doc.getLine(scriptBodyIndex).trim() === newOutputParam &&
        doc.getLine(lastLineIndex).trim() === newOutputParamDef)),
      null, dfCloneTip);

    callEditorDlg = await this.openDialog('Run the script', sampleScriptName, playBtn);
    const historyInfo = 'Find the previously entered parameter in the dialog\'s history.';
    await this.dlgInputAction(callEditorDlg, `Set "Table" to ${sampleDfName}`, 'Table',
      sampleDfName, historyInfo, true);
    await this.action('Click "OK"', grok.functions.onAfterRunAction.pipe(
      filter((c) => c.func.name.startsWith(sampleScriptName) && c.outputs.size() === 2 &&
        c.outputs.get('clone') instanceof DG.DataFrame)),
      $(callEditorDlg.root).find('button.ui-btn').filter((idx, btn) => btn.textContent === 'OK')[0]);
  }
}
