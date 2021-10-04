import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';


export class ScriptingTutorial extends Tutorial {
  get name() {
    return 'Scripting';
  }
  get description() {
    return 'Scripting is an integration mechanism with languages for statistical computing'+
    '<a href="https://datagrok.ai/help/develop/scripting" target="_blank" class="ui-link d4-link-external">Read more about scripting.</a>';
  }
  get steps() { return 1; }

  protected async _run() {
    this.header.textContent = 'Create and run a script';
    this.describe('Scripting is an integration mechanism with languages for statistical computing');

    const editorIntro = 'This is a script editor. Here, you write code and bind the parameters to the ' +
      'sample dataset (press F1 to get help on parameter format). Also, the editor lets you load ' +
      'previously saved scripts, including the samples designed to help better understand the platform.';
    const sv = await this.openViewByType(
      'Click on "Functions | Scripts | New Script" to open a script editor',
      'ScriptView', null, editorIntro
    );

  }
}
