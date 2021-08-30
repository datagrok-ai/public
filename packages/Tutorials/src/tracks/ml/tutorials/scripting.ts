import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { filter } from 'rxjs/operators';
import { Tutorial } from "../../../tutorial";


export class ScriptingTutorial extends Tutorial {
  get name() { return 'Scripting'; }
  get description() {
    return 'Scripting is an integration mechanism with languages for statistical computing';
  }

  protected async _run() {
    this.header.textContent = 'Create and run a script';
    let sv: DG.View;

    await this.action('Click on "Functions | Scripts | New Script" to open a script editor',
      grok.events.onViewAdded.pipe(filter((v: DG.View) => {
        if (v.type === 'ScriptView') {
          sv = v;
          return true;
        }
        return false;
      })));

    this.describe('This is a script editor. Here, you write code and bind the parameters to ' +
      'the sample dataset (press F1 to get help on parameter format). Also, the editor ' +
      'lets you load previously saved scripts, including the samples designed to help ' +
      'better understand the platform.');

  }
}
