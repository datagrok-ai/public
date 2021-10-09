import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {list} from "datagrok-api/ui";

class FunctionTest {
  fun: DG.Func;
  status: HTMLDivElement;
}

export async function testPackages() {

  let scripts = await grok.dapi.scripts.list();
  scripts = scripts.filter((s) => s.inputs.length == 0);

  let functions = DG.Func
    .find({tags: [DG.FUNC_TYPES.UNIT_TEST]})
    .concat(scripts);

  let testFunctions: FunctionTest[] = functions
    .map((f) => {
      return {
        fun: f,
        status: ui.divText('Not run')
      }
    });

  function runFunction(f: FunctionTest) {
    function set(status: string, color: string) {
      f.status.innerHTML = status;
      f.status.style.color = color;
    }

    set('Running...', 'orange')
    f.fun.apply()
      .then((_) => set('OK','green'))
      .catch((error) => set(`${error}`, 'red'));
  }

  function run() {
    for (let f of testFunctions)
      runFunction(f);
  }

  ui.dialog()
    .add(ui.table(testFunctions, (f) => [
      ui.icons.play(() => runFunction(f)),
      f.fun.package.name,
      f.fun,
      f.status]))
    .add(ui.buttonsInput([ui.button('RUN', run) ]))
    .show();
}