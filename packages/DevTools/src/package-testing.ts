import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {list} from "datagrok-api/ui";

class FunctionTest {
  fun: DG.Func;
  status: HTMLDivElement;
}

export async function testPackages() {

  let testFunctions: FunctionTest[] = DG.Func
    .find({tags: [DG.FUNC_TYPES.UNIT_TEST]})
    .map((f) => {
      return {
        fun: f,
        status: ui.divText('Not run')
      }
    });

  function run() {
    for (let f of testFunctions) {
      f.fun.apply()
        .then((s) => { f.status.innerHTML = 'OK'; f.status.style.color = 'green'; })
        .catch((error) => f.status.innerHTML = `${error}`);
    }
  }

  ui.dialog()
    .add(ui.table(testFunctions, (f) => [f.fun.f.fun, f.status]))
    .add(ui.buttonsInput([ui.button('RUN', run) ]))
    .show();
}