/* eslint-disable max-len */
/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FunctionView} from './function-view';

/**
 * Base class for handling Compute models (see https://github.com/datagrok-ai/public/blob/master/help/compute/compute.md).
 * In most cases, the computation is a simple {@link Func}
 * Extend it in cases where a behavior or UI not supported by the {@link FunctionView} is needed.
 *
 * It provides the following functionality out-of-the-box, where each section could be customized:
 * - a structured way to represent input and output parameters: {@link parameters}
 * - generic way to generate UI for inputs, outputs, and interactivity (running the model, etc)
 *   - persisting historical results to the db (via {@link parameters})
 * - export (to Excel and PDF): {@link export}
 * - easy loading of historical runs
 * - routing
 * - entering the real, measured (as opposed to predicted) values manually
 * - notifications for changing inputs, completion of computations, etc: {@link onInputChanged}
 * */
export class ComputationView extends FunctionView {
  constructor(func: DG.Func) {
    super(func);
  }

  /** Override to create output block. */
  override buildOutputPanel(): HTMLElement {
    console.log('ouptut updated');
    const outputs = this.lastCall?.outputs as Record<string, any>;
    const panel = ui.accordion('Output data');
    panel.addPane('Output data', () => {
      return ui.divV(
        Object.entries(outputs).map(([key, val]) => ui.span([`${key}: `, `${val}`])),
      );
    });
    return panel.root;
  }

  /** helper methods */
  override buildInputForm(call: DG.FuncCall): HTMLElement {
    return ui.wait(async () => {
      const runButton = ui.bigButton('Run', async () => {
        await this.run();
      });
      const editor = ui.div([], 'ui-form');
      const inputs: DG.InputBase[] = await call.buildEditor(editor, {condensed: true});
      editor.appendChild(ui.divH([ui.divV([this.historyPanel], {style: {'justify-content': 'center'}}), runButton], {style: {'justify-content': 'space-between'}}));
      return ui.panel([editor]);
    });
  };
}
