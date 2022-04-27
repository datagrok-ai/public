/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Component, StateService} from '../common/service-interfaces';
import {tokens} from '../common/inject-tokens';

export class DefaultInputPanelComponent implements Component {
  public root: HTMLElement = ui.div();

  /** All inputs that are bound to fields */
  public inputFields: Record<string, DG.InputBase> = {};

  public static inject = [
    tokens.stateService,
    tokens.historyPanel,
  ] as const;

  constructor(
    private state: StateService,
    private loadRunsButton: Component,
  ) {
    this.render();
  }

  render() {
    this.root = ui.div([this.renderRunSection(this.state.lastCall!)], 'ui-div');
  }

  renderRunSection(call: DG.FuncCall): HTMLElement {
    return ui.wait(async () => {
      const runButton = ui.bigButton('Run', async () => {
        await this.state.run();
      });
      const editor = ui.div([], 'ui-form');
      const inputs: DG.InputBase[] = await call.buildEditor(editor, {condensed: true});
      editor.appendChild(ui.divH([ui.divV([this.loadRunsButton.root], {style: {'justify-content': 'center'}}), runButton], {style: {'justify-content': 'space-between'}}));
      return ui.panel([editor]);
    });
  }
}
