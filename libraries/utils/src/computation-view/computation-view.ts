/* eslint-disable max-len */
/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FunctionView} from '../function-view';
import {Component, StateService} from './common/service-interfaces';
import {tokens} from './common/inject-tokens';

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
  public static inject = [
    tokens.stateService,
    tokens.inputPanel,
    tokens.outputPanel,
    tokens.exportPanel,
  ] as const;

  /** Not intented to be used directly in common case. Use {@link ComputationViewBuilder} instead. */
  constructor(
    public state: StateService,
    public inputPanel: Component,
    public outputPanel: Component,
    public exportButton: Component,
  ) {
    super(state.func);

    if (this.state.func != null)
      this.render(this.state.func);
  }

  /** Override to provide custom initialization. {@link onViewInitialized} gets fired after that. */
  public async render(func: DG.Func): Promise<void> {
    this.state.func = func;

    await super.init(this.state.func);
    this.root.innerHTML = '';
    this.root.append(ui.divH([this.inputPanel.root, this.outputPanel.root]));

    this.setRibbonPanels([[
      this.exportButton.root,
    ]]);
  }
}
