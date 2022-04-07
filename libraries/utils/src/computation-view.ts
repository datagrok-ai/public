import * as grok from 'datagrok-api/grok';
import * as rxjs from 'rxjs';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FuncCall, FuncCallParam} from "datagrok-api/dg";
import {FunctionView} from "./function-view";

/**
 * Base class for handling Compute models (see https://github.com/datagrok-ai/public/blob/master/help/compute/compute.md).
 * In most cases, the computation is a simple {@link Func}
 * Extend it in cases where a behavior or UI not supported by the {@link FunctionView} is needed.
 *
 * It provides the following functionality out-of-the-box, where each section could be customized:
 * - a structured way to represent input and output parameters: {@link parameters}
 * - generic way to generate UI for inputs, outputs, and interactivity (running the model, etc)
 *   - persisting historical results to the db (via {@link parameters})
 * - export (to Excel and PDF): {@link saveExcel}, {@link savePdf}
 * - easy loading of historical runs
 * - routing
 * - entering the real, measured (as opposed to predicted) values manually
 * - notifications for changing inputs, completion of computations, etc: {@link onInputChanged}
 * */
export class ComputationView extends FunctionView {

  onInputChanged: rxjs.Observable<void> = new rxjs.Subject();
  onComputationStarted: rxjs.Observable<FuncCall> = new rxjs.Subject();
  onComputationCompleted: rxjs.Observable<FuncCall> = new rxjs.Subject();

  /** List of parameters (both input and output) for this computation */
  get parameters(): FuncCallParam[] { return []; }

  get inputFields(): Iterable<DG.InputBase> { return []; }

  /** Saves the computation results to the historical results, returns its id. See also {@link loadRun}. */
  async saveRun(call: FuncCall): Promise<string> { return 'xxx'; }

  /** Loads the specified historical results. See also {@link saveRun}. */
  async loadRun(runId: string): Promise<void> { }

  /** The actual computation function. */
  async compute(call: FuncCall): Promise<FuncCall> { return new DG.FuncCall(null); }

  async run(): Promise<void> {
    const call = this.inputFieldsToParameters();
    const result = await this.compute(call);
    this.outputParametersToView(result);
  }

  init(): void {}

  /** Builds the complete view. Override it if completely custom implementation is needed,
   * otherwise consider overriding {@link buildInputBlock} and/or {@link buildOutputBlock} */
  build(): HTMLElement { return ui.splitH([this.buildInputBlock(), this.buildOutputBlock()]); }

  /** Override to build completely custom input block (the one with inputs, usually on the left).
   * If only specific field, consider overriding {@link buildInputField}. */
  buildInputBlock(): HTMLElement { return ui.div(); }

  /** Override to create a custom input control. */
  buildInputField(p: DG.Property): DG.InputBase { return ui.intInput('foo', 0); }

  buildOutputBlock(): HTMLElement { return ui.div(); }

  /** Override to create a custom input control. */
  buildRibbonPanel(): HTMLElement { return ui.divH([
    ui.button('RUN', () => {}),
    ui.comboPopup('Export', this.supportedExportFormats, (format) => this.export(format))]);
  }

  inputFieldsToParameters(): FuncCall { return new FuncCall(null); }
  outputParametersToView(call: FuncCall | null): void { }

  /** Override to provide supported export formats.
   * These formats are available under the "Export" popup on the ribbon panel. */
  get supportedExportFormats(): string[] { return ['Excel', 'PDF']; }

  /** Exports  */
  async export(format: string): Promise<Blob> { throw 'Not implemented'; }
}