/* eslint-disable max-len */
/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as rxjs from 'rxjs';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DataFrame, FuncCall} from 'datagrok-api/dg';
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
  lastCall?: DG.FuncCall;
  _inputFields: Map<string, DG.InputBase> = new Map<string, DG.InputBase>();

  onViewInitialized: rxjs.Subject<void> = new rxjs.Subject();
  onInputChanged: rxjs.Subject<void> = new rxjs.Subject();
  onComputationStarted: rxjs.Subject<FuncCall> = new rxjs.Subject();
  onComputationCompleted: rxjs.Subject<FuncCall> = new rxjs.Subject();
  onComputationSucceeded: rxjs.Subject<FuncCall> = new rxjs.Subject();
  onComputationError: rxjs.Subject<FuncCall> = new rxjs.Subject();

  constructor(func: DG.Func | null) {
    super(func);

    if (func != null)
      this.init(func);
  }

  /** All inputs that are bound to fields */
  get inputFields(): Map<string, DG.InputBase> {return this._inputFields;}

  /** Saves the computation results to the historical results, returns its id. See also {@link loadRun}. */
  async saveRun(call: FuncCall): Promise<string> {return 'xxx'; /* await grok.dapi.functions.calls.save(call);*/}

  /** Loads the specified historical results. See also {@link saveRun}. */
  async loadRun(runId: string): Promise<void> {/* await grok.dapi.functions.calls.find(call.id);*/}

  /** The actual computation function. */
  async compute(call: FuncCall): Promise<void> {await call.call();}

  /** Maps inputs to parameters, computes, and maps output parameters to the UI. */
  async run(): Promise<void> {
    this.lastCall = this.call!.clone();
    this.inputFieldsToParameters(this.lastCall);

    try {
      await this.compute(this.lastCall);
    } catch (e) {
      this.onComputationError.next(this.lastCall);
    }

    this.outputParametersToView(this.lastCall);
  }

  /** Override to provide custom initialization. {@link onViewInitialized} gets fired after that. */
  async init(func: DG.Func): Promise<void> {
    await super.init(func);
    // this.initTopMenu(this.ribbonMenu);
    this.initRibbonPanels();
    return this.onViewInitialized.next();
  }

  /** Override to customize top menu*/
  initTopMenu(menu: DG.Menu) {
    if (this.func) {
      menu
        .group(this.func.friendlyName)
        .item('Run', () => this.run());
    }
  }

  /** Override to customize ribbon panels */
  initRibbonPanels() {
    this.setRibbonPanels([[
      ui.bigButton('RUN', () => this.run()),
      ui.comboPopup(
        ui.icons.save(null, 'Export'),
        this.supportedExportFormats,
        async (format: string) =>
          DG.Utils.download(`${this.name} - ${new Date().toLocaleString()}.${this.supportedExportExtensions[format]}`, await this.export(format)),
      ),
    ]]);
  }

  /** Builds the complete view. Override it if completely custom implementation is needed,
   * otherwise consider overriding {@link buildInputBlock} and/or {@link buildOutputBlock} */
  build(): HTMLElement {return ui.splitH([this.buildInputBlock(), this.buildOutputBlock()]);}

  /** Override to build completely custom input block (the one with inputs, usually on the left).
   * If only specific field, consider overriding {@link buildCustomInputs}. */
  buildInputBlock(): HTMLElement {
    return super.buildInputBlock();
  }

  /** Custom inputs for the specified fields. Inputs for these fields will not be created by {@link buildInputBlock}. */
  buildCustomInputs(): Map<String, DG.InputBase> {return new Map();}

  /** Override to create output block. */
  buildOutputBlock(): HTMLElement {return super.buildOutputBlock();}

  /** Override to create a custom input control. */
  buildRibbonPanel(): HTMLElement {
    return ui.divH([
      ui.button('RUN', () => {}),
      ui.comboPopup('Export',
        this.supportedExportFormats, (format: string) => this.export(format)),
    ]);
  }

  /** Creates function parameters based on the data entered by the user. */
  inputFieldsToParameters(call: FuncCall): void {super.inputFieldsToParameters(call);}

  /** Visualizes computation results */
  outputParametersToView(call: FuncCall): void {super.outputParametersToView(call);}

  /** Override to provide custom computation error handling. */
  processComputationError(call: FuncCall) {grok.shell.error(call.errorMessage!);}

  /** Override to provide supported export formats.
   * These formats are available under the "Export" popup on the ribbon panel. */
  get supportedExportFormats(): string[] {
    return [
      'Excel', 'PDF', 'CSV',
    ];
  }

  /** Override to provide custom file extensions for exported formats.
   * These formats are available under the "Export" popup on the ribbon panel. */
  get supportedExportExtensions(): Record<string, string> {
    return {
      'Excel': 'xlsx',
      'PDF': 'pdf',
      'CSV': 'csv',
    };
  }

  /** Override to provide custom export. */
  async export(format: string): Promise<Blob> {
    if (format == 'CSV')
      return new Blob([(this.lastCall?.getOutputParamValue() as DataFrame).toCsv()]);
    throw new Error(`Format "${format}" is not supported.`);
  }
}
