import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject} from 'rxjs';

export interface Component {
  root: HTMLElement

  render(): void
}

export interface ExportService {
  get supportedExportFormats(): string[]

  /** Override to provide custom file extensions for exported formats.
   * These formats are available under the "Export" popup on the ribbon panel. */
  get supportedExportExtensions(): Record<string, string>

  filename(format: string): string

  /** Override to provide custom export. */
  export(format: string): Promise<Blob>
}

export interface ComputationViewStateService {
  /** Name of view */
  name: string

  /** Function last call info */
  lastCall: DG.FuncCall | null;

  /** Main function of ComputationView. */
  func: DG.Func | null;

  /** The actual computation function. */
  compute(call: DG.FuncCall): Promise<void>;

  /** Maps inputs to parameters, computes, and maps output parameters to the UI. */
  run(): Promise<void>;

  /** Prepare func to run */
  prepare(): void;

  onComputationError: Subject<DG.FuncCall>;
  onComputationCompleted: Subject<DG.FuncCall>;
}

export interface HistoricalRunService {
  /** Saves the computation results to the historical results, returns its id. See also {@link loadRun}. */
  saveRun(call: DG.FuncCall): Promise<string>

  /** Loads the specified historical results. See also {@link saveRun}. */
  loadRun(runId: string): Promise<DG.FuncCall>

  /** Pulls functions' historical runs list */
  pullRuns(): Promise<DG.FuncCall[]>
}
