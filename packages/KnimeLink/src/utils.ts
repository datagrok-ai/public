import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {IKnimeClient} from './knime-client';
import {KnimeJobState} from './constants';
import {KnimeExecutionInput, KnimeExecutionResult} from './types';
import {_package} from './package';

export interface PollOptions {
  timeoutMinutes?: number;
  pollingIntervalMs?: number;
  onProgress?: (state: KnimeJobState, message?: string) => void;
}

/** Poll a KNIME async job until it completes, errors, or times out. */
export async function pollJobUntilComplete(
  client: IKnimeClient,
  jobId: string,
  options?: PollOptions,
): Promise<KnimeExecutionResult> {
  const timeoutMs = ((options?.timeoutMinutes ?? getSettingInt('timeoutMinutes', 10)) * 60 * 1000);
  const intervalMs = options?.pollingIntervalMs ?? getSettingInt('pollingIntervalMs', 3000);
  const startTime = Date.now();

  const indicator = DG.TaskBarProgressIndicator.create(`KNIME job ${jobId}`);
  try {
    while (true) {
      const elapsed = Date.now() - startTime;
      if (elapsed > timeoutMs) {
        await client.cancelJob(jobId);
        throw new Error(`KNIME job ${jobId} timed out after ${options?.timeoutMinutes ?? 10} minutes`);
      }

      const status = await client.getJobStatus(jobId);
      options?.onProgress?.(status.state, status.message);
      indicator.update(Math.min(90, (elapsed / timeoutMs) * 100), `KNIME: ${status.state}`);

      if (status.state === KnimeJobState.ExecutionFinished) {
        indicator.update(100, 'KNIME: Complete');
        // Use raw data from the status response to avoid a second GET (job may be discarded)
        if (status._rawData)
          return await client.buildJobResult(jobId, status._rawData);
        return await client.getJobResult(jobId);
      }

      if (status.state === KnimeJobState.ExecutionFailed ||
        status.state === KnimeJobState.ExecutionFailedWithContent ||
        status.state === KnimeJobState.LoadError ||
        status.state === KnimeJobState.NotExecutable)
        throw new Error(`KNIME job ${jobId} failed: ${status.message ?? 'Unknown error'}`);

      if (status.state === KnimeJobState.ExecutionCancelled ||
        status.state === KnimeJobState.Discarded ||
        status.state === KnimeJobState.Vanished)
        throw new Error(`KNIME job ${jobId} was ${status.state.toLowerCase()}`);

      await delay(intervalMs);
    }
  } finally {
    indicator.close();
  }
}

function getSettingInt(name: string, defaultValue: number): number {
  const val = _package.settings?.[name];
  if (val == null)
    return defaultValue;
  const num = parseInt(String(val), 10);
  return isNaN(num) ? defaultValue : num;
}

export function delay(ms: number): Promise<void> {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

export function dataFrameFromObjects(objects: any[]): DG.DataFrame | null {
  if (!objects || objects.length === 0)
    return null;
  return DG.DataFrame.fromObjects(objects) ?? null;
}
