import {forgeDb} from '../generated/db';
import {TrainingRunRecord} from './model-fields';

export async function recordTrainingRun(record: TrainingRunRecord): Promise<string> {
  const [{id}] = await forgeDb.trainingRuns.insert(record);
  return id;
}

export async function linkTrainingRun(runId: string, modelId: string): Promise<void> {
  await forgeDb.trainingRuns.update(runId, {model_id: modelId});
}
