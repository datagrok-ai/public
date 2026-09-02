import * as DG from 'datagrok-api/dg';

// Resolves the ArtifactAlignment publish dialog if the package is installed (no npm dependency)
export function findPublishDialogFunc(): DG.Func | null {
  return DG.Func.find({name: 'publishWorkflowRunDialog'})[0] ?? null;
}

export async function publishLiveRun(func: DG.Func, sourceCall: DG.FuncCall, defaultName?: string): Promise<void> {
  await func.prepare({sourceCall, defaultName}).call();
}

export async function publishSavedRun(func: DG.Func, sourceMetaCallId: string, defaultName?: string): Promise<void> {
  await func.prepare({sourceMetaCallId, defaultName}).call();
}
