import {KnimeInputParam, KnimeWorkflowSpec, KnimeExecutionInput, KnimeExecutionResult, KnimeJobStatus, KnimeDeployment, KnimeOutputResource} from './types';

export interface IKnimeClient {
  listDeployments(): Promise<KnimeDeployment[]>;
  getWorkflowInputs(id: string): Promise<KnimeWorkflowSpec>;
  executeSyncWorkflow(id: string, input: KnimeExecutionInput): Promise<KnimeExecutionResult>;
  startAsyncJob(id: string, input: KnimeExecutionInput): Promise<string>;
  getJobStatus(jobId: string): Promise<KnimeJobStatus>;
  getJobResult(jobId: string): Promise<KnimeExecutionResult>;
  buildJobResult(jobId: string, data: any): Promise<KnimeExecutionResult>;
  fetchOutputResource(jobId: string, resourceId: string): Promise<KnimeOutputResource>;
  cancelJob(jobId: string): Promise<void>;
}
