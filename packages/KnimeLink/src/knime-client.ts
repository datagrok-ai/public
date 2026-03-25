import {KnimeInputParam, KnimeWorkflowSpec, KnimeExecutionInput, KnimeExecutionResult, KnimeJobStatus, KnimeDeployment, KnimeOutputResource} from './types';

export interface IKnimeClient {
  listDeployments(typeFilter?: string): Promise<KnimeDeployment[]>;
  getWorkflowInputs(id: string): Promise<KnimeWorkflowSpec>;
  executeSyncWorkflow(id: string, input: KnimeExecutionInput): Promise<KnimeExecutionResult>;
  startAsyncJob(id: string, input: KnimeExecutionInput): Promise<string>;
  getJobStatus(jobId: string): Promise<KnimeJobStatus>;
  getJobResult(jobId: string): Promise<KnimeExecutionResult>;
  buildJobResult(jobId: string, data: any): Promise<KnimeExecutionResult>;
  fetchOutputResource(jobId: string, resourceId: string): Promise<KnimeOutputResource>;
  cancelJob(jobId: string): Promise<void>;
  /** Build a proxy URL for the workflow SVG image, suitable for use as img.src. */
  getWorkflowImageUrl(workflowId: string): Promise<string | null>;
}
