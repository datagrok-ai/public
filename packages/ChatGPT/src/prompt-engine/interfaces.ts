import * as DG from 'datagrok-api/dg';

export interface FuncParam {
  type: string;
  description: string;
  default?: any;
}

export type ChatGptFuncParams = Record<string, FuncParam>;

export interface FunctionInput {
  type: string;
  properties: ChatGptFuncParams;
}

export interface FunctionOutput {
  name: string;
  type: string;
}

export interface FunctionMeta {
  name: string;
  description: string;
  inputs: FunctionInput;
  outputs: FunctionOutput[];
}

export interface Step {
  action: 'call_function';
  function: string;
  inputs: Record<string, any>;
  outputs: string[];
}

export interface Plan {
  goal: string;
  analysis: string[];
  steps: Step[];
}

export interface PackageInfo {
  name: string;
  description: string;
}

export interface StepExecutionOutput<T = any> {
  value: T;
  meta?: DG.Property;
}

export type ExecutionContext = Record<string, StepExecutionOutput>;

export interface ExecutePlanResult {
  context: ExecutionContext;
  finalResult: StepExecutionOutput | null;
}
