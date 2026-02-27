import {JSONSchema7} from 'ai';
import * as DG from 'datagrok-api/dg';

export interface FuncParam {
  type: string;
  description: string;
  default?: any;
  optional?: boolean;
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

export interface PackageSelection {
  selected_packages: {
    name: string;
    reason: string;
  }[];
}

export type JsonSchema = JSONSchema7;

/* Schemas defined to enable structured output handling */
export const PackageSelectionSchema: JsonSchema = {
  type: 'object',
  properties: {
    selected_packages: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          name: {
            type: 'string',
            description: 'The name of the package to use.',
          },
          reason: {
            type: 'string',
            description: 'The reason why this package is relevant to the user\'s goal.',
          }
        },
        required: ['name', 'reason'],
        additionalProperties: false,
      }
    }
  },
  required: ['selected_packages'],
  additionalProperties: false,
};

export const PlanSchema: JsonSchema = {
  type: 'object',
  properties: {
    goal: {
      type: 'string',
      description: 'The user goal restated clearly in your own words.'
    },
    analysis: {
      type: 'array',
      description: 'Reasoning explaining why the selected steps will achieve the goal.',
      items: {
        type: 'string',
      }
    },
    steps: {
      type: 'array',
      description: 'An ordered list of steps to perform, using available functions.',
      items: {
        type: 'object',
        properties: {
          action: {
            type: 'string',
            description: 'The type of action to perform in this step. Currently only "call_function" is supported.',
            enum: ['call_function']
          },
          function: {
            type: 'string',
            description: 'The exact name of the function to call in this step.'
          },
          inputs: {
            type: 'object',
            description: 'Function inputs as a dictionary of key-value pairs.',
            properties: {},
            required: [],
            patternProperties: {
              '.*': {},
            },
            additionalProperties: false,
          },
          outputs: {
            type: 'array',
            description: 'List of names for the outputs this step will produce.',
            items: {
              type: 'string',
            }
          }
        },
        additionalProperties: false,
        required: ['action', 'function', 'inputs', 'outputs']
      }
    }
  },
  additionalProperties: false,
  required: ['goal', 'analysis', 'steps'],
};

export enum LLMApiNames {
  OpenAIChatCompletions = 'OpenAI Chat Completions',
  OpenAIResponses = 'OpenAI Responses',
  AntropicMessages = 'Antropic Messages'
}

export type LLMApiName = `${LLMApiNames}`;
// export type LLMApiName = 'OpenAI Chat Completions' | 'OpenAI Responses' | 'Antropic Messages'

export type PackageSettings = {
  vectorStoreId?: string;
  defaultFastModel?: string;
  defaultDeepResearchModel?: string;
  defaultCodingModel?: string;
  APIName: LLMApiName;
  apiVersion?: string;
}
