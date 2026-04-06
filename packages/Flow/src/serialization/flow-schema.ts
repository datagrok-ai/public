/** JSON schema types for the .funcflow.json save format */

export interface FuncFlowDocument {
  version: '1.0';
  name: string;
  description: string;
  author: string;
  created: string;
  modified: string;

  /** LiteGraph's native serialization */
  graph: any;

  /** FuncFlow-specific metadata */
  metadata: FuncFlowMetadata;
}

export interface FuncFlowMetadata {
  nodes: Record<number, FuncFlowNodeMeta>;
  settings: FlowSettings;
}

export interface FuncFlowNodeMeta {
  dgFuncName: string;
  dgNodeType: 'func' | 'input' | 'output' | 'utility';
  paramName?: string;
  defaultValue?: any;
  description?: string;
  paramQualifiers?: Record<string, string>;
}

export interface FlowSettings {
  scriptName: string;
  scriptDescription: string;
  tags: string[];
}
