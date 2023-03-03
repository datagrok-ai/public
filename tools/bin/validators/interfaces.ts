export interface ValidationResult {
  message: string,
  value: boolean,
  warnings?: string[],
}

export interface FuncMetadata {
  name: string,
  inputs: FuncParam[],
  outputs: FuncParam[],
  tags?: string[],
  description?: string,
}

export interface FuncParam {
  type: string,
}

export type FuncValidator = ({}: {
  name: string,
  inputs: FuncParam[],
  outputs: FuncParam[],
  tags?: string[],
  description?: string,
}) => ValidationResult;