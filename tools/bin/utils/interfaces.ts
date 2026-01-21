export interface PackageFile {
  name: string,
  friendlyName?: string,
  fullName?: string,
  description: string,
  version: string,
  author?: {
    name?: string,
    email?: string,
  },
  repository?: {
    type?: string,
    url?: string,
    directory?: string,
  },
  servicePackage?: boolean,
  dependencies?: {
    [dependency: string]: string,
  },
  devDependencies?: {
    [devDependency: string]: string,
  },
  scripts?: {
    [script: string]: string,
  },
  properties?: [
    {
      name: string,
      propertyType: string,
      defaultValue?: any,
      nullable?: boolean,
    }
  ],
  sources?: string[],
  canEdit?: string[],
  canView?: string[],
  category?: string,
}

export interface ValidationResult {
  message: string,
  value: boolean,
  warnings?: string[],
}

export interface Indexable {
  [key: string]: any,
}


export interface FuncMetadata extends Indexable {
  name?: string,
  inputs: FuncParam[],
  outputs: FuncParam[],
  tags?: string[],
  description?: string,
  cache?: string,
  meta?: Record<string, string>,
  invalidateOn?: string,
  isInvalidateOnWithoutCache?: boolean,
  actualType?: string;
}


interface InputOptions {
  semType?: string;
  category?: string;
  optional?: boolean;
  editor?: string;
  nullable?: boolean;
  separators?: string[];
  choices?: string[] | string;
  format?: string;
  min?: string;
  max?: string;
  caption?: string;
  description?: string;
  initialValue?: string;
  viewer?: string;
  units?: string;
  type?: string;
  optionsType?: string;
  step?: string;
  'meta.url'?: boolean;
  metaUrl?: boolean;
}


export interface FuncParam extends InputOptions{
  name?: string, 
  type?: string, 
  actualType?: string, 
  defaultValue?: string, 
  options?: any[],
  optional?: boolean
}

export type FuncValidator = ({}: FuncMetadata) => ValidationResult;
