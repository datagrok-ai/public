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
  invalidateOn?: string,
  isInvalidateOnWithoutCache?: boolean,
}

export interface FuncParam {
  name: string,
  type: string,
}

export type FuncValidator = ({}: FuncMetadata) => ValidationResult;
