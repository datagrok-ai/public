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
