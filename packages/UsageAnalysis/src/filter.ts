export class UaFilter {
  date?: string = 'this week';
  groups?: string[] = [];
  packages?: string[] = ['all'];
  tags?: string[] = ['any'];
  packagesCategories?: string[] = ['any'];
  projects?: string[] = ['all'];
  // events: string[] = ['all'];
  // isExactly: boolean = false;

  public constructor(init?: Partial<UaFilter>) {
    Object.assign(this, init);
  }
}
