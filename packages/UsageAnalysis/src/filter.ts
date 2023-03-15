export class UaFilter {
  date: String = 'this week';
  users: String[] = ['all'];
  packages: String[] = ['all'];
  events: String[] = ['all'];
  isExactly: boolean = false;

  public constructor(init?:Partial<UaFilter>) {
    Object.assign(this, init);
  }
}
