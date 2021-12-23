import * as DG from "datagrok-api/dg";

export class UaFilter {
    date: String = 'today';
    users: String[] = ['all'];
    events: String[] = ['all'];
    isExactly: boolean = false;

    public constructor(init?:Partial<UaFilter>) {
        Object.assign(this, init);
    }
}