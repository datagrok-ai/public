import * as DG from "datagrok-api/dg";

export class Template {
  public name: string = '';
  public query: string = '';
  public queryParams: object = {};
  public data?: DG.DataFrame;
  public context?: DG.Context;
}