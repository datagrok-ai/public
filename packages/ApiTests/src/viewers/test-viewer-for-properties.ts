import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export class TestViewerForProperties extends DG.JsViewer {

  public field = 'field value';
  public testPropertyString: string;
  public testPropertyInt: number;

  public onTableAttachedCounter: number = 0;

  constructor() {
    super();
    this.testPropertyString = this.string('testPropertyString', '');
    this.testPropertyInt = this.int('testPropertyInt', -1);
  }

  public onTableAttached() {
    super.onTableAttached();

    this.onTableAttachedCounter++;
  }
}