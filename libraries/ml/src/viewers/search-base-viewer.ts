import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

export class SearchBaseViewer extends DG.JsViewer {
  name: string = '';
  semType: string = '';
  limit: number;
  targetColumn?: DG.Column<string>;
  targetColumnName: string;
  initialized: boolean = false;

  constructor(name: string, semType: string) {
    super();
    this.limit = this.int('limit', 10);
    this.targetColumnName = this.string('targetColumnName');
    this.name = name;
    this.semType = semType;
  }

  init(): void {
    this.initialized = true;
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  async onTableAttached(): Promise<void> {
    this.init();

    if (this.dataFrame) {
      this.subs.push(DG.debounce(this.dataFrame.onRowsRemoved, 50)
        .subscribe((_: any) => this.render(true)));
      const compute = this.name !== 'diversity';
      this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50)
        .subscribe((_: any) => this.render(compute)));
      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50)
        .subscribe((_: any) => this.render(false)));
      this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50)
        .subscribe((_: any) => this.render(false)));
      this.targetColumn = this.dataFrame.columns.bySemType(this.semType) as DG.Column<string>;
      this.targetColumnName = this.targetColumn?.name;
      this.getProperty('limit')!.fromOptions({min: 1, max: this.dataFrame.rowCount});
    }
    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    if (!this.initialized)
      return;
    if (property.name === 'targetColumnName') {
      const col = this.dataFrame.col(property.get(this))!;
      if (col.semType === this.semType)
        this.targetColumn = col;
    }
    this.render();
  }

  /** For tests */ public computeRequested: boolean = false;
  public renderPromise: Promise<void> = Promise.resolve();

  protected render(computeData = true): void {
    this.renderPromise = this.renderPromise.then(async () => {
      this.computeRequested = this.computeRequested || computeData;
      await this.renderInt(computeData);
    });
  }

  async renderInt(_computeData: boolean): Promise<void> {

  }

  beforeRender() {
    if (!this.initialized)
      return false;
    if (this.dataFrame && this.targetColumnName &&
          this.dataFrame.col(this.targetColumnName)!.semType !== this.semType) {
      grok.shell.error(`${this.targetColumnName} is not ${this.semType} type`);
      return false;
    }
    return true;
  }
}