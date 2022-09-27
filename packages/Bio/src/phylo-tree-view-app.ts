import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export class PhyloTreeViewApp {

  private tv!: DG.TableView;

  private ptv!: DG.Viewer; // PhyloTreeViewer
  private pctv!: DG.Viewer; // PhylocanvasGL

  private _df!: DG.DataFrame;

  get df() { return this._df; }

  constructor() {

  }

  async init(): Promise<void> {

    await this.loadData();
  }

  async loadData(): Promise<void> {
    const treeData: string = await grok.dapi.files.readAsText('System:AppData/Bio/data/trees/tree95.nwk');

    const df: DG.DataFrame = (await grok.functions.call(
      'PhyloTreeViewer:_newickToDf', {newick: treeData, filename: 'filename'})) as DG.DataFrame;
    // const col: DG.Column = DG.Column.fromList('string', 'id', []);
    // const df: DG.DataFrame = DG.DataFrame.fromColumns([col]);
    // df.setTag('.newick', treeData);

    await this.setData(df);
  }

  async setData(df: DG.DataFrame): Promise<void> {
    await this.destroyView();

    this._df = df;

    await this.buildView();
  }

  //#region -- View --

  async destroyView(): Promise<void> {

  }

  async buildView(): Promise<void> {
    if (!this.tv) {
      this.tv = grok.shell.addTableView(this.df);
    }

    if (!this.ptv) {
      this.ptv = this.tv.addViewer('PhyloTree');
    }

    if (!this.pctv) {
      this.pctv = this.tv.addViewer('PhylocanvasGlTree');
    }
  }

  //#endregion -- View --
}