import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {FilterDesc} from '../types';
import {_package} from '../package';
import {DataLoader} from '../utils/data-loader';

export class PtmFilterApp {
  private viewed: boolean = false;

  private readonly dataLoader: DataLoader;

  private dataDf: DG.DataFrame;
  private predictedPtmDf: DG.DataFrame;
  private observedPtmDf: DG.DataFrame;

  private antigenName: string = 'IAPW8';

  constructor(dataLoader: DataLoader) {
    this.dataLoader = dataLoader;
  }

  async init(): Promise<void> {
    await this.loadData();
  }

  async loadData(): Promise<void> {
    const dataDf: DG.DataFrame = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id1', ['id1_0001', 'id1_0002']),
      DG.Column.fromStrings('id2', ['id2_001', 'id2_002'])]);

    const [predictedPtmDf, observedPtmDf]: [DG.DataFrame, DG.DataFrame] = await Promise.all([
      this.dataLoader.getPredictedPtmByAntigen(this.antigenName),
      this.dataLoader.getObservedPtmByAntigen(this.antigenName)]);

    // const predictedPtmDf: DG.DataFrame = DG.DataFrame.fromColumns([]);
    // const observedPtmDf: DG.DataFrame = DG.DataFrame.fromColumns([]);

    await this.setData(dataDf, predictedPtmDf, observedPtmDf);
  }

  async setData(dataDf: DG.DataFrame, predictedPtmDf: DG.DataFrame, observedPtmDf: DG.DataFrame): Promise<void> {
    if (this.viewed) {
      await this.destroyView();
      this.viewed = false;
    }

    this.dataDf = dataDf;
    this.predictedPtmDf = predictedPtmDf;
    this.observedPtmDf = observedPtmDf;

    if (!this.viewed) {
      await this.buildView();
      this.viewed = true;
    }
  }

  // -- View --

  private tableView?: DG.TableView;
  private filterView?: DG.FilterGroup;

  async destroyView(): Promise<void> {
    if (this.tableView) {
      this.tableView.close();
      delete this.tableView;
    }
  }

  async buildView(): Promise<void> {
    if (!this.tableView && this.dataDf) {
      //
      this.tableView = grok.shell.addTableView(this.dataDf, DG.DOCK_TYPE.FILL);
      //this.tableView.basePath = `func/${_package.name}.ptmFilterApp`;
      this.tableView.path = this.tableView.basePath = `func/${_package.name}.ptmFilterApp`;

      const predictedPtmCsv: string = this.predictedPtmDf.toCsv();
      const observedPtmCsv: string = this.observedPtmDf.toCsv();


      const filterList: FilterDesc[] = [
        {column: 'id1', type: DG.FILTER_TYPE.CATEGORICAL, label: 'flt id1'},
        ...PtmFilterApp.buildFilterList(
          'chothia', predictedPtmCsv, observedPtmCsv),
        {column: 'id2', type: DG.FILTER_TYPE.CATEGORICAL, label: 'flt id2'}];

      // const filterList: FilterDesc[] = [
      //   {column: 'v id', type: DG.FILTER_TYPE.CATEGORICAL, label: 'v id'}];

      // const filterList: FilterDesc[] = [];

      this.filterView = this.tableView.filters({filters: filterList}) as DG.FilterGroup;
      this.tableView.dockManager.dock(this.filterView, DG.DOCK_TYPE.LEFT);
    }
  }

  // -- Filters --

  static buildFilterList(cdrName: string, predictedPtmCsv: string, observedPtmCsv: string): FilterDesc[] {
    const res: FilterDesc[] = [];

    // const filterType = `${_package.name}:ptmFilter`;
    const filterType = `${_package.name}:mlbPtmFilter`;
    res.push({
      type: filterType,
      currentCdr: cdrName,
      predictedPtmCsv: predictedPtmCsv,
      observerdPtm: observedPtmCsv,
    });
    return res;
  }
}
