import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Logo} from './ca-viewer-logo';
import {MlbEvents} from '../const';
import {Subscription} from 'rxjs';


export class CompostionPviewer {
  root: HTMLElement;

  cdrChoice: DG.InputBase;

  panelNode: DG.DockNode;
  logoNode: DG.DockNode;
  barCharNode: DG.DockNode;
  openPanels: DG.DockNode[];

  logo: Logo;
  //barChart: StackedBarChart;

  logoHost: HTMLElement;
  barChartHost: HTMLElement;

  aligned: DG.DataFrame;

  isOpen: boolean;

  subs: Subscription[];

  #splitAlignedPeptides(peptideColumn: DG.Column) {
    let splitPeptidesArray: string[][] = [];
    let isFirstRun = true;
    let splitted: string[];

    for (const peptideStr of peptideColumn.toList()) {
      splitted = peptideStr.split('-');

      if (isFirstRun) {
        for (let i = 0; i < splitted.length; i++)
          splitPeptidesArray.push([]);

        isFirstRun = false;
      }

      splitted.forEach((value, index) => {
        splitPeptidesArray[index].push(value === '' ? '-' : value);
      });
    }

    //create column names list
    let columnNames = ['N'];
    columnNames = columnNames.concat(splitPeptidesArray.map((_, index) => `${index + 1 < 10 ? 0 : ''}${index + 1}`));
    columnNames.push('C');

    // filter out the columns with the same values
    splitPeptidesArray = splitPeptidesArray.filter((positionArray, index) => {
      const isRetained = new Set(positionArray).size > 1;
      if (!isRetained)
        columnNames.splice(index, 1);

      return isRetained;
    });

    const columnsArray = splitPeptidesArray.map((positionArray, index) => {
      return DG.Column.fromList('string', columnNames[index], positionArray);
    });

    return DG.DataFrame.fromColumns(columnsArray);
  }

  open(mlbView: DG.TableView, mlbTable: DG.DataFrame) {
    // ---- SIDEPANEL REMOVAL ----
    const windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showHelp = false;
    windows.showConsole = false;

    // ---- INPUTS ----
    const cdrOptions = ['clothia'];
    this.cdrChoice = ui.choiceInput('Representation', 'cartoon', cdrOptions);

    // ---- INPUTS PANEL ----
    this.root = ui.div();
    const accOptions = ui.accordion();
    accOptions.addPane('CDR3 Scheme', () => ui.inputs([this.cdrChoice]));
    this.root.append(accOptions.root);

    this.aligned = this.#splitAlignedPeptides(mlbTable.columns.byName('CDR Clothia')!);
    //Composition analysis Logo view
    this.logo = new Logo(this.aligned, mlbTable);
    this.logoHost = ui.divV([this.logo.root]);
    this.logo.render(this.aligned);

    //Composition analysis barchart view
    // this.barChart = new StackedBarChart(this.aligned, mlbTable);
    // this.barChartHost = ui.divV([this.logo.root]);
    // this.barChart.render();

    // ---- DOCKING ----
    this.panelNode = mlbView.dockManager.dock(this.root, 'right', null, 'Composition');
    this.logoNode = mlbView.dockManager.dock(this.logoHost, 'left', this.panelNode, 'caLogo');
    //this.barCharNode = mlbView.dockManager.dock(this.barChartHost, 'down', this.logoNode, 'caBarChart');

    this.openPanels = [this.panelNode, this.logoNode];

    this.subs.push(grok.events.onCustomEvent(MlbEvents.CdrChanged).subscribe((value: string) => {
      const key = cdrOptions.find((v) => v.toUpperCase() == value.toUpperCase());
      console.debug(`MLB: CompositionPviewer.onCustomEvent(${MlbEvents.CdrChanged}) ` +
        `value ="${value}" -> key="${key}".`);
      this.cdrChoice.value = key ? key : 'default';
    }));

    this.isOpen = true;
  }

  async close(mlbView: DG.TableView) {
    if (!!this.openPanels)
      this.openPanels.forEach((p) => mlbView.dockManager.close(p));
    this.isOpen = false;

    this.subs.forEach((sub) => sub.unsubscribe());
  }

  do(view: DG.TableView) {
    if (this.isOpen) {
      const a = view.grid.dataFrame.clone(view.grid.dataFrame.filter);
      const b = a.columns.byName('CDR Clothia')!;

      this.aligned = this.#splitAlignedPeptides(b);
      this.logo.render(this.aligned);
    }
  }
}
