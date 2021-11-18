import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import { _package } from "../package";
import { Logo } from "./ca-viewer-logo";

export class CompostionPviewer {
  root: HTMLElement;

  cdrChoice: DG.InputBase;

  panelNode: DG.DockNode;
  logoNode: DG.DockNode;
  openPanels: DG.DockNode[];

  logo: Logo;

  logoHost: HTMLElement;

  #splitAlignedPeptides(peptideColumn: DG.Column) {
    let splitPeptidesArray: string[][] = [];
    let isFirstRun = true;
    let splitted: string[];

    for (const peptideStr of peptideColumn.toList()) {
      splitted = peptideStr.split('-');

      if (isFirstRun) {
        for (let i = 0; i < splitted.length; i++) {
          splitPeptidesArray.push([]);
        }
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
      if (!isRetained) {
        columnNames.splice(index, 1);
      }
      return isRetained;
    });

    const columnsArray = splitPeptidesArray.map((positionArray, index) => {
      return DG.Column.fromList('string', columnNames[index], positionArray);
    });

    return DG.DataFrame.fromColumns(columnsArray);
  }

  open(mlbView: DG.TableView, mlbTable: DG.DataFrame) {
    // ---- SIDEPANEL REMOVAL ----
    let windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showHelp = false;
    windows.showConsole = false;

    // ---- INPUTS ----
    const cdrOptions = ['clothia'];
    this.cdrChoice = ui.choiceInput('Representation', 'cartoon', cdrOptions);

    // ---- INPUTS PANEL ----
    this.root = ui.div();
    let accOptions = ui.accordion();
    accOptions.addPane('CDR3 Scheme', () => ui.inputs([this.cdrChoice]));
    this.root.append(accOptions.root);

    //Logo view
    let aligned = this.#splitAlignedPeptides(mlbTable.columns.byName('CDR Clothia')!);
    this.logo = new Logo(aligned, mlbTable);
    this.logoHost = ui.divV([this.logo.root]);
    this.logo.render();


    // ---- DOCKING ----
    this.panelNode = mlbView.dockManager.dock(this.root, 'right', null, 'Composition');
    this.logoNode = mlbView.dockManager.dock(this.logoHost, 'left', this.panelNode, 'NGL');

    this.openPanels = [this.panelNode, this.logoNode];
  }

  async close(mlbView: DG.TableView) {
    if (!!this.openPanels)
      this.openPanels.forEach((p) => mlbView.dockManager.close(p));
  }


}