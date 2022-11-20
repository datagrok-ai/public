import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import { TwinPviewer } from "./viewers/twin-p-viewer"
import { _package } from "./package";

export class BioStructureViewer {
  bsView: DG.TableView;
  twinPviewer: TwinPviewer;
  ligandSelection: {[key: string]: any};

  private changeId = async (): Promise<void> => {
    const chains = ['A', 'B']; 
    let pi = DG.TaskBarProgressIndicator.create('Creating 3D view');

    let pdbStr: string = await _package.files.readAsText( 'samples/1bdq.pdb');

    if (!this.twinPviewer) {
      this.twinPviewer = new TwinPviewer();
      this.twinPviewer.init(pdbStr, this.bsView, this.ligandSelection, chains);
    } else {
      this.twinPviewer.reset(pdbStr);
    }
    this.twinPviewer.show(this.bsView);

    pi.close();
  };

  private setView = async (): Promise<void> => {

    let table = (await grok.data.loadTable(_package.webRoot + 'files/samples/dock.csv'));

    this.bsView = grok.shell.addTableView(table);

    this.ligandSelection = {};
    let ligands = ["R", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "S", "T", "U", "V", "W"];
    for(let i = 0; i < ligands.length; i++)
      this.ligandSelection[ligands[i]] = [false, 400 + i];
    
    this.bsView.grid.columns.byName("ligand")!.cellType = 'html';
    this.bsView.grid.onCellPrepare((gc) => {
      if (gc.isTableCell && gc.gridColumn.name === "ligand") {
        //debugger;
        gc.style.element = ui.divV([
          ui.boolInput("", false, () => {
            this.ligandSelection[gc.cell.value][0] = !this.ligandSelection[gc.cell.value][0];
            this.twinPviewer.changeLigands(this.bsView, this.ligandSelection);
          })
        ]);
      }
    });

    grok.shell.windows.showProperties = false;
    grok.shell.windows.showHelp = false;
    this.bsView.ribbonMenu.clear();
  }

  public async init() {
    await this.setView();
    this.changeId();
  }
}
