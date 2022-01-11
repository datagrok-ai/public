import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import { PdbEntry } from './pdb-entry';
import { TwinPviewer } from "./viewers/twin-p-viewer"
import { _package } from "./package";

const EXAMPLE_ID: string = '1BDQ';//'1A2P';//'2v0a';//

export class BioStructureViewer {
  bsView: DG.TableView;
  idInput: DG.InputBase;
  twinPviewer: TwinPviewer;
  ligandSelection: {[key: string]: any};

  private changeId = async (): Promise<void> => {

    //TODO check if ID is present in PDB
    // if (isPresent(this.idInput.value)) {
    //   grok.shell.warning("No PDB data data for associated id");
    //   return;
    // }

    let pi = DG.TaskBarProgressIndicator.create('Creating 3D view');
    const entry = new PdbEntry(this.idInput.value);
    await entry.fetchInfo();

    let pdbStr: string = _package.webRoot + 'pdb/1bdq.pdb';
    entry.sbody = pdbStr;

    if (!this.twinPviewer) {
      this.twinPviewer = new TwinPviewer();
      this.twinPviewer.init(entry, this.bsView, this.ligandSelection);
    } else {
      this.twinPviewer.reset(entry);
    }
    this.twinPviewer.show(this.bsView);

    pi.close();
  };

  private setView = async (): Promise<void> => {

    let table = (await grok.data.loadTable(_package.webRoot + 'src/examples/dock.csv'));

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

  private setIdInput = (id: string): void => {
    this.idInput = ui.stringInput('ID', '');
    this.idInput.value = id;

    this.idInput.root.addEventListener("keyup", (event) => {
      if (event.key === 'Enter')
        this.changeId();
    });
  }

  public async init() {

    await this.setView();
    this.setIdInput(EXAMPLE_ID);
    this.bsView.setRibbonPanels([[this.idInput.root]]);
    this.changeId();
  }
}
