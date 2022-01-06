import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import { PdbEntry } from './pdb-entry';
import { TwinPviewer } from "./viewers/twin-p-viewer"
import { _package } from "./package";

const EXAMPLE_ID: string = '1A2P';//'2v0a';//

export class BioStructureViewer {
  bsView: DG.TableView;
  idInput: DG.InputBase;
  twinPviewer: TwinPviewer;

  private changeId = async (): Promise<void> => {

    //TODO check if ID is present in PDB
    // if (isPresent(this.idInput.value)) {
    //   grok.shell.warning("No PDB data data for associated id");
    //   return;
    // }

    let pi = DG.TaskBarProgressIndicator.create('Creating 3D view');
    const entry = new PdbEntry(this.idInput.value);
    await entry.fetchInfo();

    if (!this.twinPviewer) {
      this.twinPviewer = new TwinPviewer();
      this.twinPviewer.init(entry, this.bsView);
    } else {
      this.twinPviewer.reset(entry);
    }
    this.twinPviewer.show(this.bsView);

    pi.close();
  };

  private setView = (): void => {
    this.bsView = grok.shell.addTableView(DG.DataFrame.fromColumns([
      DG.Column.fromList('int', 'empty', []),
    ]));

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

    this.setView();
    this.setIdInput(EXAMPLE_ID);
    this.bsView.setRibbonPanels([[this.idInput.root]]);
    this.changeId();
  }
}
