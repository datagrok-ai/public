import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import { createDynamicForm, getQueryParams, performChemicalPropertyPredictions } from './admetica-utils';

interface ISplash {
  close: () => void;
  el: HTMLElement;
}

export class AdmeticaViewApp {
  view?: DG.View;
  parentCall: DG.FuncCall;
  tableView?: DG.TableView;
  mode: 'single' | 'file' = 'single';
  container: HTMLElement;
  formContainer: HTMLElement;
  modeContainer: HTMLElement;

  constructor(parentCall: DG.FuncCall) {
    this.parentCall = parentCall;

    this.container = ui.divH([], { style: { width: '100%', height: '100%' } });

    this.formContainer = ui.box(null, { 
      style: { 
        width: '50%',
        height: '100%', 
        border: '1px solid #ccc', 
        borderRadius: '5px', 
        boxShadow: '0 2px 5px rgba(0,0,0,0.2)', 
        padding: '10px', 
        backgroundColor: '#fff'
      } 
    });

    this.modeContainer = ui.box(null, { 
      style: { 
        width: '50%',
        height: '100%', 
        padding: '10px', 
        backgroundColor: '#f9f9f9'
      } 
    });

    this.container.appendChild(this.modeContainer);
    this.container.appendChild(this.formContainer);
  }

  async init() {
    const smiles = 'CC(=O)Oc1ccccc1C(=O)O';
    
    const sketcherDiv = ui.box(null, { 
      style: { 
        margin: '0', 
        height: '100%', 
        border: '1px solid #ccc', 
        borderRadius: '5px', 
        boxShadow: '0 2px 5px rgba(0,0,0,0.2)', 
        padding: '10px', 
        backgroundColor: '#fff'
      } 
    });

    const sketcher = grok.chem.sketcher(async (smiles: string) => await this.onChanged(smiles), smiles);
    sketcherDiv.appendChild(sketcher);

    const modeTabs = ui.tabControl();
    modeTabs.addPane('Single Molecule', () => {
      this.mode = 'single';
      return ui.divV([sketcherDiv]);
    });
    modeTabs.addPane('File Input', () => {
      this.mode = 'file';
      return this.createFileInputPane();
    });

    this.modeContainer.appendChild(modeTabs.root);

    const table = DG.DataFrame.fromCsv(`smiles\n${smiles}`);
    await grok.data.detectSemanticTypes(table);
    this.tableView = DG.TableView.create(table, false);
    this.tableView.parentCall = this.parentCall;

    setTimeout(async () => {
      this.tableView!._onAdded();
      await this.refresh(table, this.container);
      (grok.shell.view(DG.VIEW_TYPE.BROWSE) as DG.BrowseView).preview = this.tableView!;
    }, 300);
  }

  createFileInputPane() {
    return ui.divV([
      ui.input.file('Choose File', {
        onValueChanged: async (file) => {
          const csvData = await file.readAsString();
          const table = DG.DataFrame.fromCsv(csvData);
          await grok.data.detectSemanticTypes(table);
          this.tableView!.dataFrame = table;

          const splashScreen = this.buildSplash(this.formContainer, 'Processing...');
          const models = await getQueryParams();
          await performChemicalPropertyPredictions(
            table.columns.bySemType(DG.SEMTYPE.MOLECULE)!,
            table,
            models,
            undefined,
            true,
            false,
            true
          );
          splashScreen.close();

          const form = createDynamicForm(this.tableView!.dataFrame, models.split(','), 'smiles', false);
          this.formContainer.innerHTML = ''; 
          this.formContainer.appendChild(form.root);
          this.refresh(table, this.container);
        }
      })
    ]);
  }

  async onChanged(smiles: string) {
    if (this.mode === 'single') {
      this.tableView?.dataFrame.set('smiles', 0, smiles);
      const models = await getQueryParams();
      const dataFrame = this.tableView!.dataFrame;

      const splashScreen = this.buildSplash(this.formContainer, 'Processing...');
      await performChemicalPropertyPredictions(
        dataFrame.getCol('smiles'),
        dataFrame,
        models,
        undefined,
        true,
        false,
        true
      );

      splashScreen.close();
      const form = createDynamicForm(this.tableView!.dataFrame, models.split(','), 'smiles', false);
      this.formContainer.innerHTML = ''; 
      this.formContainer.appendChild(form.root);
    }
  }

  refresh(table: DG.DataFrame, modeContainer: HTMLElement) {
    if (table.rowCount > 0) {
      this.tableView?.dockManager.dock(modeContainer, DG.DOCK_TYPE.TOP, null, '', 0.5);
    }
  }

  buildSplash(root: HTMLElement, description: string): ISplash {
    const indicator = ui.loader();
    indicator.style.cssText = 'margin: 0 auto; padding-right: 50px;';

    const panel = ui.divV([
      indicator,
      ui.p(description),
    ], { style: { textAlign: 'center', color: 'white' } });

    const loaderEl = ui.div([panel], 'bsv-modal-background');
    loaderEl.style.cssText = 'display: flex; justify-content: center; align-items: center; color: white; background: rgba(0, 0, 0, 0.7); position: fixed; top: 0; left: 0; width: 100%; height: 100%; z-index: 9999;';

    root.append(loaderEl);

    return new class implements ISplash {
      constructor(public readonly el: HTMLElement) {};
      close(): void { this.el.remove(); }
    }(loaderEl);
  }
}