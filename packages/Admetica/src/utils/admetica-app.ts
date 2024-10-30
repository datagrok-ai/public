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
  mode: string = 'single';
  container: HTMLElement;
  formContainer: HTMLElement;
  modeContainer: HTMLElement;
  sketcherDiv: HTMLElement;
  sketcher?: HTMLElement; // Keep reference to sketcher as HTMLElement

  hintShownForMode: boolean = false;
  hintShownForSketcher: boolean = false;
  hintShownForForm: boolean = false;
  indicators: HTMLElement[] = [];
  popups: HTMLElement[] = [];
  
  placeholder: HTMLElement;

  contentContainerStyle = {
    flex: '1',
    height: '100%',
    border: '1px solid #ccc',
    borderRadius: '5px',
    boxShadow: '0 2px 5px rgba(0,0,0,0.2)',
    backgroundColor: '#fff',
    overflow: 'hidden',
    padding: '10px',
    display: 'flex',
    flexDirection: 'column',
  };

  constructor(parentCall: DG.FuncCall) {
    this.parentCall = parentCall;

    this.container = ui.divH([], { style: { width: '100%', height: '100%', display: 'flex', gap: '20px' } });
    this.formContainer = ui.box(null, { style: { ...this.contentContainerStyle } });
    this.sketcherDiv = ui.div([], { style: { ...this.contentContainerStyle, border: 'none' } });
    this.modeContainer = ui.box(null, {
      style: {
        flex: '1',
        height: '100%',
        border: '1px solid #ccc',
        borderRadius: '5px',
        boxShadow: '0 2px 5px rgba(0,0,0,0.2)',
        backgroundColor: '#fff',
        display: 'flex',
        flexDirection: 'column',
        justifyContent: 'space-between',
      }
    });
    
    this.placeholder = ui.divText('Please select a molecule or upload a file to see the predictions.', {
      style: {
        textAlign: 'center',
        padding: '40px',
        color: '#333',
        fontSize: '18px',
        fontWeight: 'bold',
        backgroundColor: '#f9f9f9',
        borderRadius: '10px',
        margin: '20px',
        boxShadow: '0 4px 15px rgba(0, 0, 0, 0.1)',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        height: '100%',
      }
    });

    // Add keyframes for fade-in animation
    const styleSheet = document.styleSheets[0];
    styleSheet.insertRule(`
      @keyframes fadeIn {
        from { opacity: 0; }
        to { opacity: 1; }
      }
    `, styleSheet.cssRules.length);

    this.formContainer.appendChild(this.placeholder);
    this.container.appendChild(this.modeContainer);
    this.container.appendChild(this.formContainer);
  }

  async init() {
    //const smiles = 'CC(=O)Oc1ccccc1C(=O)O';
    this.sketcher = grok.chem.sketcher(async (smiles: string) => await this.onChanged(smiles));
    this.sketcherDiv.appendChild(this.sketcher);

    const modeTabs = ui.tabControl();
    modeTabs.onTabChanged.subscribe((tab: DG.TabPane) => {
      this.clearTable();
      this.clearForm();
      this.mode = tab.header.innerText.split(' ')[0].toLowerCase();
      //this.resetSketcher(); // Reset sketcher on mode change
    });

    modeTabs.addPane('Single Molecule', () => {
      this.mode = 'single';
      return this.sketcherDiv;
    });
    
    modeTabs.addPane('File Input', () => {
      this.mode = 'file';
      return this.createFileInputPane();
    });

    this.modeContainer.appendChild(modeTabs.root);

    const table = DG.DataFrame.create(1);
    this.tableView = DG.TableView.create(table, false);
    this.tableView.parentCall = this.parentCall;
    this.tableView.dataFrame.name = 'Table';

    setTimeout(async () => {
      this.tableView!._onAdded();
      await this.refresh(table, this.container);
      (grok.shell.view(DG.VIEW_TYPE.BROWSE) as DG.BrowseView).preview = this.tableView!;
      this.showSketcherHint();
    }, 300);
  }

  clearTable() {
    if (this.tableView)
      this.tableView.dataFrame = DG.DataFrame.create(1);
  }

  clearForm() {
    this.formContainer.innerHTML = '';
    this.formContainer.appendChild(this.placeholder);
  }

  showSketcherHint() {
    if (!this.hintShownForSketcher) {
      const hint = ui.hints.addHintIndicator(this.sketcherDiv, false, 7000);
      hint.style.cssText = `
        background-color: rgba(255, 255, 255, 0.9); /* Light background for better visibility */
      `;
      
      const popup = ui.hints.addHint(this.sketcherDiv, ui.divText('Draw a structure'), ui.hints.POSITION.LEFT);
      this.indicators.push(hint);
      this.popups.push(popup);
    }
  }

  closeHints() {
    this.indicators.forEach(indicator => indicator.click());
    this.popups.forEach(popup => popup.remove());
    this.indicators = [];
    this.popups = [];
  }

  async onChanged(smiles: string) {
    if (this.mode === 'single') {
      this.closeHints();
      const col = this.tableView?.dataFrame.columns.getOrCreate('smiles', 'string', 1);
      col!.semType = DG.SEMTYPE.MOLECULE;
      this.tableView?.dataFrame.set('smiles', 0, smiles);
      await grok.data.detectSemanticTypes(this.tableView!.dataFrame);
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
      
      const form = createDynamicForm(this.tableView!.dataFrame, models.split(','), 'smiles', true);
      this.tableView!.grid.invalidate();
      this.formContainer.innerHTML = ''; 
      this.formContainer.appendChild(form.root);
      this.closePlaceholder(); // Hide placeholder when form is populated
    }
  }

  createFileInputPane() {
    return ui.divV([
      ui.input.file('Choose File', {
        onValueChanged: async (file) => {
          const csvData = await file.readAsString();
          const table = DG.DataFrame.fromCsv(csvData);
          this.tableView!.dataFrame = table;
          await grok.data.detectSemanticTypes(this.tableView!.dataFrame);

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

          const form = createDynamicForm(this.tableView!.dataFrame, models.split(','), table.columns.bySemType(DG.SEMTYPE.MOLECULE)!.name, true);
          this.tableView!.grid.invalidate();
          this.formContainer.innerHTML = ''; 
          this.formContainer.appendChild(form.root);
          this.closePlaceholder(); // Hide placeholder when form is populated
          this.refresh(table, this.container);
        }
      })
    ], { style: { padding: '10px', width: '100%' } });
  }

  closePlaceholder() {
    if (this.placeholder.parentElement) {
      this.formContainer.removeChild(this.placeholder);
    }
  }

  refresh(table: DG.DataFrame, modeContainer: HTMLElement) {
    if (table.rowCount > 0) {
      this.tableView?.dockManager.dock(modeContainer, DG.DOCK_TYPE.TOP, null, '', 0.65);
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