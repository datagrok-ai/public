import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import { addColorCoding, addSparklines, createDynamicForm, getQueryParams, performChemicalPropertyPredictions, properties, setProperties } from './admetica-utils';
import { UiUtils } from '@datagrok-libraries/compute-utils/shared-components';
import '../css/admetica.css';

interface ISplash {
  close: () => void;
  el: HTMLElement;
}

export class AdmeticaViewApp {
  view?: DG.View;
  parentCall: DG.FuncCall;
  tableView?: DG.TableView;
  mode: string = 'sketch';
  container: HTMLElement;
  formContainer: HTMLElement;
  modeContainer: HTMLElement;
  sketcherDiv: HTMLElement;
  sketcher?: HTMLElement;
  sketcherInstance: grok.chem.Sketcher;

  constructor(parentCall: DG.FuncCall) {
    this.parentCall = parentCall;
    this.container = ui.divH([], { classes: 'app-container' });
    this.formContainer = ui.box(null, { classes: 'content-container', style: { height: '100%' } });
    this.sketcherDiv = ui.div([], { classes: 'content-container', style: { border: 'none' } });
    this.modeContainer = ui.box(null, { classes: 'mode-container', style: { height: '100%' } });
    this.sketcherInstance = new grok.chem.Sketcher();

    const styleSheet = document.styleSheets[0];
    styleSheet.insertRule(`
      @keyframes fadeIn {
        from { opacity: 0; }
        to { opacity: 1; }
      }
    `, styleSheet.cssRules.length);

    this.container.appendChild(this.modeContainer);
    this.container.appendChild(this.formContainer);
  }

  async init(): Promise<void> {
    this.sketcherInstance.setMolecule('CC(Oc1ccccc1C(O)=O)=O');
    this.sketcherInstance.onChanged.subscribe(async () => {
      const smiles: string = this.sketcherInstance.getSmiles();
      await this.onChanged(smiles);
    });
    this.sketcherInstance.molInput.value = 'aspirin';
    this.sketcher = this.sketcherInstance.root;
    this.sketcherDiv.appendChild(this.sketcher);
  
    const modeTabs: DG.TabControl = ui.tabControl();
    modeTabs.onTabChanged.subscribe((tab: DG.TabPane) => this.handleModeChange(tab));
    modeTabs.addPane('Sketch', () => this.createSketchPane());
    modeTabs.addPane('File', () => this.fileInputPane());
    this.modeContainer.appendChild(modeTabs.root);

    this.prepareTableView();
  }
  
  private async handleModeChange(tab: DG.TabPane): Promise<void> {
    this.mode = tab.header.innerText.split(' ')[0].toLowerCase();
    this.clearTable();
    this.clearForm();
  
    if (this.mode === 'sketch')
      this.activateSketchMode();
    else
      this.activateFileMode();
  }
  
  private activateSketchMode(): void {
    this.formContainer.style.visibility = 'visible';
    this.formContainer.style.position = 'relative';
    
    if (this.tableView?.grid) {
      this.tableView.grid.root.style.visibility = 'hidden';
      this.tableView.dockManager.dock(this.container, DG.DOCK_TYPE.TOP, null, '', 0.99);
    }

    this.sketcherInstance.onChanged.next();
  }
  
  private createSketchPane(): HTMLElement {
    this.mode = 'sketch';
    this.formContainer.style.visibility = 'visible';
    this.formContainer.style.position = 'relative';
    
    if (this.tableView?.grid) {
      this.tableView.grid.root.style.visibility = 'hidden';
      this.tableView.dockManager.dock(this.container, DG.DOCK_TYPE.TOP, null, '', 0.99);
    }
  
    return this.sketcherDiv;
  }
  
  private fileInputPane(): HTMLElement {
    this.mode = 'file';
    this.formContainer.style.visibility = 'hidden';
    this.formContainer.style.position = 'absolute';
    return this.createFileInputPane();
  }  

  prepareTableView() {
    const table = DG.DataFrame.create(1);
    this.tableView = DG.TableView.create(table, false);
    this.tableView.parentCall = this.parentCall;
    this.tableView.dataFrame.name = 'Table';
    
    setTimeout(async () => {
      this.tableView!._onAdded();
      this.tableView!.grid.root.style.visibility = 'hidden';
      await this.refresh(table, this.container, 0.99);
      (grok.shell.view(DG.VIEW_TYPE.BROWSE) as DG.BrowseView).preview = this.tableView!;
    }, 300);
  }

  clearTable() {
    if (this.tableView) this.tableView.dataFrame = DG.DataFrame.create(1);
  }

  clearForm() {
    this.formContainer.innerHTML = '';
  }

  async onChanged(smiles: string) {
    if (!smiles) {
      this.clearForm();
      return;
    }

    this.clearTable();
    const col = this.tableView?.dataFrame.columns.getOrCreate('smiles', 'string', 1);
    col!.semType = DG.SEMTYPE.MOLECULE;
    this.tableView?.dataFrame.set('smiles', 0, smiles);
    await grok.data.detectSemanticTypes(this.tableView!.dataFrame);

    const models = await getQueryParams();
    const splashScreen = this.buildSplash(this.formContainer, 'Calculating...');

    await performChemicalPropertyPredictions(
      this.tableView!.dataFrame.getCol('smiles'),
      this.tableView!.dataFrame,
      models,
      undefined,
      false,
      false,
      true
    );

    const molIdx = this.tableView?.dataFrame.columns.names().indexOf('smiles');
    await addSparklines(this.tableView!.dataFrame, models.split(','), molIdx! + 1);
    splashScreen.close();
    
    const form = createDynamicForm(this.tableView!.dataFrame, models.split(','), 'smiles', true);
    this.clearForm();
    this.formContainer.appendChild(form.root);
    this.tableView?.grid.invalidate();  
  }

  createFileInputPane() {
    const fileInputEditor = this.initializeFileInputEditor();
    this.removeLabels(fileInputEditor.root);
    this.styleInputEditor(fileInputEditor.root);
    this.setupDragAndDrop(fileInputEditor.root);
    this.removeOptionsIcon(fileInputEditor.root);
    fileInputEditor.root.classList.add('file-input');
    return ui.divV([fileInputEditor], { classes: 'file-input-container' });
  }
  
  private initializeFileInputEditor() {
    const fileInputEditor = UiUtils.fileInput('', null, async (file: File) => {
      await this.processFile(file);
    }, null);
    return fileInputEditor;
  }
  
  private removeLabels(root: HTMLElement): void {
    const labelSelectors = '.ui-label, .ui-input-label';
    const labels = root.querySelectorAll<HTMLElement>(labelSelectors);
    labels.forEach(label => label.remove());
  }
  
  private styleInputEditor(root: HTMLElement): void {
    const inputEditor = root.querySelector<HTMLElement>('.ui-input-editor');
    if (!inputEditor) return;
  
    Object.assign(inputEditor.style, {
      width: '100%',
      height: '100%',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      backgroundColor: '#ffffff',
      color: '#007bff',
      fontSize: '22px',
      cursor: 'pointer',
      textAlign: 'center',
    });
  }
  
  private setupDragAndDrop(root: HTMLElement): void {
    const inputEditor = root.querySelector<HTMLElement>('.ui-input-editor');
    if (!inputEditor) return;
  
    const highlightColor = '#e0f7fa';
    const defaultColor = '#ffffff';
  
    const changeBackgroundColor = (color: string) => {
      inputEditor.style.backgroundColor = color;
    };
  
    inputEditor.addEventListener('dragenter', () => changeBackgroundColor(highlightColor));
    inputEditor.addEventListener('dragleave', () => changeBackgroundColor(defaultColor));
    inputEditor.addEventListener('dragover', (event) => {
      event.preventDefault();
      changeBackgroundColor(highlightColor);
    });
    inputEditor.addEventListener('drop', () => changeBackgroundColor(defaultColor));
  }
  
  private removeOptionsIcon(root: HTMLElement): void {
    const optionsIcon = root.querySelector<HTMLElement>('.ui-input-options .grok-icon.fal.fa-cloud-upload');
    optionsIcon?.remove();
  }  

  private async activateFileMode(): Promise<void> {
    this.formContainer.style.visibility = 'hidden';
    this.formContainer.style.position = 'absolute';
  
    if (this.tableView?.grid) {
      this.tableView.grid.root.style.visibility = 'visible';
      this.tableView.dockManager.dock(this.container, DG.DOCK_TYPE.TOP, null, '', 0.2);
  
      const df = await openMoleculeDataset('System:AppData/Admetica/demo_files/sar-small_app.csv');
      await grok.data.detectSemanticTypes(df);
      this.tableView!.dataFrame = df;
      await this.processFileData();
    }
  }
  
  async processFile(file: File) {
    if (!file) return;
    const csvData = await file.text();
    const table = DG.DataFrame.fromCsv(csvData);
    await grok.data.detectSemanticTypes(table);
    this.tableView!.dataFrame = table;
  
    const splashScreen = this.buildSplash(this.formContainer, 'Calculating...');
    const models = await getQueryParams();
    await performChemicalPropertyPredictions(
      this.tableView!.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE)!,
      this.tableView!.dataFrame,
      models,
      undefined,
      false,
      false,
      true
    );
    splashScreen.close();
  
    await this.processFileData();
  }
  
  private async processFileData(): Promise<void> {  
    await grok.data.detectSemanticTypes(this.tableView!.dataFrame);
    const models = await getQueryParams();
    const molColName = this.tableView!.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE)!.name;
    const molIdx = this.tableView!.dataFrame.columns.names().findIndex(c => c === molColName);
    let i = molIdx + 2;
    
    await addColorCoding(this.tableView!.dataFrame, models.split(','));
    await addSparklines(this.tableView!.dataFrame, models.split(','), i, 'ADMET');
    i += 1;
    await setProperties();
  
    const uniqueSubgroupNames: string[] = Array.from(
      new Set(properties.subgroup.map((subg: any) => subg.name))
    );
    
    for (const subgroupName of uniqueSubgroupNames) {
      const subgroupModels = properties.subgroup
        .filter((subg: any) => subg.name === subgroupName)
        .flatMap((subg: any) => subg.models.map((model: any) => model.name));
  
      await addSparklines(this.tableView!.dataFrame, subgroupModels, i, subgroupName);
      i += 1;
    }
  
    this.tableView!.grid.scrollToCell(molColName, 0);
  }  

  refresh(table: DG.DataFrame, modeContainer: HTMLElement, ratio: number) {
    if (table.rowCount > 0)
      this.tableView?.dockManager.dock(modeContainer, DG.DOCK_TYPE.TOP, null, '', ratio);
  }
  
  buildSplash(root: HTMLElement, description: string): ISplash {
    const indicator = ui.loader();
    const panel = ui.divV([indicator, ui.p(description)], { classes: 'splash-panel' });
    const loaderEl = ui.div([panel], { classes: 'splash-container' });
    root.append(loaderEl);
    return { el: loaderEl, close: () => loaderEl.remove() };
  }
}

async function openMoleculeDataset(name: string): Promise<DG.DataFrame> {
  const table = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(name));
  return table;
}