import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import { addSparklines, createDynamicForm, getQueryParams, performChemicalPropertyPredictions, properties, setProperties } from './admetica-utils';
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
  mode: string = 'single';
  container: HTMLElement;
  formContainer: HTMLElement;
  modeContainer: HTMLElement;
  sketcherDiv: HTMLElement;
  sketcher?: HTMLElement;
  placeholder: HTMLElement;

  private taskQueue: (() => Promise<void>)[] = [];
  private isTaskActive: boolean = false;

  constructor(parentCall: DG.FuncCall) {
    this.parentCall = parentCall;
    
    // Create main container with unified style
    this.container = ui.divH([], { classes: 'app-container' });
    
    // Create form and mode containers with reusable styling
    this.formContainer = ui.box(null, { classes: 'content-container', style: { height: '100%' } });
    this.sketcherDiv = ui.div([], { classes: 'content-container', style: { border: 'none' } });
    this.modeContainer = ui.box(null, { classes: 'mode-container', style: { height: '100%' } });
    
    // Placeholder with unified style
    this.placeholder = ui.divText('Please select a molecule or upload a file to see the predictions.', {
      classes: 'placeholder-text'
    });

    // Insert keyframes for animations
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
    const sketcherInstance = new grok.chem.Sketcher();
    sketcherInstance.setMolecule('CC(Oc1ccccc1C(O)=O)=O');
    sketcherInstance.onChanged.subscribe(async () => {
      const smiles = sketcherInstance.getSmiles();
      await this.onChanged(smiles)
    });
    this.sketcher = sketcherInstance.root;
    this.sketcherDiv.appendChild(this.sketcher);

    // Mode selection with unified tabs style
    const modeTabs = ui.tabControl();
    modeTabs.onTabChanged.subscribe((tab: DG.TabPane) => {
      this.clearTable();
      this.clearForm();
      this.mode = tab.header.innerText.split(' ')[0].toLowerCase();
    
      if (this.mode === 'file') {
        // Hide formContainer in "File Input" mode
        this.formContainer.style.visibility = 'hidden';
        this.formContainer.style.position = 'absolute';
      } else {
        // Show formContainer in "Single Molecule" mode
        this.formContainer.style.visibility = 'visible';
        this.formContainer.style.position = 'relative'; // Restore original layout positioning
      }
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
    this.prepareTableView();
  }

  prepareTableView() {
    const table = DG.DataFrame.create(1);
    this.tableView = DG.TableView.create(table, false);
    this.tableView.parentCall = this.parentCall;
    this.tableView.dataFrame.name = 'Table';
    
    setTimeout(async () => {
      this.tableView!._onAdded();
      await this.refresh(table, this.container);
      (grok.shell.view(DG.VIEW_TYPE.BROWSE) as DG.BrowseView).preview = this.tableView!;
    }, 300);
  }

  clearTable() {
    if (this.tableView) this.tableView.dataFrame = DG.DataFrame.create(1);
  }

  clearForm() {
    this.formContainer.innerHTML = '';
    this.formContainer.appendChild(this.placeholder);
  }

  async onChanged(smiles: string) {
    if (this.isTaskActive) {
      this.isTaskActive = false;
      return;
    }

    this.isTaskActive = true; // Mark the task as active
    try {
      if (this.mode === 'single') {
        this.clearTable();
        const col = this.tableView?.dataFrame.columns.getOrCreate('smiles', 'string', 1);
        col!.semType = DG.SEMTYPE.MOLECULE;
        this.tableView?.dataFrame.set('smiles', 0, smiles);
        await grok.data.detectSemanticTypes(this.tableView!.dataFrame);

        const models = await getQueryParams();
        const splashScreen = this.buildSplash(this.formContainer, 'Processing...');

        // Perform chemical property predictions without passing signal
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
        await addSparklines(this.tableView!.dataFrame, models.split(','), molIdx!);
        splashScreen.close();
        
        const form = createDynamicForm(this.tableView!.dataFrame, models.split(','), 'smiles', true);
        this.formContainer.innerHTML = ''; 
        this.formContainer.appendChild(form.root);
        this.hidePlaceholder();
      }
    } catch (error) {
      console.error('Error in prediction task:', error);
    } finally {
      this.isTaskActive = false; // Reset the task flag
    }
  }

  addTaskToQueue(task: () => Promise<void>) {
    // Add the task to the queue and immediately attempt to execute it
    this.taskQueue.push(task);
    this.processNextTask();
  }

  processNextTask() {
    if (this.taskQueue.length === 0 || this.isTaskActive) return;
    const nextTask = this.taskQueue.shift();
    if (nextTask) {
      nextTask().finally(() => this.processNextTask());
    }
  }

  createFileInputPane() {
    const fileInputEditor = UiUtils.fileInput('', null, async (file: File) => {
      await this.processFile(file);
    }, null);

    const labels = fileInputEditor.root.querySelectorAll('.ui-label, .ui-input-label');
    labels.forEach(label => label.remove());

    const inputEditor = fileInputEditor.root.querySelector('.ui-input-editor') as HTMLElement;
    if (inputEditor) {
      Object.assign(inputEditor.style, {
        width: '100%',
        height: '100%',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        backgroundColor: '#f9f9f9',
        color: '#007bff',
        fontSize: '22px',
        cursor: 'pointer',
        textAlign: 'center'
      });

      inputEditor.addEventListener('dragenter', () => {
        inputEditor.style.backgroundColor = '#e0f7fa';
      });
      
      inputEditor.addEventListener('dragleave', () => {
        inputEditor.style.backgroundColor = '#f9f9f9';
      });
      
      inputEditor.addEventListener('dragover', (event) => {
        inputEditor.style.backgroundColor = '#e0f7fa';
      });
      
      inputEditor.addEventListener('drop', () => {
        inputEditor.style.backgroundColor = '#f9f9f9';
      });
    }

    const optionsIcon = fileInputEditor.root.querySelector('.ui-input-options .grok-icon.fal.fa-cloud-upload') as HTMLElement;
    if (optionsIcon)
      optionsIcon.remove();
  
    fileInputEditor.root.classList.add('file-input');
    return ui.divV([fileInputEditor], { classes: 'file-input-container' });
  }

  async processFile(file: File) {
    const csvData = await file.text();
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
      false,
      false,
      true
    );
    splashScreen.close();

    const molColName = table.columns.bySemType(DG.SEMTYPE.MOLECULE)!.name;
    const molIdx = table.columns.names().findIndex(c => c === molColName);
    await setProperties();
    const uniqueSubgroupNames: string[] = Array.from(
      new Set(properties.subgroup.map((subg: any) => subg.name))
    );
    
    uniqueSubgroupNames.forEach(async (subgroupName: string) => {
      const models = properties.subgroup
        .filter((subg: any) => subg.name === subgroupName) // Filter subgroups with the current name
        .flatMap((subg: any) => subg.models.map((model: any) => model.name)); // Map to model names and flatten
      
      await addSparklines(this.tableView!.dataFrame, models, molIdx, subgroupName);
    });
  }   

  hidePlaceholder() {
    if (this.placeholder.parentElement) this.formContainer.removeChild(this.placeholder);
  }

  refresh(table: DG.DataFrame, modeContainer: HTMLElement) {
    if (table.rowCount > 0) {
      this.tableView?.dockManager.dock(modeContainer, DG.DOCK_TYPE.TOP, null, '', 0.65);
    }
  }
  
  buildSplash(root: HTMLElement, description: string): ISplash {
    const indicator = ui.loader();
    const panel = ui.divV([indicator, ui.p(description)], { classes: 'splash-panel' });
    const loaderEl = ui.div([panel], { classes: 'splash-container' });
    root.append(loaderEl);
  
    return { el: loaderEl, close: () => loaderEl.remove() };
  }
  
}