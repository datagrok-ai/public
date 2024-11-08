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
  placeholder: HTMLElement;
  sketcherInstance: grok.chem.Sketcher;

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
    this.sketcherInstance = new grok.chem.Sketcher();
    
    // Placeholder with unified style
    this.placeholder = ui.divText('Please sketch a molecule to see the predictions.', {
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

    // /this.showPlaceholder();
    this.formContainer.appendChild(this.placeholder);
    this.container.appendChild(this.modeContainer);
    this.container.appendChild(this.formContainer);
  }

  async init() {
    this.sketcherInstance.setMolecule('CC(Oc1ccccc1C(O)=O)=O');
    this.sketcherInstance.onChanged.subscribe(async () => {
      const smiles = this.sketcherInstance.getSmiles();
      await this.onChanged(smiles)
    });
    this.sketcherInstance.molInput.value = 'aspirin';
    this.sketcher = this.sketcherInstance.root;
    this.sketcherDiv.appendChild(this.sketcher);

    // Mode selection with unified tabs style
    const modeTabs = ui.tabControl();
    modeTabs.onTabChanged.subscribe((tab: DG.TabPane) => {
      this.mode = tab.header.innerText.split(' ')[0].toLowerCase();
      this.clearTable();
      this.clearForm();
      if (this.mode === 'sketch') {
        this.formContainer.style.visibility = 'visible';
        this.formContainer.style.position = 'relative'; // Restore original layout positioning
        if (this.tableView?.grid) {
          this.tableView.grid.root.style.visibility = 'hidden';
          this.tableView.dockManager.dock(this.container, DG.DOCK_TYPE.TOP, null, '', 0.99);
        }
        this.sketcherInstance.onChanged.next();
      } else {
        this.formContainer.style.visibility = 'hidden';
        this.formContainer.style.position = 'absolute';
        if (this.tableView?.grid) {
          this.tableView.grid.root.style.visibility = 'visible';
          this.tableView.dockManager.dock(this.container, DG.DOCK_TYPE.TOP, null, '', 0.2 );
          this.tableView.grid.root.style.visibility = 'visible';
          this.tableView.dockManager.dock(this.container, DG.DOCK_TYPE.TOP, null, '', 0.2 );
          openMoleculeDataset('System:AppData/Admetica/demo_files/sar-small_app.csv').then((df) => {
            this.tableView!.dataFrame = df;
            grok.data.detectSemanticTypes(this.tableView!.dataFrame).then(() => {
            const molColName = this.tableView!.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE)!.name;
            const molIdx = this.tableView!.dataFrame.columns.names().findIndex(c => c === molColName);
            let i = molIdx + 2;
            getQueryParams().then((models) => {
              addColorCoding(this.tableView!.dataFrame, models.split(','));
              setProperties().then(() => {
                addSparklines(this.tableView!.dataFrame, models.split(','), i, 'ADMET');
                i += 1;
                const uniqueSubgroupNames: string[] = Array.from(
                  new Set(properties.subgroup.map((subg: any) => subg.name))
                );
        
                uniqueSubgroupNames.forEach(async (subgroupName: string) => {
                  const models = properties.subgroup
                    .filter((subg: any) => subg.name === subgroupName) // Filter subgroups with the current name
                    .flatMap((subg: any) => subg.models.map((model: any) => model.name)); // Map to model names and flatten
                    
                  await addSparklines(this.tableView!.dataFrame, models, i, subgroupName);
                  i += 1;
                });
                this.tableView!.grid.scrollToCell(molColName, 0);
              });
            })
          });
          })
        }
      }
    });    

    modeTabs.addPane('Sketch', () => {
      this.mode = 'sketch';
      this.formContainer.style.visibility = 'visible';
      this.formContainer.style.position = 'relative'; // Restore original layout positioning
      if (this.tableView?.grid) {
        this.tableView.grid.root.style.visibility = 'hidden';
        this.tableView.dockManager.dock(this.container, DG.DOCK_TYPE.TOP, null, '', 0.99);
      }
      return this.sketcherDiv;
    });
    
    modeTabs.addPane('File', () => {
      this.mode = 'file';
      this.formContainer.style.visibility = 'hidden';
      this.formContainer.style.position = 'absolute';
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
      //this.tableView!.grid.root.style.display = 'none';
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
    //this.formContainer.appendChild(this.placeholder);
  }

  async onChanged(smiles: string) {
    if (this.isTaskActive) {
      this.isTaskActive = false;
      return;
    }

    this.isTaskActive = true; // Mark the task as active
    try {
      if (this.mode === 'sketch') {

        if (!smiles) {
          this.clearForm();
          //this.showPlaceholder();
          return;
        }

        /*if (smiles && this.formContainer.innerHTML === '') {
          this.showPlaceholder();
        }*/

        this.clearTable();
        const col = this.tableView?.dataFrame.columns.getOrCreate('smiles', 'string', 1);
        col!.semType = DG.SEMTYPE.MOLECULE;
        this.tableView?.dataFrame.set('smiles', 0, smiles);
        await grok.data.detectSemanticTypes(this.tableView!.dataFrame);

        const models = await getQueryParams();
        const splashScreen = this.buildSplash(this.formContainer, 'Calculating...');

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
        await addSparklines(this.tableView!.dataFrame, models.split(','), molIdx! + 1);
        splashScreen.close();
        
        const form = createDynamicForm(this.tableView!.dataFrame, models.split(','), 'smiles', true);
        this.formContainer.innerHTML = ''; 
        this.formContainer.appendChild(form.root);
        this.tableView?.grid.invalidate();
        //this.hidePlaceholder();
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
        backgroundColor: '#ffffff',
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
  
    const splashScreen = this.buildSplash(this.formContainer, 'Calculating...');
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
    let i = molIdx + 2;
    await addSparklines(this.tableView!.dataFrame, models.split(','), i, 'ADMET');
    i += 1;
    await setProperties();
    const uniqueSubgroupNames: string[] = Array.from(
      new Set(properties.subgroup.map((subg: any) => subg.name))
    );
    console.log('unique subgroups');
    console.log(uniqueSubgroupNames);
    uniqueSubgroupNames.forEach(async (subgroupName: string) => {
      const models = properties.subgroup
        .filter((subg: any) => subg.name === subgroupName) // Filter subgroups with the current name
        .flatMap((subg: any) => subg.models.map((model: any) => model.name)); // Map to model names and flatten
      
      await addSparklines(this.tableView!.dataFrame, models, i, subgroupName);
      i += 1;
    });
    this.tableView!.grid.scrollToCell(molColName, 0);
  }

  showPlaceholder() {
    this.formContainer.appendChild(this.placeholder);
  }

  hidePlaceholder() {
    if (this.placeholder.parentElement) this.formContainer.removeChild(this.placeholder);
  }

  refresh(table: DG.DataFrame, modeContainer: HTMLElement, ratio: number) {
    if (table.rowCount > 0) {
      this.tableView?.dockManager.dock(modeContainer, DG.DOCK_TYPE.TOP, null, '', ratio);
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

async function openMoleculeDataset(name: string): Promise<DG.DataFrame> {
  const table = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(name));
  return table;
}