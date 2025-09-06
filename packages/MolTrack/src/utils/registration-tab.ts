/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ErrorHandlingLabels, ScopeLabels } from './constants';
import '../../css/moltrack.css';
import { renderMappingEditor, TargetProperty } from '../components/mapping_editor';
import { _package, fetchCompoundProperties } from '../package';
import { MolTrackDockerService } from './moltrack-docker-service';
import {UiUtils} from '@datagrok-libraries/compute-utils';

let openedView: DG.ViewBase | null = null;

export class RegistrationView {
  view: DG.View;
  datasetInput: DG.InputBase | null = null;
  entityTypeInput: DG.InputBase | null = null;
  errorStrategyInput: DG.InputBase | null = null;
  mappingEditorDiv: HTMLDivElement = ui.div();
  uploadedPreviewDiv: HTMLDivElement;
  uploadedDf: DG.DataFrame | null = null;
  mappingDict: Record<string, string> = {};
  summaryDiv: HTMLDivElement = ui.div([]);

  constructor() {
    this.view = DG.View.create();
    this.view.name = 'Bulk Registration';

    this.createInputs();
    this.createSummaryDiv();
    this.uploadedPreviewDiv = ui.div('', 'moltrack-register-preview-div');
    this.buildUI();
    this.addRibbonButtons();
  }

  createSummaryDiv() {
    const summaryDiv: HTMLDivElement = ui.div([]);

    const placeholderText = ui.divText('Summary will appear here once you register entities.', {
      style: {
        textAlign: 'center',
        marginBottom: '8px',
      },
    });
    summaryDiv.appendChild(placeholderText);

    this.summaryDiv.style.cssText = `
  border: 1px solid #ccc;
  border-radius: 8px;
  background-color: #f9f9f9;
  min-height: 80px;
  max-width: 240px;
  margin-top: 40px;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;

`;

    this.summaryDiv.appendChild(summaryDiv);
  }

  private async createInputs() {
    this.entityTypeInput = this.createChoiceInput('Register', Object.keys(ScopeLabels));
    // this.datasetInput = ui.input.file('File', /*, {
    //   value: (await grok.dapi.files.list('System:AppData/MolTrack/samples/'))[1],
    // }*/);

    this.errorStrategyInput = this.createChoiceInput('On error', Object.keys(ErrorHandlingLabels));
  }

  private buildUI() {
    this.createMapping();

    const inputRow = ui.wideForm([
      this.datasetInput!,
      this.entityTypeInput!,
      this.errorStrategyInput!,
    ], 'moltrack-input-row');

    const leftPanel = ui.divV([inputRow, this.mappingEditorDiv, this.summaryDiv]);

    const dragDropInput = this.createFileInputPane();

    const rightPanel = ui.divV([
      dragDropInput,
      this.uploadedPreviewDiv,
    ], 'moltrack-right-panel');

    leftPanel.style.flex = '1';
    rightPanel.style.flex = '2';

    // const splitterBar = new SplitterBar(leftPanel, rightPanel, this.stackedVertical);

    // const container = ui.divH([leftPanel, rightPanel], 'moltrack-container');
    const container = ui.splitH([
      leftPanel,
      rightPanel,
    ], {}, true);
    container.classList.add('moltrack-container');
    this.view.root.append(container);

    // TODO: No need to preload default table - we gonna add just a drag-n-drop with an action item
    // this.preloadDefaultTable();
  }

  /** Copy pasted from the Admetica demo app */
  private createFileInputPane() {
    const fileInputEditor = this.initializeFileInputEditor();
    this.removeLabels(fileInputEditor.root);
    this.styleInputEditor(fileInputEditor.root);
    this.setupDragAndDrop(fileInputEditor.root);
    this.removeOptionsIcon(fileInputEditor.root);
    fileInputEditor.root.classList.add('demo-file-input');
    return ui.divV([fileInputEditor], {classes: 'demo-file-input-container'});
  }

  private initializeFileInputEditor() {
    const fileInputEditor = UiUtils.fileInput('', null, async (file: File) => {
      await this.loadFile(file);
    }, null);
    fileInputEditor.stringValue = 'Drag and drop a CSV file here, or click to select a file.';
    return fileInputEditor;
  }

  private async loadFile(file: File) {
    if (!file) return;

    const extension = file.name.split('.').pop()!;
    if (!['csv'].includes(extension)) {
      grok.shell.info(`The file extension ${extension} is not supported!`);
      return;
    }

    const csvData = await file.text();
    this.uploadedDf = DG.DataFrame.fromCsv(csvData);
    await grok.data.detectSemanticTypes(this.uploadedDf);
    this.createMapping();
    this.uploadedPreviewDiv.append(this.uploadedDf.plot.grid().root);
  }

  private removeLabels(root: HTMLElement): void {
    const labelSelectors = '.ui-label, .ui-input-label';
    const labels = root.querySelectorAll<HTMLElement>(labelSelectors);
    labels.forEach((label) => label.remove());
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
      fontSize: '14px',
      cursor: 'pointer',
      textAlign: 'center',
      borderBottom: 'none',
    });
  }

  private setupDragAndDrop(root: HTMLElement): void {
    const inputEditor = root.querySelector<HTMLElement>('.ui-input-editor');
    if (!inputEditor) return;

    const highlightColor = '#e0f7fa';
    const defaultColor = '#ffffff';
    inputEditor.style.border = '1px dashed #007bff';

    const setHighlightedStyle = () => inputEditor.style.backgroundColor = highlightColor;
    const resetStyle = () => inputEditor.style.backgroundColor = defaultColor;

    inputEditor.addEventListener('dragenter', setHighlightedStyle);
    inputEditor.addEventListener('dragover', (event) => {
      event.preventDefault();
      setHighlightedStyle();
    });

    inputEditor.addEventListener('dragleave', resetStyle);
    inputEditor.addEventListener('drop', (event) => {
      event.preventDefault();
      resetStyle();
    });
  }

  private removeOptionsIcon(root: HTMLElement): void {
    const optionsIcon = root.querySelector<HTMLElement>('.ui-input-options .grok-icon.fal.fa-cloud-upload');
    optionsIcon?.remove();
  }

  private addRibbonButtons() {
    const registerButton = ui.bigButton('REGISTER', async () => await this.registerEntities());

    const addToWorkspaceButton = ui.icons.add(() => {
      if (this.uploadedDf) grok.shell.addTablePreview(this.uploadedDf);
    }, 'Add uploaded dataset to workspace');

    this.view.setRibbonPanels([[registerButton, addToWorkspaceButton]]);
  }

  private createChoiceInput<T>(label: string, choices: T[]): DG.ChoiceInput<T | null> {
    return ui.input.choice(label, {
      value: choices[0],
      items: choices,
    });
  }

  // private async preloadDefaultTable() {
  //   try {
  //     this.uploadedDf = await grok.data.loadTable(_package.webRoot + 'files/samples/compounds.csv');
  //     await grok.data.detectSemanticTypes(this.uploadedDf);

  //     ui.empty(this.uploadedPreviewDiv);
  //     this.uploadedPreviewDiv.append(this.uploadedDf.plot.grid().root);
  //   } catch (e: any) {
  //     grok.shell.error(`Failed to load dataset: ${e.message}`);
  //   }
  // }

  private async registerEntities() {
    if (!this.uploadedDf) return;

    const csv = await this.uploadedDf.toCsv();
    const csvFile = DG.FileInfo.fromString('data.csv', csv);
    try {
      ui.empty(this.summaryDiv);
      this.summaryDiv.appendChild(ui.loader());
      const df: DG.DataFrame = await grok.functions.call('Moltrack:registerBulk', {
        csvFile: csvFile,
        scope: ScopeLabels[this.entityTypeInput?.value],
        mapping: Object.keys(this.mappingDict).length === 0 ? '' : JSON.stringify(this.mappingDict),
        errorHandling: ErrorHandlingLabels[this.errorStrategyInput?.value],
      });

      // TODO: Add ui.loader while we are waiting for the results
      this.uploadedDf.join(df, ['smiles'], ['smiles'], null, ['registration_status', 'registration_error_message'], DG.JOIN_TYPE.INNER, true);

      this.createSummary();
      // TODO: Add summary with successful/failed rows

      grok.shell.info('Entities registered successfully!');
    } catch (err: any) {
      grok.shell.error(`Registration failed: ${err.message}`);
    }
  }

  private createSummary() {
    if (!this.uploadedDf) return;
    ui.empty(this.summaryDiv);
    const statuses = this.uploadedDf.getCol('registration_status').toList();

    const successCount = statuses.filter((v) => v === 'success').length;
    const failedCount = statuses.filter((v) => v === 'failed').length;
    const notProcessedCount = statuses.filter((v) => v === 'not_processed').length;

    // const successText = ui.divText('Successful', {
    //   style: { color: '#286344' },
    // });

    // const failedText = ui.divText('Failed', {
    //   style: { color: '#763434' },
    // });

    // const notProcessedText = ui.divText('Not processed', {
    //   style: { color: '#805125' },
    // });

    const summary: { [key: string]: number } = {
      'Successful': successCount,
      'Failed': failedCount,
      'Not processed': notProcessedCount,
    };
    const summaryTable = ui.tableFromMap(summary);
    this.summaryDiv.appendChild(summaryTable);
  }

  private async createMapping() {
    const parsedProps = JSON.parse(await fetchCompoundProperties());
    const targetProperties: TargetProperty[] = [
      ...parsedProps.map((p: any) => ({
        name: p.friendly_name ?? p.name,
        type: p.value_type,
        required: false,
      })),
      { name: 'smiles', required: true },
    ];

    let autoMapping: Map<string, string> = new Map();
    if (this.uploadedDf) {
      await MolTrackDockerService.init();
      autoMapping = await MolTrackDockerService.getAutoMapping(this.uploadedDf!.columns.names(), 'COMPOUND');
    }

    const sourceColumns = this.uploadedDf ? this.uploadedDf.columns.names() : [];
    const mappings = new Map<string, string>();
    for (const [source, target] of Object.entries(autoMapping)) {
      const cleanTarget = target.replace(/^.*?_details\./, '');
      mappings.set(source, cleanTarget);
    }

    const handleMap = (target: string, source: string) => {
      if (!this.uploadedDf) return;

      const col = this.uploadedDf.col(source);
      if (col) {
        if (target === 'smiles') {
          this.mappingDict[source] = target;
          return;
        }
        const newName = `${this.entityTypeInput?.value.toLowerCase().replace(/(es|s)$/, '')}_details.${target}`;
        this.mappingDict[source] = newName;
      }
    };

    const handleUndo = () => grok.shell.info('Mapping undone');

    renderMappingEditor(this.mappingEditorDiv, {
      targetProperties,
      sourceColumns,
      mappings,
      onMap: handleMap,
      onUndo: handleUndo,
    },
    this.uploadedDf!,
    );

    // TODO: near the render mapping editor there should be a separate column with the mapping status (warning, error, success)
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
