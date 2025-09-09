/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ErrorHandlingLabels, ScopeLabels, ScopeLabelsReduced } from './constants';
import '../../css/moltrack.css';
import { renderMappingEditor, TargetProperty } from '../components/mapping_editor';
import { MolTrackDockerService } from './moltrack-docker-service';
import {UiUtils} from '@datagrok-libraries/compute-utils';
import { fetchSchema } from '../package';
import { Subscription } from 'rxjs';

let openedView: DG.ViewBase | null = null;

export class RegistrationView {
  view: DG.View;

  entityTypeInput: DG.InputBase | null = null;
  errorStrategyInput: DG.InputBase | null = null;
  registerButton: HTMLButtonElement | null = null;

  mappingEditorDiv: HTMLDivElement = ui.div();
  uploadedPreviewDiv: HTMLDivElement;
  // summaryDiv: HTMLDivElement = ui.div([]);
  rightPanel: HTMLDivElement | null = null;

  uploadedDf: DG.DataFrame | null = null;
  mappingDict: Record<string, string> = {};
  subs: Subscription[] = [];

  hasErrors: boolean = false;
  registrationStarted: boolean = false;

  title: string = 'Bulk registration';
  messageContainer: HTMLDivElement = ui.div([], 'moltrack-info-container');

  constructor() {
    this.view = DG.View.create();
    this.view.name = 'Bulk Registration';

    this.createInputs();
    this.addTitle();
    this.uploadedPreviewDiv = ui.div('', 'moltrack-register-preview-div');
    // this.summaryDiv.style.marginTop = '50px';
    this.buildUI();
    this.addRibbonButtons();
    this.addSubs();
  }

  private async addTitle() {
    const titleText = ui.divText(this.title, 'moltrack-title');
    this.messageContainer.append(titleText);
  }

  private async createInputs() {
    this.entityTypeInput = this.createChoiceInput(
      'Register',
      Object.keys(ScopeLabels).filter((k) => k !== 'Assays'),
      (value) => this.createMapping(value),
    );

    this.errorStrategyInput = this.createChoiceInput('On error', Object.keys(ErrorHandlingLabels));
  }

  private buildUI() {
    this.createMapping();

    const inputRow = ui.wideForm([
      this.entityTypeInput!,
      this.errorStrategyInput!,
    ], 'moltrack-input-row');

    const leftPanel = ui.divV([inputRow, this.mappingEditorDiv/*, this.summaryDiv*/]);

    const dragDropInput = this.createFileInputPane();

    this.rightPanel = ui.divV([
      dragDropInput,
    ]);

    const container = ui.splitH([
      ui.box(leftPanel, {style: {display: 'flex', flex: '1'}}),
      ui.box(this.rightPanel, {style: {display: 'flex', flex: '2'}}),
      // leftPanel,
      // this.rightPanel,
    ], {}, true);
    container.classList.add('moltrack-container');
    this.view.root.append(ui.divV([this.messageContainer, container]));
  }

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

    if (this.rightPanel?.firstElementChild)
      (this.rightPanel.firstElementChild as HTMLElement).style.height = '5%';

    this.uploadedPreviewDiv.style.height = '95%';
    this.rightPanel?.appendChild(this.uploadedPreviewDiv);

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
    this.registerButton = ui.bigButton('REGISTER', () => {});

    this.registerButton!.addEventListener('click', async (e) => {
      if (this.hasErrors) {
        e.preventDefault();
        e.stopPropagation();
      } else
        await this.registerEntities();
    });

    this.registerButton.onmouseenter = (e) => {
      let tooltipMessage: string | null = null;

      if (this.hasErrors) {
        tooltipMessage = 'Form validation failed. Please check your input';
        this.registerButton!.style.cursor = 'not-allowed';
      } else if (this.registrationStarted) {
        tooltipMessage = 'Registration is already in progress. Please wait…';
        this.registerButton!.style.cursor = 'not-allowed';
      } else
        this.registerButton!.style.cursor = 'pointer';

      if (tooltipMessage)
        ui.tooltip.show(tooltipMessage, e.clientX + 15, e.clientY);
    };

    this.registerButton.onmouseleave = (e) => ui.tooltip.hide();

    const addToWorkspaceButton = ui.icons.add(() => {
      if (this.uploadedDf) grok.shell.addTablePreview(this.uploadedDf);
    }, 'Add uploaded dataset to workspace');

    this.view.setRibbonPanels([[this.registerButton, addToWorkspaceButton]]);
  }

  private createChoiceInput<T>(
    label: string,
    choices: T[],
    onChanged?: (value: T) => void,
  ): DG.ChoiceInput<T | null> {
    if (!choices || choices.length === 0)
      throw new Error('Choices array cannot be empty');

    return ui.input.choice(label, {
      value: choices[0],
      items: choices,
      onValueChanged: (value: T) => onChanged?.(value),
    });
  }

  private async registerEntities() {
    if (!this.uploadedDf) return;

    const csv = await this.uploadedDf.toCsv();
    const csvFile = DG.FileInfo.fromString('data.csv', csv);
    const loader = ui.loader();
    try {
      // ui.empty(this.summaryDiv);
      // this.summaryDiv.appendChild(ui.loader());
      this.registrationStarted = true;
      this.registerButton?.classList.add('dim');

      ui.empty(this.messageContainer);
      this.messageContainer.appendChild(loader);

      const df: DG.DataFrame = await grok.functions.call('Moltrack:registerBulk', {
        csvFile: csvFile,
        scope: ScopeLabels[this.entityTypeInput?.value],
        mapping: Object.keys(this.mappingDict).length === 0 ? '' : JSON.stringify(this.mappingDict),
        errorHandling: ErrorHandlingLabels[this.errorStrategyInput?.value],
      });

      let corporateId;
      if (Object.keys(ScopeLabelsReduced).includes(this.entityTypeInput?.value))
        corporateId = `corporate_${this.entityTypeInput?.value.toLowerCase().replace(/(es|s)$/, '')}_id`;

      this.uploadedDf.join(df, ['smiles'], ['smiles'], null, ['registration_status', 'registration_error_message', corporateId!], DG.JOIN_TYPE.INNER, true);

      this.createSummary();
    } catch (err: any) {
      grok.shell.error(`Registration failed: ${err.message}`);
    } finally {
      this.registerButton?.classList.remove('dim');
      this.registrationStarted = false;
      this.messageContainer.removeChild(loader);
    }
  }

  // private createSummary() {
  //   if (!this.uploadedDf) return;
  //   ui.empty(this.summaryDiv);
  //   const statuses = this.uploadedDf.getCol('registration_status').toList();

  //   const successCount = statuses.filter((v) => v === 'success').length;
  //   const failedCount = statuses.filter((v) => v === 'failed').length;
  //   const notProcessedCount = statuses.filter((v) => v === 'not_processed').length;

  //   const summary: { [key: string]: number } = {
  //     'Successful': successCount,
  //     'Failed': failedCount,
  //     'Not processed': notProcessedCount,
  //   };

  //   const summaryTable = ui.table(
  //     [summary],
  //     (item) => [item['Successful'], item['Failed'], item['Not processed']],
  //     Object.keys(summary),
  //   );

  //   this.summaryDiv.appendChild(summaryTable);
  // }

  private createSummary() {
    if (!this.uploadedDf) return;

    const statuses = this.uploadedDf.getCol('registration_status').toList();
    const successCount = statuses.filter((v) => v === 'success').length;
    const failedCount = statuses.filter((v) => v === 'failed').length;
    const notProcessedCount = statuses.filter((v) => v === 'not_processed').length;

    const header = failedCount === 0 ?
      `✅ Bulk registration completed` :
      `⚠ Bulk registration completed with errors`;
    const status = failedCount === 0 ? 'success' : 'failed';

    const messageParts: string[] = [];

    if (successCount > 0) messageParts.push(`${successCount} registered`);
    if (failedCount > 0) messageParts.push(`${failedCount} failed`);
    if (notProcessedCount > 0) messageParts.push(`${notProcessedCount} skipped`);

    const message = messageParts.join(', ');

    ui.empty(this.messageContainer);
    const infoDiv = ui.info(message, header, true);
    const bar = infoDiv.querySelector('.grok-info-bar') as HTMLElement;
    if (bar) {
      bar.classList.toggle('moltrack-bar-success', status === 'success');
      bar.classList.toggle('moltrack-bar-error', status !== 'success');
    }
    this.messageContainer.appendChild(infoDiv);
  }

  private async createMapping(value?: string) {
    const parsedProps = JSON.parse(await fetchSchema());
    const targetProperties: TargetProperty[] = [
      ...parsedProps
        .filter((p: any) => ScopeLabels[value ?? this.entityTypeInput?.value].includes(p.entity_type.toLowerCase()))
        .map((p: any) => ({
          name: p.friendly_name ?? p.name,
          min: p.min,
          max: p.max,
          semType: p.semantic_type ? p.semantic_type['name'] : null,
          type: p.value_type,
          required: false,
        })),
      { name: 'smiles', required: true, semType: DG.SEMTYPE.MOLECULE },
    ];

    // TODO: Auto-mapping should handle cases where some properties may not exist in the database
    let autoMapping: Map<string, string> = new Map();
    if (this.uploadedDf) {
      await MolTrackDockerService.init();
      autoMapping = await MolTrackDockerService.getAutoMapping(this.uploadedDf!.columns.names(), 'COMPOUND');
    }

    const sourceColumns = this.uploadedDf ? this.uploadedDf.columns.names() : [''];
    const mappings = new Map<string, string>();
    for (const [source, target] of Object.entries(autoMapping)) {
      const cleanTarget = target.replace(/^.*?_details\./, '');
      const exists = targetProperties.some((tp: TargetProperty) => tp.name === cleanTarget);
      if (exists) mappings.set(source, cleanTarget);
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
  }

  private addSubs() {
    this.subs.push(grok.events.onCustomEvent('mappingValidationChanged').subscribe((args) => {
      this.hasErrors = args.hasErrors;
      this.registerButton?.classList.toggle('dim', this.hasErrors);
    }));
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
