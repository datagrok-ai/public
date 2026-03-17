/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FileInputUtils} from '@datagrok-libraries/tutorials/src/utils/file-input-utils';

import {ErrorHandlingLabels, MOLTRACK_MAPPING_VALIDATION_CHANGED, MOLTRACK_REQUEST_TITLE_UPDATE, Scope, ScopeLabels, ScopeLabelsReduced, ScopeToEntityType} from '../utils/constants';
import {renderMappingEditor, TargetProperty} from '../components/mapping_editor';
import {MolTrackDockerService} from '../services/moltrack-docker-service';
import {fetchSchema} from '../package';
import {Subscription} from 'rxjs';

import '../../css/moltrack.css';

let openedView: DG.ViewBase | null = null;

export class RegistrationView {
  view: DG.View;

  entityTypeInput: DG.InputBase | null = null;
  errorStrategyInput: DG.InputBase | null = null;
  registerButton: HTMLButtonElement | null = null;

  mappingEditorDiv: HTMLDivElement = ui.div();
  uploadedPreviewDiv: HTMLDivElement;
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
    this.buildUI();
    this.addRibbonButtons();
    this.addSubs();
  }

  private async addTitle() {
    ui.empty(this.messageContainer);
    const titleText = ui.divText(this.title, 'moltrack-title');
    this.messageContainer.append(titleText);
  }

  private async createInputs() {
    this.entityTypeInput = this.createChoiceInput(
      'Register',
      Object.keys(ScopeLabels).filter((k) => k !== 'Assays'),
      (value) => {
        this.createMapping(value);
        requestTitleUpdate();
      },
    );

    this.errorStrategyInput = this.createChoiceInput(
      'On error',
      Object.keys(ErrorHandlingLabels),
      (value) => requestTitleUpdate(),
    );
  }

  private buildUI() {
    this.createMapping();

    const inputRow = ui.wideForm([
      this.entityTypeInput!,
      this.errorStrategyInput!,
    ], 'moltrack-input-row');

    this.mappingEditorDiv.classList.add('moltrack-mapping-editor');
    const leftPanel = ui.divV([inputRow, this.mappingEditorDiv]);
    leftPanel.classList.add('moltrack-left-panel');

    const dragDropInput = FileInputUtils.createFileInputPane(async (file: File) => {
      await this.loadFile(file);
      requestTitleUpdate();
    });
    const inputEditor = dragDropInput.querySelector('.ui-input-editor') as HTMLElement;
    inputEditor.style.border = '';

    this.rightPanel = ui.divV([
      dragDropInput,
    ]);

    const container = ui.splitH([
      ui.box(leftPanel, {style: {display: 'flex', flex: '1'}}),
      ui.box(this.rightPanel, {style: {display: 'flex', flex: '2'}}),
    ], {}, true);
    container.classList.add('moltrack-container');
    this.view.root.append(ui.divV([this.messageContainer, container]));
  }

  private async loadFile(file: File) {
    ui.empty(this.uploadedPreviewDiv);
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

      const value = this.entityTypeInput?.value ?? '';
      const scope = ScopeLabels[value];
      const corporateId = value in ScopeLabelsReduced ? `corporate_${ScopeToEntityType[scope]}_id` : '';

      const joinColumns = [
        'registration_status',
        'registration_error_message',
        ...(corporateId && df.columns.names().includes(corporateId) ? [corporateId] : []),
      ];

      const renamedColumns: string[] = [];
      for (const col of joinColumns) {
        const newName = this.uploadedDf.columns.getUnusedName(col);
        df.col(col)!.name = newName;
        renamedColumns.push(newName);
      }

      const firstColName = this.uploadedDf.columns.names()[0];

      this.uploadedDf.join(
        df,
        [firstColName],
        [firstColName],
        null,
        renamedColumns,
        DG.JOIN_TYPE.INNER,
        true,
      );

      const statusColName = renamedColumns[0];
      this.createSummary(statusColName);
    } catch (err: any) {
      grok.shell.error(`Registration failed: ${err.message}`);
    } finally {
      this.registerButton?.classList.remove('dim');
      loader.remove();
      this.registrationStarted = false;
    }
  }

  private createSummary(statusColName: string = 'registration_status') {
    if (!this.uploadedDf) return;

    const statuses = this.uploadedDf.getCol(statusColName).toList();
    const successCount = statuses.filter((v) => v === 'success').length;
    const failedCount = statuses.filter((v) => v === 'failed').length;
    const notProcessedCount = statuses.filter((v) => v === 'not_processed').length;
    const totalFailures = failedCount + notProcessedCount;

    const header = totalFailures === 0 ?
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

  private getScope(value?: string): Scope {
    return ScopeLabels[value ?? this.entityTypeInput?.value];
  }

  private buildTargetProperties(parsedProps: any[], scope: Scope): TargetProperty[] {
    const needsAssayName = scope === Scope.ASSAY_RUNS || scope === Scope.ASSAY_RESULTS;
    const isAssayResults = scope === Scope.ASSAY_RESULTS;
    const needsSmiles = scope === Scope.COMPOUNDS || scope === Scope.BATCHES;
    const entityScopes = isAssayResults ? [Scope.BATCHES, Scope.ASSAY_RUNS, Scope.ASSAY_RESULTS] : [scope];

    return [
      ...(needsAssayName ? [{name: 'assay.name', required: true, type: 'string'}] : []),
      ...entityScopes.flatMap((entityScope) => parsedProps
        .filter((p: any) => entityScope.includes(p.entity_type.toLowerCase()))
        .map((p: any) => {
          const propName = p.friendly_name ?? p.name;
          // For assay results, prefix cross-entity properties so the backend groups them correctly
          const prefix = isAssayResults && entityScope !== scope ? `${ScopeToEntityType[entityScope]}_details.` : '';
          return {
            name: `${prefix}${propName}`,
            min: p.min,
            max: p.max,
            semType: p.semantic_type ? p.semantic_type['name'] : null,
            type: p.value_type,
            required: false,
          };
        }),
      ),
      ...(needsSmiles ? [{name: 'smiles', required: true, semType: DG.SEMTYPE.MOLECULE}] : []),
    ];
  }

  private resolveAutoMapping(
    autoMapping: Record<string, string>, targetProperties: TargetProperty[],
  ): Map<string, string> {
    const mappings = new Map<string, string>();
    for (const [source, target] of Object.entries(autoMapping)) {
      const unprefixed = target.replace(/^.*?_details\./, '');
      const matched = targetProperties.find((tp) => tp.name === target || tp.name === unprefixed);
      if (matched)
        mappings.set(matched.name, source);
    }
    return mappings;
  }

  private async createMapping(value?: string) {
    this.mappingDict = {};
    const parsedProps = JSON.parse(await fetchSchema());
    const scope = this.getScope(value);
    const targetProperties = this.buildTargetProperties(parsedProps, scope);

    let autoMapping: Record<string, string> = {};
    if (this.uploadedDf)
      autoMapping = await MolTrackDockerService.getAutoMapping(this.uploadedDf.columns.names(), ScopeToEntityType[scope]);

    const sourceColumns = this.uploadedDf ? this.uploadedDf.columns.names() : [''];
    const mappings = this.resolveAutoMapping(autoMapping, targetProperties);

    // Passthrough targets: dotted keys (assay.name, batch_details.X) and smiles are sent as-is
    const isPassthrough = (target: string) => target.includes('.') || target === 'smiles';
    const handleMap = (target: string, source: string) => {
      if (!this.uploadedDf)
        return;
      this.mappingDict[source] = isPassthrough(target) ? target : `${ScopeToEntityType[scope]}_details.${target}`;
    };

    const handleUndo = (target: string) => {
      for (const [key, val] of Object.entries(this.mappingDict)) {
        if (val.endsWith(`.${target}`) || val === target)
          delete this.mappingDict[key];
      }
    };

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
    this.subs.push(grok.events.onCustomEvent(MOLTRACK_MAPPING_VALIDATION_CHANGED).subscribe((args) => {
      this.hasErrors = args.hasErrors;
      this.registerButton?.classList.toggle('dim', this.hasErrors);
    }));
    this.subs.push(grok.events.onCustomEvent(MOLTRACK_REQUEST_TITLE_UPDATE).subscribe((_) => this.addTitle()));
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}

export function requestTitleUpdate(): void {
  grok.events.fireCustomEvent(MOLTRACK_REQUEST_TITLE_UPDATE, {});
}
