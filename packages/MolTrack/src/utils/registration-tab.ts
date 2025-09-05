/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ErrorHandling, ErrorHandlingLabels, Scope, ScopeLabels } from './constants';
import '../../css/moltrack.css';
import { renderMappingEditor, TargetProperty } from '../components/mapping_editor';
import { _package, fetchCompoundProperties } from '../package';
import { MolTrackDockerService } from './moltrack-docker-service';

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

  constructor() {
    this.view = DG.View.create();
    this.view.name = 'Entity Registration';

    this.createInputs();
    this.uploadedPreviewDiv = ui.div('', 'moltrack-register-preview-div');
    this.buildUI();
    this.addRibbonButtons();
  }

  private async createInputs() {
    this.entityTypeInput = this.createChoiceInput('Register', Object.keys(ScopeLabels));
    this.datasetInput = ui.input.file('File'/*, {
      value: (await grok.dapi.files.list('System:AppData/MolTrack/samples/'))[1],
    }*/);

    this.errorStrategyInput = this.createChoiceInput('On error', Object.keys(ErrorHandlingLabels));
  }

  private buildUI() {
    this.createMapping();
    const leftPanel = ui.divV([this.mappingEditorDiv], 'moltrack-left-panel');

    const inputRow = ui.divH([
      this.datasetInput!.root,
      this.entityTypeInput!.root,
      this.errorStrategyInput!.root,
    ], 'moltrack-input-row');

    const rightPanel = ui.divV([
      inputRow,
      this.uploadedPreviewDiv,
    ], 'moltrack-right-panel');

    const container = ui.divH([leftPanel, rightPanel], 'moltrack-container');
    this.view.root.append(container);

    this.preloadDefaultTable();
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

  private async preloadDefaultTable() {
    try {
      this.uploadedDf = await grok.data.loadTable(_package.webRoot + 'files/samples/compounds.csv');
      await grok.data.detectSemanticTypes(this.uploadedDf);

      ui.empty(this.uploadedPreviewDiv);
      this.uploadedPreviewDiv.append(this.uploadedDf.plot.grid().root);
    } catch (e: any) {
      grok.shell.error(`Failed to load dataset: ${e.message}`);
    }
  }

  private async registerEntities() {
    if (!this.uploadedDf) return;

    const csv = await this.uploadedDf.toCsv();
    const csvFile = DG.FileInfo.fromString('data.csv', csv);
    try {
      const df: DG.DataFrame = await grok.functions.call('Moltrack:registerBulk', {
        csvFile: csvFile,
        scope: ScopeLabels[this.entityTypeInput?.value],
        mapping: Object.keys(this.mappingDict).length === 0 ? '' : JSON.stringify(this.mappingDict),
        errorHandling: ErrorHandlingLabels[this.errorStrategyInput?.value],
      });

      this.uploadedDf.join(df, ['smiles'], ['smiles'], null, ['registration_status', 'registration_error_message'], DG.JOIN_TYPE.INNER, true);

      grok.shell.info('Entities registered successfully!');
    } catch (err: any) {
      grok.shell.error(`Registration failed: ${err.message}`);
    }
  }


  private async createMapping() {
    const parsedProps = JSON.parse(await fetchCompoundProperties());
    const targetProperties: TargetProperty[] = [
      ...parsedProps.map((p: any) => ({
        name: p.friendly_name ?? p.name,
        required: false,
      })),
      { name: 'smiles', required: true },
    ];

    await MolTrackDockerService.init();
    const autoMapping: Map<string, string> = await MolTrackDockerService.getAutoMapping(this.uploadedDf!.columns.names(), 'COMPOUND');

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
    });
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
