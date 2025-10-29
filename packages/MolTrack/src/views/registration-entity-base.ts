/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { ErrorHandling, MolTrackProp, Scope } from '../utils/constants';
import { createPath, createPropertySection } from '../utils/view-utils';
import { buildPropertyOptions, getCorporateCompoundIdByExactStructure } from '../utils/utils';

import { fetchBatchProperties, fetchCompoundProperties, registerBulk } from '../package';

import RandExp from 'randexp';

let openedView: DG.ViewBase | null = null;

export class EntityBaseView {
  view: DG.View;
  sketcherInstance: grok.chem.Sketcher;

  inputs: DG.InputBase[] = [];
  batchInputs: DG.InputBase<any>[] = [];
  isBatchSectionExpanded: boolean = false;
  formBackingObject: Record<string, any> = {};

  title: string = 'Register a new compound';
  initialSmiles: string = 'CC(=O)Oc1ccccc1C(=O)O';
  singleRetrieved: boolean = false;
  compoundExists: boolean = false;
  invalidForm: boolean = false;
  path: string = 'Compound';

  protected messageContainer: HTMLDivElement = ui.div([], 'moltrack-info-container');
  registerButton: HTMLButtonElement | undefined;

  constructor(buildUI: boolean = true) {
    this.view = DG.View.create();
    this.sketcherInstance = new grok.chem.Sketcher();

    const validationFunc = (s: string) => {
      const valFunc = DG.Func.find({package: 'Chem', name: 'validateMolecule'})[0];
      if (valFunc) {
        const funcCall: DG.FuncCall = valFunc.prepare({s});
        funcCall.callSync();
        const res = funcCall.getOutputParamValue();
        return res;
      }
    };
    this.sketcherInstance._validationFunc = (s) => validationFunc(s);
    this.sketcherInstance.onChanged.subscribe(async () => {
      ui.empty(this.messageContainer);
      this.messageContainer.appendChild(ui.divText(this.title, 'moltrack-title'));

      const smiles = this.sketcherInstance.getSmiles();
      this.compoundExists = !!(await getCorporateCompoundIdByExactStructure(smiles));
      this.invalidForm = !smiles || smiles.trim().length === 0;
      this.registerButton?.classList.toggle('dim', (this.compoundExists || this.invalidForm));
    });
    DG.chem.currentSketcherType = 'Ketcher';

    if (buildUI)
      this.buildUIMethod();
  }

  protected convertToDGProperty(p: MolTrackProp): DG.Property {
    return DG.Property.fromOptions(
      buildPropertyOptions(p, {reserved: reservedProperties, skipReservedCheck: true}),
    );
  }

  private clearAll() {
    this.sketcherInstance.setValue('');
    for (const input of this.inputs)
      input.value = null;
  }

  private async registerButtonHandler() {
    try {
      const smiles = this.sketcherInstance.getSmiles();
      const propValues = this.collectNonEmptyInputValues();
      const batchInputNames = this.batchInputs.map((input) => input.property.name);
      const isBatch = Object.keys(propValues).find((name) => batchInputNames.includes(name)) !== undefined;
      const scope = (isBatch || this.isBatchSectionExpanded) ? Scope.BATCHES : Scope.COMPOUNDS;
      const singularScope = scope.replace(/(es|s)$/, '');
      const csvFile = this.createCsvFile(smiles, propValues, scope);

      const resultDf = await registerBulk(csvFile, scope, '', ErrorHandling.REJECT_ROW);

      ui.empty(this.messageContainer);
      if (!resultDf)
        return;

      const { status, compoundId, batchId, errorMsg } = this.extractResultData(resultDf);
      this.showRegistrationMessage(status, batchId?.trim() || compoundId, errorMsg, singularScope);

      for (const [propName, value] of [
        ['corporate_compound_id', compoundId],
        ['corporate_batch_id', batchId],
      ]) {
        if (value) {
          const input = this.inputs.find((inp) => inp.property.name === propName);
          if (input) input.value = value;
        }
      }

      if (openedView) {
        const basePath = createPath(this.path);
        openedView.path = status === 'success' ?
          `${basePath}?corporate_${singularScope}_id=${encodeURIComponent(batchId?.trim() || compoundId)}` :
          basePath;
      }

      if (status === 'success') {
        this.compoundExists = true;
        this.registerButton?.classList.toggle('dim', this.compoundExists);
      }
    } catch {
      ui.empty(this.messageContainer);
    }
  }

  private collectNonEmptyInputValues(): Record<string, any> {
    return this.inputs
      .filter((input) => input.value !== null && input.value !== undefined && input.value !== '')
      .reduce((acc, input) => {
        acc[input.property.name] = input.stringValue;
        return acc;
      }, {} as Record<string, any>);
  }

  private createCsvFile(smiles: string, propValues: Record<string, any>, scope: string): DG.FileInfo {
    const headers = ['smiles', ...Object.keys(propValues)];
    const row = [smiles, ...Object.values(propValues)];
    const csvString = `${headers.join(',')}\n${row.join(',')}\n`;
    return DG.FileInfo.fromString(`${scope}.csv`, csvString);
  }

  private extractResultData(df: DG.DataFrame): { status: string, compoundId: string, batchId: string, errorMsg: string } {
    const getValue = (colName: string, fallback = '') => df.col(colName)?.get(0) ?? fallback;

    return {
      status: getValue('registration_status', 'Unknown'),
      compoundId: getValue('corporate_compound_id', ''),
      batchId: getValue('corporate_batch_id', ''),
      errorMsg: getValue('registration_error_message', 'Unknown error'),
    };
  }

  private showRegistrationMessage(status: string, id: string, errorMsg: string, scope: string) {
    const upperCasedScope = scope.charAt(0).toUpperCase() + scope.slice(1);
    const header = status === 'success' ?
      `${id} successfully registered!` :
      `${upperCasedScope} registration failed!`;

    const errorMatch = errorMsg.match(/^\d+:\s*(.*)$/);
    const message = status === 'success' ? '' : `${errorMatch ? errorMatch[1] : ''}`;

    const infoDiv = ui.info(message, header, true);
    const bar = infoDiv.querySelector('.grok-info-bar') as HTMLElement;
    if (bar) {
      bar.classList.toggle('moltrack-bar-success', status === 'success');
      bar.classList.toggle('moltrack-bar-error', status !== 'success');
    }

    this.messageContainer.appendChild(infoDiv);
  }

  private async buildSketcherOrPreview(): Promise<HTMLElement> {
    const titleText = ui.divText(this.title, 'moltrack-title');
    this.messageContainer.append(titleText);

    this.sketcherInstance.setSmiles(this.initialSmiles);

    if (this.singleRetrieved) {
      const image: HTMLCanvasElement = ui.canvas();
      image.width = 600;
      image.height = 300;

      const container = ui.div();
      container.appendChild(image);

      await grok.chem.canvasMol(0, 0, image.width, image.height, image, this.initialSmiles);

      return container;
    } else
      return this.sketcherInstance.root;
  }


  public async buildUIMethod() {
    const container = await this.buildSketcherOrPreview();
    if (!this.singleRetrieved) {
      this.registerButton = ui.bigButton('REGISTER', () => {});

      this.registerButton!.addEventListener('click', async (e) => {
        if (this.compoundExists || this.invalidForm) {
          e.preventDefault();
          e.stopPropagation();
        } else
          await this.registerButtonHandler();
      });

      this.registerButton.onmouseenter = (e) => {
        let tooltipMessage: string | null = null;

        if (this.compoundExists) {
          tooltipMessage = 'This compound already exists.';
          this.registerButton!.style.cursor = 'not-allowed';
        } else if (this.invalidForm) {
          tooltipMessage = 'Please check the form. Values are invalid or do not match the required patterns.';
          this.registerButton!.style.cursor = 'not-allowed';
        } else
          this.registerButton!.style.cursor = 'pointer';

        if (tooltipMessage)
          ui.tooltip.show(tooltipMessage, e.clientX + 15, e.clientY);
      };

      this.registerButton.onmouseleave = (e) => ui.tooltip.hide();

      const clearAllIcon = ui.iconFA('eraser', () => this.clearAll(), 'Clear all');
      this.view.setRibbonPanels([[this.registerButton, clearAllIcon]]);
    }

    const {
      section: compoundSection,
      inputs: compoundInputs,
      formBackingObject: compoundFormBackingObject,
    } = await createPropertySection(
      'Compound properties',
      fetchCompoundProperties,
      this.convertToDGProperty.bind(this),
      {
        disableNames: this.singleRetrieved ? ['*'] : reservedProperties,
        initiallyOpen: true,
        reservedProperties: reservedProperties,
        generateExample: generateExample,
        onValidationChange: (invalid) => {
          this.invalidForm = invalid;
          this.registerButton?.classList.toggle('dim', this.invalidForm);
        },
      },
    );

    const {
      section: batchSection,
      inputs: batchInputs,
      formBackingObject: batchFormBackingObject,
    } = await createPropertySection(
      'Batch properties',
      fetchBatchProperties,
      this.convertToDGProperty.bind(this),
      {
        disableNames: this.singleRetrieved ? ['*'] : reservedProperties,
        initiallyOpen: this.isBatchSectionExpanded,
        reservedProperties: reservedProperties,
        generateExample: generateExample,
        onValidationChange: (invalid) => {
          this.invalidForm = invalid;
          this.registerButton?.classList.toggle('dim', this.invalidForm);
        },
      },
    );

    this.batchInputs = batchInputs;
    this.inputs = [...compoundInputs, ...batchInputs];
    this.formBackingObject = {
      ...compoundFormBackingObject,
      ...batchFormBackingObject,
    };

    const element = this.singleRetrieved ? container : this.sketcherInstance.root;
    const topSections = ui.divV([compoundSection, batchSection], 'moltrack-top-sections');
    const sketcherFormContainer = ui.divH(
      [element, topSections],
      'moltrack-top-container',
    );

    const registerContainer = ui.divV([this.messageContainer, sketcherFormContainer], 'moltrack-register-container');

    this.sketcherInstance.root.classList.add('moltrack-sketcher-root');
    this.view.root.append(registerContainer);
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}

function generateExample(pattern: string): string {
  try {
    return new RandExp(new RegExp(pattern)).gen();
  } catch {
    return '';
  }
}

const reservedProperties = ['corporate_compound_id', 'corporate_batch_id'];
