/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ErrorHandling, MolTrackProp, Scope } from './constants';
import { createPath, fetchBatchProperties, fetchCompoundProperties, registerBulk } from '../package';
import RandExp from 'randexp';

let openedView: DG.ViewBase | null = null;

export abstract class EntityBaseView {
  view: DG.View;
  sketcherInstance: grok.chem.Sketcher;
  previewDf: DG.DataFrame | undefined;

  inputs: DG.InputBase[] = [];
  formBackingObject: Record<string, any> = {};

  title: string = 'Register new chemical compounds';
  initialSmiles: string = 'CC(=O)Oc1ccccc1C(=O)O';
  singleRetrieved: boolean = false;

  protected messageContainer: HTMLDivElement = ui.div([], 'moltrack-info-container');

  protected abstract scope: Scope;

  constructor(buildUI: boolean = true) {
    this.view = DG.View.create();
    this.sketcherInstance = new grok.chem.Sketcher();
    DG.chem.currentSketcherType = 'Ketcher';

    if (buildUI)
      this.buildUIMethod();
  }

  protected convertToDGProperty(p: MolTrackProp): DG.Property {
    const options: any = {
      name: p.name,
      type: p.value_type,
    };

    if (p.description)
      options.description = p.description;

    if (p.pattern) {
      const regex = new RegExp(p.pattern);
      options.valueValidators = [
        (val: any) =>
          val == null || val === '' ?
            null :
            regex.test(String(val)) ?
              null :
              `Value does not match pattern: ${p.pattern}`,
      ];
    }

    return DG.Property.fromOptions(options);
  }

  private async createPropertySection(
    title: string,
    fetchPropsFn: () => Promise<any>,
    options?: {
        disableNames?: string[];
        initiallyOpen?: boolean;
  },
  ): Promise<{ section: HTMLElement; inputs: DG.InputBase[]; formBackingObject: Record<string, any> }> {
    const { disableNames = [], initiallyOpen = false } = options ?? {};
    const disableAll = disableNames.includes('*');

    let props: DG.Property[] = [];
    let propArray: MolTrackProp[] = [];
    const formBackingObject: Record<string, any> = {};

    try {
      const rawProps = await fetchPropsFn();
      const parsed: any = JSON.parse(rawProps);

      propArray = parsed.properties ?? parsed;
      props = propArray.map(this.convertToDGProperty.bind(this));
      for (const p of props)
        formBackingObject[p.name] = null;
    } catch (err) {
      grok.shell.error(`Failed to fetch properties for "${title}": ${err}`);
    }

    const inputs = props.map((p) => {
      const input = DG.InputBase.forProperty(p, formBackingObject);
      const rawProp = propArray.find((rp) => rp.name === p.name);
      if (rawProp?.pattern)
        (input.input as HTMLInputElement).placeholder = `e.g., ${generateExample(rawProp?.pattern)}`;

      if (disableAll || disableNames.includes(p.name))
        input.enabled = false;
      return input;
    });

    const form = ui.wideForm(inputs);
    form.classList.add('moltrack-compound-form', 'moltrack-form');
    const section = createCollapsibleSection(title, form, initiallyOpen);

    return { section, inputs, formBackingObject };
  }

  private async registerButtonHandler() {
    try {
      const smiles = this.sketcherInstance.getSmiles();
      const propValues = this.collectNonEmptyInputValues();
      const csvFile = this.createCsvFile(smiles, propValues);

      const resultDf = await registerBulk(csvFile, this.scope, '', ErrorHandling.REJECT_ROW);

      ui.empty(this.messageContainer);
      if (!resultDf)
        return;

      const { status, compoundId, errorMsg } = this.extractResultData(resultDf);
      this.showRegistrationMessage(status, compoundId, errorMsg);

      if (status === 'success' && openedView)
        openedView.path = `${createPath('Compound')}?corporate_compound_id=${encodeURIComponent(compoundId)}`;
    } catch {
      ui.empty(this.messageContainer);
    }
  }

  private collectNonEmptyInputValues(): Record<string, any> {
    return this.inputs
      .filter((input) => input.value !== null && input.value !== undefined && input.value !== '')
      .reduce((acc, input) => {
        acc[input.property.name] = input.value;
        return acc;
      }, {} as Record<string, any>);
  }

  private createCsvFile(smiles: string, propValues: Record<string, any>): DG.FileInfo {
    const headers = ['smiles', ...Object.keys(propValues)];
    const row = [smiles, ...Object.values(propValues)];
    const csvString = `${headers.join(',')}\n${row.join(',')}\n`;
    return DG.FileInfo.fromString(`${this.scope}.csv`, csvString);
  }

  private extractResultData(df: DG.DataFrame): { status: string, compoundId: string, errorMsg: string } {
    const getValue = (colName: string, fallback = '') => df.col(colName)?.get(0) ?? fallback;

    return {
      status: getValue('registration_status', 'Unknown'),
      compoundId: getValue('property_corporate_compound_id', ''),
      errorMsg: getValue('registration_error_message', 'Unknown error'),
    };
  }

  private showRegistrationMessage(status: string, compoundId: string, errorMsg: string) {
    const header = status === 'success' ?
      'Compound successfully registered!' :
      'Compound registration failed!';
    const message = status === 'success' ?
      `The compound has been added to the database with ID: ${compoundId}` :
      `Error: ${errorMsg}`;

    const infoDiv = ui.info(message, header, true);
    const bar = infoDiv.querySelector('.grok-info-bar') as HTMLElement;
    if (bar)
      bar.style.setProperty('background-color', status === 'success' ? '#d4edda' : '#f8d7da', 'important');

    this.messageContainer.appendChild(infoDiv);
  }

  public async buildUIMethod() {
    const titleText = ui.divText(this.title, 'moltrack-title');
    this.messageContainer.append(titleText);
    this.sketcherInstance.setSmiles(this.initialSmiles);

    if (!this.singleRetrieved) {
      const registerButton = ui.bigButton('REGISTER', async () => {
        await this.registerButtonHandler();
      });
      registerButton.classList.add('moltrack-run-register-button');
      this.view.setRibbonPanels([[registerButton]]);
    }

    const {
      section: compoundSection,
      inputs: compoundInputs,
      formBackingObject: compoundFormBackingObject,
    } = await this.createPropertySection(
      'Compound properties',
      fetchCompoundProperties,
      {
        disableNames: this.singleRetrieved ? ['*'] : reservedProperties,
        initiallyOpen: true,
      },
    );

    const {
      section: batchSection,
      inputs: batchInputs,
      formBackingObject: batchFormBackingObject,
    } = await this.createPropertySection(
      'Batch properties',
      fetchBatchProperties,
      {
        disableNames: this.singleRetrieved ? ['*'] : reservedProperties,
        initiallyOpen: false,
      },
    );

    this.inputs = [...compoundInputs, ...batchInputs];
    this.formBackingObject = {
      ...compoundFormBackingObject,
      ...batchFormBackingObject,
    };

    const topSections = ui.divV([compoundSection, batchSection], 'moltrack-top-sections');
    const sketcherFormContainer = ui.divH(
      [this.sketcherInstance.root, topSections],
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

function createCollapsibleSection(
  title: string,
  form: HTMLElement,
  initiallyOpen: boolean = false,
): HTMLElement {
  const chevron = ui.iconFA(
    initiallyOpen ? 'chevron-down' : 'chevron-right',
    () => toggleSection(form, chevron, title),
    `Toggle ${title}`,
  );
  chevron.classList.add('moltrack-chevron');

  const headerText = ui.divText(title, 'moltrack-section-title');
  const header = ui.divH([chevron, headerText], 'moltrack-section-header');

  form.style.display = initiallyOpen ? 'block' : 'none';
  form.classList.add('moltrack-section-form');

  const section = ui.divV([header, form], 'moltrack-section');
  return section;
}

function toggleSection(form: HTMLElement, chevron: HTMLElement, title: string) {
  const isHidden = form.style.display === 'none';
  form.style.display = isHidden ? 'block' : 'none';
  chevron.className = `fa fa-chevron-${isHidden ? 'down' : 'right'}`;
}

function generateExample(pattern: string, index: number = 1): string {
  const template = pattern.match(/\{\:0?(\d+)d\}/);
  if (template) {
    const width = Number(template[1]);
    return pattern.replace(/\{\:0?(\d+)d\}/, index.toString().padStart(width, '0'));
  }

  try {
    return new RandExp(new RegExp(pattern)).gen();
  } catch {
    return '';
  }
}

const reservedProperties = ['corporate_compound_id', 'corporate_batch_id'];
