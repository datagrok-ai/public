import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ErrorHandling, MolTrackProp, Scope } from './constants';
import { registerBulk } from '../package';

let openedView: DG.ViewBase | null = null;

export abstract class EntityBaseView {
  view: DG.View;
  gridDiv: HTMLDivElement;
  sketcherInstance: grok.chem.Sketcher;
  previewDf: DG.DataFrame | undefined;

  protected props: DG.Property[] = [];
  protected formBackingObject: Record<string, any> = {};

  protected abstract scope: Scope;

  constructor(initialSmiles = 'CC(=O)Oc1ccccc1C(=O)O') {
    this.view = DG.View.create();
    this.gridDiv = ui.div('', 'moltrack-register-res-div');

    this.sketcherInstance = new grok.chem.Sketcher();
    DG.chem.currentSketcherType = 'Ketcher';
    this.sketcherInstance.setSmiles(initialSmiles);
    this.sketcherInstance.onChanged.subscribe(() => {
      this.previewDf?.col('smiles')?.set(0, this.sketcherInstance.getSmiles());
    });

    this.buildUI();
  }

  protected convertToDGProperty(p: MolTrackProp): DG.Property {
    const options: any = {
      name: p.name,
      type: p.value_type,
      description: p.description ?? '',
    };

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

  protected abstract getMolTrackDGProperties(): Promise<DG.Property[]>;

  private async registerButtonHandler() {
    try {
      ui.setUpdateIndicator(this.gridDiv, true);

      const smiles = this.sketcherInstance.getSmiles();
      const filteredProps = this.props.filter((p) => {
        const value = p.get(this.formBackingObject);
        return value !== null && value !== undefined && value !== '';
      });

      const propValues = filteredProps.reduce((acc, p) => {
        const value = p.get(this.formBackingObject);
        acc[p.name] = value;
        return acc;
      }, {} as Record<string, any>);

      const headers = ['smiles', ...Object.keys(propValues)];
      const row = [smiles, ...Object.values(propValues)];

      const csvString = headers.join(',') + '\n' + row.join(',') + '\n';
      const file = DG.FileInfo.fromString(`${this.scope}.csv`, csvString);

      const resultDf = await registerBulk(file, this.scope, '', ErrorHandling.REJECT_ROW);
      await grok.data.detectSemanticTypes(resultDf);

      ui.empty(this.gridDiv);
      this.gridDiv.appendChild(resultDf.plot.grid().root);
    } catch (err) {
      grok.shell.error(`Registration failed: ${err}`);
      return null;
    }
    ui.setUpdateIndicator(this.gridDiv, false);
  }

  private async buildUI() {
    this.props = await this.getMolTrackDGProperties();

    this.formBackingObject = Object.fromEntries(this.props.map((p) => [p.name, null]));
    const form = ui.input.form(this.formBackingObject, this.props, {
      onValueChanged: (value: any, input: DG.InputBase) => {
        this.previewDf?.col(input.property.name)?.set(0, value);
      },
    });

    const registerButton = ui.bigButton('REGISTER', async () => {
      await this.registerButtonHandler();
    });
    registerButton.classList.add('moltrack-run-register-button');

    this.view.setRibbonPanels([[registerButton]]);

    const topContainer = ui.divH([this.sketcherInstance.root, form], 'moltrack-top-container');
    this.sketcherInstance.root.classList.add('moltrack-sketcher-root');
    form.classList.add('moltrack-form');

    this.view.root.append(ui.divV([topContainer, this.gridDiv]));

    await this.createPreview();
  }

  private async createPreview() {
    const smiles = this.sketcherInstance.getSmiles();

    const rowObj: Record<string, any> = { smiles };
    this.props.forEach((p) => {
      rowObj[p.name] = this.formBackingObject[p.name];
    });

    this.previewDf = DG.DataFrame.fromObjects([rowObj]);

    if (this.previewDf) {
      await grok.data.detectSemanticTypes(this.previewDf);
      ui.empty(this.gridDiv);
      this.gridDiv.appendChild(this.previewDf.plot.grid().root);
    }
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
