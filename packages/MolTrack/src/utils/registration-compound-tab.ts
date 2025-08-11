import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { fetchMolTrackProperties } from '../package';

let openedView: DG.ViewBase | null = null;

export class RegistrationCompoundView {
  view: DG.View;
  gridDiv: HTMLDivElement;
  sketcherInstance: grok.chem.Sketcher;

  private props: DG.Property[] = [];
  private formBackingObject: Record<string, any> = {};

  constructor(initialSmiles = 'CC(=O)Oc1ccccc1C(=O)O') {
    this.view = DG.View.create();
    this.view.name = 'Register a compound';
    this.gridDiv = ui.div('', 'moltrack-register-res-div');

    this.sketcherInstance = new grok.chem.Sketcher();
    DG.chem.currentSketcherType = 'Ketcher';
    this.sketcherInstance.setSmiles(initialSmiles);

    this.buildUI();
  }

  private async getMolTrackDGProperties(): Promise<DG.Property[]> {
    try {
      type MolTrackProp = {
        name: string;
        value_type: string;
        description?: string;
        pattern?: string;
      };

      const propsJson: MolTrackProp[] = (JSON.parse(await fetchMolTrackProperties()) as MolTrackProp[])
        .filter((p) => p.name.toLowerCase() !== 'corporate_compound_id');

      return propsJson.map((p) => {
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
      });
    } catch (err) {
      grok.shell.error(`Failed to fetch properties: ${err}`);
      return [];
    }
  }

  private async buildUI() {
    this.props = await this.getMolTrackDGProperties();

    this.formBackingObject = Object.fromEntries(this.props.map((p) => [p.name, null]));
    const form = ui.input.form(this.formBackingObject, this.props);

    const registerButton = ui.bigButton('REGISTER', async () => {
      const results = this.registerCompound();
      this.gridDiv.innerText = JSON.stringify(results, null, 2);
    });
    registerButton.classList.add('moltrack-run-register-button');

    this.view.setRibbonPanels([[registerButton]]);

    const topContainer = ui.divH([this.sketcherInstance.root, form], 'moltrack-top-container');
    topContainer.style.gap = '20px';
    topContainer.style.alignItems = 'flex-start';

    this.sketcherInstance.root.style.flex = '0 0 auto';
    form.style.flex = '1 1 auto';
    form.style.minWidth = '200px';

    this.view.root.append(ui.divV([topContainer/*, this.gridDiv*/], 'moltrack-register-main'));
  }

  private registerCompound() {
    const results = this.props.reduce((acc, p) => {
      const value = p.get(this.formBackingObject);
      if (value != null && value !== '') acc.push({ property: p.name, value });
      return acc;
    }, [] as { property: string; value: any }[]);
    return results;
  }

  show() {
    openedView?.close();
    openedView = this.view;
    grok.shell.addPreview(this.view);
  }
}
