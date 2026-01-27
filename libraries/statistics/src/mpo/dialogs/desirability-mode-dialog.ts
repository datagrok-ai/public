import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PropertyDesirability, DesirabilityMode} from '../mpo';
import {MpoDesirabilityLineEditor} from '../editors/mpo-line-editor';

const DESIRABILITY_MODES = ['freeform', 'gaussian', 'sigmoid'];

export class DesirabilityModeDialog {
  constructor(
    private readonly propertyName: string,
    private readonly prop: PropertyDesirability,
    private readonly onUpdate: (patch: Partial<PropertyDesirability>) => void,
  ) {}

  show(): void {
    const original = {...this.prop};
    const working = {...this.prop};

    const dialog = ui.dialog({
      title: `Desirability Settings: ${this.propertyName}`,
    });

    const previewEditor = new MpoDesirabilityLineEditor(working, 300, 80);

    const modeInput = ui.input.choice('Mode', {
      items: DESIRABILITY_MODES,
      value: working.mode ?? 'freeform',
      onValueChanged: (v) => {
        working.mode = v as DesirabilityMode;
        updateParams();
        previewEditor.redrawAll();
      },
    });

    const min = this.float('Min', working.min ?? previewEditor.getMinX(), (v) => {
      working.min = v;
      previewEditor.redrawAll();
    });

    const max = this.float('Max', working.max ?? previewEditor.getMaxX(), (v) => {
      working.max = v;
      previewEditor.redrawAll();
    });

    const mean = this.float('Mean', working.mean ?? 0, (v) => {
      working.mean = v;
      previewEditor.redrawAll();
    });

    const sigma = this.float('Sigma', working.sigma ?? 1, (v) => {
      working.sigma = Math.max(0.01, v);
      previewEditor.redrawAll();
    });

    const x0 = this.float('x0', working.x0 ?? 0, (v) => {
      working.x0 = v;
      previewEditor.redrawAll();
    });

    const k = this.float('k', working.k ?? 10, (v) => {
      working.k = Math.max(0.1, v);
      previewEditor.redrawAll();
    });

    const paramPanel = ui.divV([]);

    const updateParams = () => {
      ui.empty(paramPanel);

      const inputs: DG.InputBase[] = [min, max];

      if (working.mode === 'gaussian')
        inputs.push(mean, sigma);
      if (working.mode === 'sigmoid')
        inputs.push(x0, k);

      paramPanel.append(ui.h3('Parameters'));
      paramPanel.append(ui.form(inputs));
    };

    previewEditor.onParamsChanged = (p) => {
      Object.assign(working, p);

      mean.value = working.mean ?? previewEditor.getDefaultMean();
      sigma.value = working.sigma ?? previewEditor.getDefaultSigma();
      x0.value = working.x0 ?? previewEditor.getDefaultX0();
      k.value = working.k ?? previewEditor.getDefaultK();
      min.value = working.min ?? previewEditor.getMinX();
      max.value = working.max ?? previewEditor.getMaxX();
    };

    previewEditor.onChanged.subscribe((line) => {
      working.line = line;
    });

    updateParams();

    dialog.add(ui.divV([
      modeInput.root,
      paramPanel,
      previewEditor.root,
    ]));

    dialog.onOK(() => {
      Object.assign(this.prop, working);
      this.onUpdate(working);
      dialog.close();
    });

    dialog.onCancel(() => {
      Object.assign(this.prop, original);
      dialog.close();
    });

    dialog.show();
  }

  private float(
    label: string,
    value: number,
    onChange: (v: number) => void,
  ): DG.InputBase {
    return ui.input.float(label, {
      value,
      format: '#0.000',
      onValueChanged: (v) => onChange(v ?? value),
    });
  }
}
