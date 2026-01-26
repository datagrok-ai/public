import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PropertyDesirability, DesirabilityMode} from '../mpo';
import {MpoDesirabilityLineEditor} from '../editors/mpo-line-editor';

const DESIRABILITY_MODES = ['freeform', 'gaussian', 'sigmoid', 'categorical'];

export class DesirabilityModeDialog {
  constructor(
    private readonly propertyName: string,
    private readonly prop: PropertyDesirability,
    private readonly onUpdate: (patch: Partial<PropertyDesirability>) => void,
  ) {}

  show(): void {
    const dialog = ui.dialog({
      title: `Desirability Settings: ${this.propertyName}`,
    });

    const previewEditor = new MpoDesirabilityLineEditor(this.prop, 300, 80);

    const modeInput = ui.input.choice('Mode', {
      items: DESIRABILITY_MODES,
      value: this.prop.mode ?? 'freeform',
      onValueChanged: (v) => {
        this.prop.mode = v as DesirabilityMode;
        this.onUpdate({mode: v as DesirabilityMode});
        updateParams();
        previewEditor.redrawAll();
      },
    });

    const min = this.float('Min', this.prop.min ?? previewEditor.getMinX(), (v) => {
      this.prop.min = v;
      this.onUpdate({min: v});
      previewEditor.redrawAll();
    });

    const max = this.float('Max', this.prop.max ?? previewEditor.getMaxX(), (v) => {
      this.prop.max = v;
      this.onUpdate({max: v});
      previewEditor.redrawAll();
    });

    const mean = this.float('Mean', this.prop.mean ?? 0, (v) => {
      this.prop.mean = v;
      this.onUpdate({mean: v});
      previewEditor.redrawAll();
    });

    const sigma = this.float('Sigma', this.prop.sigma ?? 1, (v) => {
      this.prop.sigma = Math.max(0.01, v);
      this.onUpdate({sigma: this.prop.sigma});
      previewEditor.redrawAll();
    });

    const x0 = this.float('x0', this.prop.x0 ?? 0, (v) => {
      this.prop.x0 = v;
      this.onUpdate({x0: v});
      previewEditor.redrawAll();
    });

    const k = this.float('k', this.prop.k ?? 10, (v) => {
      this.prop.k = Math.max(0.1, v);
      this.onUpdate({k: this.prop.k});
      previewEditor.redrawAll();
    });

    const paramPanel = ui.divV([]);

    const updateParams = () => {
      ui.empty(paramPanel);

      const inputs: DG.InputBase[] = [min, max];

      if (this.prop.mode === 'gaussian')
        inputs.push(mean, sigma);
      else
        inputs.push(x0, k);

      paramPanel.append(ui.h3('Parameters'));
      paramPanel.append(ui.form(inputs));
    };

    previewEditor.onParamsChanged = (p) => {
      Object.assign(this.prop, p);

      mean.value = this.prop.mean ?? 0;
      sigma.value = this.prop.sigma ?? 1;
      x0.value = this.prop.x0 ?? 0;
      k.value = this.prop.k ?? 10;
      min.value = this.prop.min ?? previewEditor.getMinX();
      max.value = this.prop.max ?? previewEditor.getMaxX();

      this.onUpdate(p);
    };

    previewEditor.onChanged.subscribe((line) => {
      this.prop.line = line;
      this.onUpdate({line});
    });

    updateParams();

    dialog.add(ui.divV([
      modeInput.root,
      paramPanel,
      previewEditor.root,
    ]));

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
