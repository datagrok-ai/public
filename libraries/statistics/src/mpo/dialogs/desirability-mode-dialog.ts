/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NumericalDesirability, DesirabilityMode} from '../mpo';
import {MpoDesirabilityLineEditor} from '../editors/mpo-line-editor';

const DESIRABILITY_MODES = ['freeform', 'gaussian', 'sigmoid'];

type ParamConfig = {
  key: 'min' | 'max' | 'mean' | 'sigma' | 'x0' | 'k';
  label: string;
  fallback: () => number;
  transform?: (v: number) => number;
};

export class DesirabilityModeDialog {
  constructor(
    private propertyName: string,
    private prop: NumericalDesirability,
    private onUpdate: (patch: Partial<NumericalDesirability>) => void,
  ) {}

  show(): void {
    const original = structuredClone(this.prop);

    const dialog = ui.dialog({
      title: `Desirability Settings: ${this.propertyName}`,
    });

    const previewEditor = new MpoDesirabilityLineEditor(this.prop, 300, 80);

    const modeInput = ui.input.choice('Mode', {
      items: DESIRABILITY_MODES,
      value: this.prop.mode ?? 'freeform',
      onValueChanged: (v) => {
        this.prop.mode = v as DesirabilityMode;
        this.onUpdate({mode: this.prop.mode});
        updateParams();
        previewEditor.redrawAll();
      },
    });

    const configs: ParamConfig[] = [
      {key: 'min', label: 'Min', fallback: () => previewEditor.getMinX()},
      {key: 'max', label: 'Max', fallback: () => previewEditor.getMaxX()},
      {key: 'mean', label: 'Mean', fallback: () => previewEditor.getDefaultMean()},
      {key: 'sigma', label: 'Sigma', fallback: () => previewEditor.getDefaultSigma(), transform: (v) => Math.max(0.01, v)},
      {key: 'x0', label: 'x0', fallback: () => previewEditor.getDefaultX0()},
      {key: 'k', label: 'k', fallback: () => previewEditor.getDefaultK(), transform: (v) => Math.max(0.1, v)},
    ];

    const inputs = new Map<string, DG.InputBase>();
    for (const cfg of configs) {
      inputs.set(cfg.key, ui.input.float(cfg.label, {
        value: this.prop[cfg.key] ?? cfg.fallback(),
        format: '#0.000',
        onValueChanged: (v) => {
          const value = cfg.transform ? cfg.transform(v ?? cfg.fallback()) : (v ?? cfg.fallback());
          this.prop[cfg.key] = value;
          this.onUpdate({[cfg.key]: value});
          previewEditor.redrawAll();
        },
      }));
    }

    const syncInputs = () => {
      for (const cfg of configs)
        inputs.get(cfg.key)!.value = this.prop[cfg.key] ?? cfg.fallback();
    };

    const paramPanel = ui.divV([]);

    const updateParams = () => {
      ui.empty(paramPanel);

      const form: DG.InputBase[] = [inputs.get('min')!, inputs.get('max')!];

      if (this.prop.mode === 'gaussian')
        form.push(inputs.get('mean')!, inputs.get('sigma')!);
      if (this.prop.mode === 'sigmoid')
        form.push(inputs.get('x0')!, inputs.get('k')!);

      paramPanel.append(ui.h3('Parameters'));
      paramPanel.append(ui.form(form));
    };

    previewEditor.onParamsChanged = (p) => {
      Object.assign(this.prop, p);
      syncInputs();
      this.onUpdate(p);
      previewEditor.redrawAll();
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

    dialog.onOK(() => {
      this.onUpdate(this.prop);
      dialog.close();
    });

    dialog.onCancel(() => {
      Object.assign(this.prop, original);
      syncInputs();
      previewEditor.redrawAll();
      dialog.close();
    });

    dialog.show();
  }
}
