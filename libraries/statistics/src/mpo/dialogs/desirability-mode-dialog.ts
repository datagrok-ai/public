/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  PropertyDesirability, NumericalDesirability, CategoricalDesirability, DesirabilityMode,
  createDefaultCategorical, createDefaultNumerical, isNumerical,
} from '../mpo';
import {MpoDesirabilityLineEditor} from '../editors/mpo-line-editor';
import {MpoCategoricalEditor} from '../editors/mpo-categorical-editor';

const DESIRABILITY_MODES: DesirabilityMode[] = ['freeform', 'gaussian', 'sigmoid'];
const PROPERTY_TYPES = ['numerical', 'categorical'] as const;

type ParamConfig = {
  key: 'min' | 'max' | 'mean' | 'sigma' | 'x0' | 'k';
  label: string;
  fallback: () => number;
  transform?: (v: number) => number;
};

export class DesirabilityModeDialog {
  constructor(
    private propertyName: string,
    private prop: PropertyDesirability,
    private onUpdate: (patch: Partial<PropertyDesirability>) => void,
    private onTypeChanged?: (newProp: PropertyDesirability) => void,
  ) {}

  show(): void {
    const original = structuredClone(this.prop);

    const dialog = ui.dialog({
      title: `Desirability Settings: ${this.propertyName}`,
    });

    const contentPanel = ui.divV([]);
    const cached: Record<string, PropertyDesirability> = {[original.functionType]: structuredClone(original)};

    const typeInput = ui.input.choice('Type', {items: [...PROPERTY_TYPES], value: this.prop.functionType, onValueChanged: (v) => {
      if (v === this.prop.functionType)
        return;

      cached[this.prop.functionType] = structuredClone(this.prop);

      if (cached[v])
        this.prop = structuredClone(cached[v]);
      else if (v === 'categorical')
        this.prop = createDefaultCategorical(original.weight);
      else
        this.prop = createDefaultNumerical(original.weight);

      buildContent();
    }});

    const buildNumericalContent = () => {
      const prop = this.prop as NumericalDesirability;
      prop.mode ??= 'freeform';

      const previewEditor = new MpoDesirabilityLineEditor(prop, 300, 80);

      const modeInput = ui.input.choice('Mode', {items: DESIRABILITY_MODES, value: prop.mode, onValueChanged: (v) => {
        prop.mode = v as DesirabilityMode;
        this.onUpdate({mode: prop.mode} as any);
        updateParams();
        previewEditor.redrawAll();
      }});

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
        inputs.set(cfg.key, ui.input.float(cfg.label, {value: prop[cfg.key] ?? cfg.fallback(), format: '#0.000', onValueChanged: (v) => {
          const value = cfg.transform ? cfg.transform(v ?? cfg.fallback()) : (v ?? cfg.fallback());
          prop[cfg.key] = value;
          this.onUpdate({[cfg.key]: value} as any);
          previewEditor.redrawAll();
        }}));
      }

      const syncInputs = () => {
        for (const cfg of configs)
          inputs.get(cfg.key)!.value = prop[cfg.key] ?? cfg.fallback();
      };

      const paramPanel = ui.divV([]);

      const updateParams = () => {
        ui.empty(paramPanel);

        const form: DG.InputBase[] = [inputs.get('min')!, inputs.get('max')!];

        if (prop.mode === 'gaussian')
          form.push(inputs.get('mean')!, inputs.get('sigma')!);
        if (prop.mode === 'sigmoid')
          form.push(inputs.get('x0')!, inputs.get('k')!);

        paramPanel.append(ui.h3('Parameters'));
        paramPanel.append(ui.form(form));
      };

      previewEditor.onParamsChanged = (p) => {
        Object.assign(prop, p);
        syncInputs();
        this.onUpdate(p as any);
        previewEditor.redrawAll();
      };

      previewEditor.onChanged.subscribe((line) => {
        prop.line = line;
        this.onUpdate({line} as any);
      });

      updateParams();

      contentPanel.append(modeInput.root, paramPanel, previewEditor.root);
    };

    const buildCategoricalContent = () => {
      const prop = this.prop as CategoricalDesirability;
      const catEditor = new MpoCategoricalEditor(prop, true, true);
      catEditor.onChanged.subscribe(() => this.onUpdate(this.prop));
      contentPanel.append(catEditor.root);
    };

    const buildContent = () => {
      ui.empty(contentPanel);

      if (isNumerical(this.prop))
        buildNumericalContent();
      else
        buildCategoricalContent();
    };

    buildContent();

    dialog.add(ui.divV([typeInput.root, contentPanel]));

    dialog.onOK(() => {
      if (this.prop.functionType !== original.functionType)
        this.onTypeChanged?.(this.prop);
      else
        this.onUpdate(this.prop);
      dialog.close();
    });

    dialog.onCancel(() => {
      this.prop = original;
      this.onUpdate(original);
      dialog.close();
    });

    dialog.show();
  }
}
