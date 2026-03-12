/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subscription} from 'rxjs';
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

function sectionHeader(text: string): HTMLElement {
  return ui.divText(text, 'statistics-mpo-section-header');
}

export class DesirabilityModeDialog {
  constructor(
    private propertyName: string,
    private prop: PropertyDesirability,
    private onUpdate: (patch: Partial<PropertyDesirability>) => void,
    private onTypeChanged?: (newProp: PropertyDesirability) => void,
    private mappedCol?: DG.Column | null,
  ) {}

  private static strategyToLabel(strategy?: string): string {
    switch (strategy) {
    case 'default': return 'Use default score';
    case 'skip': return 'Skip property';
    default: return 'Exclude row';
    }
  }

  private buildDefaultScoreInput(prop: PropertyDesirability): HTMLElement {
    const mv = prop.missingValues;
    const scoreInput = ui.input.float('Default score', {value: mv?.strategy === 'default' ? mv.score : 0, min: 0, max: 1, format: '#0.000',
      onValueChanged: (v) => {
        prop.missingValues = {strategy: 'default', score: v};
        this.onUpdate({missingValues: prop.missingValues});
      },
    });
    scoreInput.root.style.display = mv?.strategy === 'default' ? '' : 'none';

    const choiceInput = ui.input.choice('If missing', {
      items: ['Exclude row', 'Use default score', 'Skip property'],
      value: DesirabilityModeDialog.strategyToLabel(mv?.strategy),
      onValueChanged: (v) => {
        const use = v === 'Use default score';
        if (v === 'Skip property')
          prop.missingValues = {strategy: 'skip'};
        else if (use)
          prop.missingValues = {strategy: 'default', score: scoreInput.value ?? 0};
        else
          prop.missingValues = {strategy: 'exclude'};
        scoreInput.root.style.display = use ? '' : 'none';
        this.onUpdate({missingValues: prop.missingValues} as any);
      },
    });
    choiceInput.setTooltip(
      'How to handle missing values:<br>' +
      '• <b>Exclude row</b> — assign an MPO score of 0 to the row<br>' +
      '• <b>Use default score</b> — use a fixed fallback desirability value<br>' +
      '• <b>Skip property</b> — calculate the score from the remaining properties',
    );
    scoreInput.setTooltip('Desirability score (0–1) to use as fallback.');
    return ui.form([choiceInput, scoreInput]);
  }

  show(): void {
    const original = structuredClone(this.prop);

    const dialog = ui.dialog({
      title: this.propertyName,
    });
    dialog.root.classList.add('statistics-mpo-desirability-dialog');

    const contentPanel = ui.divV([]);
    const cached: Record<string, PropertyDesirability> = {[original.functionType]: structuredClone(original)};
    const subs: Subscription[] = [];

    let modeInputRoot: HTMLElement | undefined;

    const typeInput = ui.input.choice('Type', {items: [...PROPERTY_TYPES], value: this.prop.functionType, onValueChanged: (v) => {
      if (v === this.prop.functionType)
        return;

      cached[this.prop.functionType] = structuredClone(this.prop);

      if (cached[v])
        this.prop = structuredClone(cached[v]);
      else if (v === 'categorical') {
        const cats = this.mappedCol?.isCategorical ?
          this.mappedCol.categories.map((c: string) => ({name: c, desirability: 1})) : undefined;
        this.prop = createDefaultCategorical(original.weight, cats);
      }
      else
        this.prop = createDefaultNumerical(original.weight);

      buildContent();
    }});

    const buildNumericalContent = () => {
      const prop = this.prop as NumericalDesirability;
      prop.mode ??= 'freeform';

      const previewEditor = new MpoDesirabilityLineEditor(prop, 355, 103);
      if (this.mappedCol?.isNumerical)
        previewEditor.setColumn(this.mappedCol);

      const modeInput = ui.input.choice('Mode', {items: DESIRABILITY_MODES, value: prop.mode, onValueChanged: (v) => {
        prop.mode = v as DesirabilityMode;
        this.onUpdate({mode: prop.mode} as any);
        updateParams();
        previewEditor.redrawAll(false);
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
          previewEditor.redrawAll(false);
        }}));
      }

      const syncInputs = () => {
        for (const cfg of configs)
          inputs.get(cfg.key)!.value = prop[cfg.key] ?? cfg.fallback();
      };

      inputs.get('min')!.setTooltip('Minimum property value');
      inputs.get('max')!.setTooltip('Maximum property value');

      const paramForm = ui.form([
        inputs.get('min')!, inputs.get('max')!,
        inputs.get('mean')!, inputs.get('sigma')!,
        inputs.get('x0')!, inputs.get('k')!,
      ]);
      paramForm.classList.add('statistics-mpo-param-grid');

      const updateParams = () => {
        inputs.get('mean')!.root.classList.toggle('statistics-mpo-hidden', prop.mode !== 'gaussian');
        inputs.get('sigma')!.root.classList.toggle('statistics-mpo-hidden', prop.mode !== 'gaussian');
        inputs.get('x0')!.root.classList.toggle('statistics-mpo-hidden', prop.mode !== 'sigmoid');
        inputs.get('k')!.root.classList.toggle('statistics-mpo-hidden', prop.mode !== 'sigmoid');
      };

      updateParams();

      previewEditor.onParamsChanged = (p) => {
        Object.assign(prop, p);
        syncInputs();
        this.onUpdate(p as any);
      };

      subs.push(previewEditor.onChanged.subscribe((line) => {
        prop.line = line;
        if (prop.mode === 'freeform')
          this.onUpdate({line} as any);
      }));

      modeInputRoot = modeInput.root;
      previewEditor.root.classList.add('statistics-mpo-plot');
      contentPanel.append(previewEditor.root, sectionHeader('PARAMETERS'), paramForm);
    };

    const buildCategoricalContent = () => {
      const prop = this.prop as CategoricalDesirability;
      const catEditor = new MpoCategoricalEditor(prop, true, true);
      if (this.mappedCol?.isCategorical)
        catEditor.setChoices([...this.mappedCol.categories]);
      subs.push(catEditor.onChanged.subscribe(() => this.onUpdate(this.prop)));
      contentPanel.append(sectionHeader('CATEGORIES'), catEditor.root);
    };

    const dispose = () => {
      for (const s of subs)
        s.unsubscribe();
      subs.length = 0;
    };

    const buildContent = () => {
      dispose();
      ui.empty(contentPanel);
      modeInputRoot = undefined;

      if (isNumerical(this.prop))
        buildNumericalContent();
      else
        buildCategoricalContent();

      const headerItems: HTMLElement[] = [];
      if (!this.mappedCol)
        headerItems.push(typeInput.root);
      if (modeInputRoot)
        headerItems.push(modeInputRoot);
      if (headerItems.length > 0)
        contentPanel.prepend(ui.divH(headerItems, 'statistics-mpo-dialog-header-row'));

      contentPanel.append(sectionHeader('MISSING VALUES'), this.buildDefaultScoreInput(this.prop));
    };

    buildContent();

    dialog.add(contentPanel);

    dialog.onOK(() => {
      dispose();
      if (this.prop.functionType !== original.functionType)
        this.onTypeChanged?.(this.prop);
      else
        this.onUpdate(this.prop);
      dialog.close();
    });

    dialog.onCancel(() => {
      dispose();
      this.prop = original;
      this.onUpdate(original);
      dialog.close();
    });

    dialog.show();
  }
}
