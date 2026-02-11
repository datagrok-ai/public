import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {DesirabilityProfile} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MPO_SCORE_CHANGED_EVENT} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

export class FormulaEditor {
  readonly root: HTMLDivElement = ui.div([], 'chem-mpo-d-none');
  private formulaColumn: DG.Column | null = null;
  private observer: MutationObserver | null = null;

  get visible(): boolean { return !this.root.classList.contains('chem-mpo-d-none'); }
  set visible(v: boolean) { this.root.classList.toggle('chem-mpo-d-none', !v); }

  getExpression(): string | undefined {
    if (!this.visible)
      return undefined;
    return this.formulaColumn?.getTag(DG.Tags.Formula) ?? '';
  }

  async build(profile: DesirabilityProfile, df: DG.DataFrame): Promise<void> {
    this.disconnect();
    const properties = profile.properties;
    const propNames = Object.keys(properties);

    const resultCol = df.columns.addNewFloat('~formula~result');
    const defaultFormula = this.buildDefaultFormula(propNames, properties);
    resultCol.setTag(DG.Tags.Formula, defaultFormula);
    this.formulaColumn = resultCol;

    ui.empty(this.root);
    try {
      const widget = await grok.functions.call('PowerPack:formulaWidget', {col: resultCol}) as DG.Widget;
      this.root.append(widget.root);
      this.setupWidget(widget.root);
    }
    catch (e) {
      this.root.append(ui.divText('Formula editor is not available. Install PowerPack.'));
    }
    finally {
      df.columns.remove(resultCol.name);
    }
  }

  private buildDefaultFormula(
    propNames: string[],
    properties: DesirabilityProfile['properties'],
  ): string {
    const terms = propNames.map((n) => '${' + n + '} * ' + properties[n].weight);
    const totalWeight = propNames.reduce((sum, n) => sum + properties[n].weight, 0);
    return '(' + terms.join(' + ') + ') / ' + totalWeight;
  }

  private setupWidget(widgetRoot: HTMLElement): void {
    let retries = 0;
    const trySetup = () => {
      const cmContent = widgetRoot.querySelector('.cm-content');
      if (!cmContent) {
        if (++retries < 100)
          requestAnimationFrame(trySetup);
        return;
      }

      for (const btn of Array.from(widgetRoot.querySelectorAll('button'))) {
        if (btn.textContent === 'Edit in dialog' || btn.textContent === 'Apply') {
          btn.parentElement?.remove();
          break;
        }
      }

      let lastText = (cmContent as HTMLElement).textContent ?? '';
      let debounceTimer: ReturnType<typeof setTimeout> | null = null;
      this.observer = new MutationObserver(() => {
        if (debounceTimer)
          clearTimeout(debounceTimer);
        debounceTimer = setTimeout(() => {
          const text = (cmContent as HTMLElement).textContent ?? '';
          if (text === lastText)
            return;
          lastText = text;
          this.formulaColumn?.setTag(DG.Tags.Formula, text);
          grok.events.fireCustomEvent(MPO_SCORE_CHANGED_EVENT, {});
        }, 500);
      });
      this.observer.observe(cmContent, {characterData: true, childList: true, subtree: true});
    };
    requestAnimationFrame(trySetup);
  }

  disconnect(): void {
    this.observer?.disconnect();
    this.observer = null;
  }
}
