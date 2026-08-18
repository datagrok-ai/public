import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {EnumeratorConfig} from './config';
import {PerRoundOverride} from './enumerate';
import {
  CHANGED_DOT_STYLE, clampRounds, combinationLimitsChanged, DataKey, estimateProductCount, MAX_ROUNDS, Mode,
  MODE_LABEL, panelHeader, productFiltersChangedCount, roundsLabel, tabPanel,
} from './shared';

export interface StrategySummaryDeps {
  getConfig: () => EnumeratorConfig;
  currentMode: () => Mode;
  currentRounds: () => number;
  templatesInput: DG.InputBase<DG.DataFrame | null>;
  bbsInput: DG.InputBase<DG.DataFrame | null>;
  reagentsInput: DG.InputBase<DG.DataFrame | null>;
  exclusionInput: DG.InputBase<DG.DataFrame | null>;
  buildPerRoundOverrides: (cfg: EnumeratorConfig) => PerRoundOverride[] | undefined;
  overrideCountFor: (overrides: PerRoundOverride[] | undefined, mode: Mode, r: number, key: DataKey) => number | null;
}

/** Right-tab recap for the "How to combine" section: mode, round-by-round chain per component, and
 * the product-count estimate with its caveats. */
export class StrategySummary {
  readonly panel: HTMLElement;
  private readonly host: HTMLElement;

  constructor(private readonly deps: StrategySummaryDeps) {
    this.host = ui.div([], {style: {padding: '16px'}});
    this.panel = tabPanel(
      panelHeader('How the current strategy and round count combine the reaction templates and building blocks.'),
      this.host, true);
  }

  render(): void {
    this.host.innerHTML = '';
    const config = this.deps.getConfig();
    const combChanged = combinationLimitsChanged(config);
    const prodChangedCount = productFiltersChangedCount(config);
    const tDf = this.deps.templatesInput.value;
    const bDf = this.deps.bbsInput.value;
    const mode = this.deps.currentMode();
    // Raw for display (so the number matches what's typed), clamped for the per-round row loop.
    const rounds = this.deps.currentRounds();
    const displayRounds = clampRounds(rounds);
    const n = estimateProductCount(tDf, bDf);

    const overrides = this.deps.buildPerRoundOverrides(config);
    const overrideCount = (r: number, key: DataKey): number | null =>
      this.deps.overrideCountFor(overrides, mode, r, key);

    const card = ui.div([], {style: {maxWidth: '480px'}});
    card.appendChild(ui.divText(`${MODE_LABEL[mode]} · ${roundsLabel(rounds)}`,
      {style: {fontWeight: 'bold', fontSize: '13px', marginBottom: '10px'}}));
    if (rounds > MAX_ROUNDS) {
      card.appendChild(ui.divText(
        `Showing the first ${MAX_ROUNDS} rounds — Number of rounds is capped at ${MAX_ROUNDS}.`,
        {style: {fontSize: '11px', color: 'var(--grey-5)', marginBottom: '8px'}}));
    }

    if (tDf && bDf) {
      const componentSection = (
        title: string, total: number, key: DataKey,
      ): HTMLElement => {
        const section = ui.div([], {style: {marginTop: '10px'}});
        section.appendChild(ui.divText(title,
          {style: {fontWeight: 'bold', fontSize: '12px', marginBottom: '4px'}}));
        for (let r = 1; r <= displayRounds; r++) {
          const oc = overrideCount(r, key);
          const rowChildren: HTMLElement[] = [
            ui.divText(`Round ${r}`, {style: {color: 'var(--grey-6)', width: '64px', flex: '0 0 auto'}}),
            ui.divText(oc != null ? `${oc} of ${total} (custom subset)` : `all ${total}`,
              oc != null ? {style: {fontWeight: '600'}} : undefined),
          ];
          if (oc != null)
            rowChildren.push(ui.div([], {style: {...CHANGED_DOT_STYLE, flex: '0 0 auto'}}));
          section.appendChild(ui.divH(rowChildren, {style: {gap: '8px', alignItems: 'center', padding: '2px 0'}}));
        }
        return section;
      };

      card.appendChild(componentSection('Reactions', tDf.rowCount, 'templates'));
      card.appendChild(componentSection('Building blocks', bDf.rowCount, 'buildingBlocks'));

      const rDf = this.deps.reagentsInput.value;
      if (mode === 'reagents' && rDf)
        card.appendChild(componentSection('Reagents', rDf.rowCount, 'reagents'));
    } else {
      card.appendChild(ui.divText('Pick reaction templates and building blocks to see round-by-round details.',
        {style: {color: 'var(--grey-5)', fontSize: '12px', marginTop: '4px'}}));
    }

    if (n > 0) {
      card.appendChild(ui.divText(`≈ ${n.toLocaleString('en-US')} estimated products`,
        {style: {marginTop: '12px', fontWeight: 'bold', fontSize: '13px', color: 'var(--blue-2)'}}));

      // The estimate is a naive multiplication — flag anything that would shrink the real output.
      const changedFilters = (combChanged ? 1 : 0) + prodChangedCount;
      const xDf = this.deps.exclusionInput.value;
      const hasExclusion = !!xDf && xDf.rowCount > 0;
      if (changedFilters > 0 || hasExclusion) {
        const bits: string[] = [];
        if (changedFilters > 0)
          bits.push(`${changedFilters} limit${changedFilters > 1 ? 's' : ''} changed from defaults`);
        if (hasExclusion) bits.push('exclusion substructures active');
        const caveatEl = ui.divText(`${bits.join(', ')} — actual output may be lower than this estimate.`,
          {style: {marginTop: '4px', fontSize: '11px', color: 'var(--grey-5)'}});
        if (changedFilters > 0) {
          // inline-block so the dot flows with the text instead of centering on the wrapped block.
          const dot = ui.span([], {style: {...CHANGED_DOT_STYLE, display: 'inline-block', marginRight: '6px'}});
          caveatEl.prepend(dot);
        }
        card.appendChild(caveatEl);
      }
    }

    this.host.appendChild(card);
  }
}
