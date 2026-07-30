import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {EnumeratorConfig} from './config';
import {PerRoundOverride} from './enumerate';
import {
  CHANGED_DOT_STYLE, combinationLimitsChanged, DataKey, estimateProductCount, Mode, MODE_LABEL,
  OVERRIDE_DOT_COLOR, panelHeader, productFiltersChangedCount, roundsLabel, tabPanel,
} from './enumerator-app';

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

/** Right-tab recap for the "How to combine" section: mode, round-by-round chain per component,
 * and the product-count estimate with its caveats. Gives Strategy its own relevant content instead
 * of leaving whatever data grid was last shown; refreshed from refreshCfgRibbon(), same trigger as
 * the ribbon chips and accordion subtitles. */
export class StrategySummary {
  readonly panel: HTMLElement;
  private readonly host: HTMLElement;

  constructor(private readonly deps: StrategySummaryDeps) {
    this.host = ui.div([], {style: {padding: '16px'}});
    this.panel = tabPanel(
      panelHeader('How the current strategy and round count combine the reaction templates and building blocks.'),
      this.host, true);
  }

  // Optional params let refreshCfgRibbon() pass in already-computed values; falls back to
  // computing fresh otherwise.
  render(
    combChanged: boolean = combinationLimitsChanged(this.deps.getConfig()),
    prodChangedCount: number = productFiltersChangedCount(this.deps.getConfig()),
  ): void {
    this.host.innerHTML = '';
    const config = this.deps.getConfig();
    const tDf = this.deps.templatesInput.value;
    const bDf = this.deps.bbsInput.value;
    const mode = this.deps.currentMode();
    const rounds = this.deps.currentRounds();
    const n = estimateProductCount(tDf, bDf);

    // Per-round subset overrides, computed once for both the round diagram and the per-component
    // sections below.
    const overrides = this.deps.buildPerRoundOverrides(config);
    const overrideCount = (r: number, key: DataKey): number | null =>
      this.deps.overrideCountFor(overrides, mode, r, key);

    const card = ui.div([], {style: {maxWidth: '480px'}});
    card.appendChild(ui.divText(`${MODE_LABEL[mode]} · ${roundsLabel(rounds)}`,
      {style: {fontWeight: 'bold', fontSize: '13px', marginBottom: '10px'}}));

    if (tDf && bDf) {
      // One section per component, each listing what every round uses.
      const componentSection = (
        title: string, total: number, key: DataKey,
      ): HTMLElement => {
        const section = ui.div([], {style: {marginTop: '10px'}});
        section.appendChild(ui.divText(title,
          {style: {fontWeight: 'bold', fontSize: '12px', marginBottom: '4px'}}));
        for (let r = 1; r <= rounds; r++) {
          const oc = overrideCount(r, key);
          const rowChildren: HTMLElement[] = [
            ui.divText(`Round ${r}`, {style: {color: 'var(--grey-6)', width: '64px', flex: '0 0 auto'}}),
            ui.divText(oc != null ? `${oc} of ${total} (custom subset)` : `all ${total}`,
              oc != null ? {style: {fontWeight: '600'}} : undefined),
          ];
          if (oc != null) {
            rowChildren.push(ui.div([], {style: {width: '6px', height: '6px', borderRadius: '50%',
              background: OVERRIDE_DOT_COLOR, flex: '0 0 auto'}}));
          }
          section.appendChild(ui.divH(rowChildren, {style: {gap: '8px', alignItems: 'center', padding: '2px 0'}}));
        }
        return section;
      };

      card.appendChild(componentSection('Reactions', tDf.rowCount, 'templates'));
      card.appendChild(componentSection('Building blocks', bDf.rowCount, 'buildingBlocks'));

      // Reagents mode has a third data source just as central to the round math — show it too.
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

      // The estimate above is a naive multiplication — flag when active filters/limits (vs.
      // platform defaults) would actually shrink the real output.
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
          // inline-block so the dot flows with the text instead of centering against the wrapped block.
          const dot = ui.span([], {style: {...CHANGED_DOT_STYLE, display: 'inline-block', marginRight: '6px'}});
          caveatEl.prepend(dot);
        }
        card.appendChild(caveatEl);
      }
    }

    this.host.appendChild(card);
  }
}
