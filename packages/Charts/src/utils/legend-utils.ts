import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';

import '../../css/legend.css';

export enum VISIBILITY_MODE {
  ALWAYS = 'Always',
  AUTO = 'Auto',
  NEVER = 'Never',
}

export type VisibilityMode = `${VISIBILITY_MODE}`;

export class LegendHelper {
  legendDiv: HTMLDivElement = ui.div();
  legend: DG.Legend | null = null;
  column: DG.Column | null = null;
  onCategoriesChanged: ((filteredIdxs: number[], categories: string[]) => void) | null = null;

  get selectedCategories(): string[] | null {
    const idxs = this.legend?.selectedCategories;
    if (!idxs || !this.column || idxs.length === 0)
      return null;
    const categories = this.column.categories;
    return idxs.map((i) => categories[i]);
  }

  update(column: DG.Column): void {
    this.column = column;
    $(this.legendDiv).empty();
    this.legend = DG.Legend.create(column);
    const {categories} = column;
    this.legend.onViewerLegendChanged = () => {
      let filteredIdxs = this.legend!.selectedCategories;
      if (!filteredIdxs)
        return;
      if (filteredIdxs.length === 0)
        filteredIdxs = Array.from({length: categories.length}, (_, i) => i);
      this.onCategoriesChanged?.(filteredIdxs, categories);
    };
    this.legendDiv.appendChild(this.legend.root);
    $(this.legend.root).addClass('charts-legend');
  }

  show(chartDom?: HTMLElement): void {
    $(this.legendDiv).show();
    if (chartDom)
      $(chartDom).css('marginRight', '100px');
  }

  hide(chartDom?: HTMLElement): void {
    $(this.legendDiv).hide();
    if (chartDom)
      $(chartDom).css('marginRight', '');
  }

  switchVisibility(mode: VisibilityMode, column: DG.Column, chartDom?: HTMLElement): void {
    const autoShow = column.matches(DG.TYPE.CATEGORICAL) && column.categories.length < 100;
    if (mode === VISIBILITY_MODE.ALWAYS || (mode === VISIBILITY_MODE.AUTO && autoShow))
      this.show(chartDom);
    else
      this.hide(chartDom);
  }
}
