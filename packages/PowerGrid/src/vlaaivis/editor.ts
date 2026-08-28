import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {createMpoRow, MpoRow} from '@datagrok-libraries/statistics/src/mpo/editors/mpo-property-row';

import {Subscription} from 'rxjs';

import {PieChartSettings, Sector} from '../sparklines/piechart';
import {VlaaiVisChange, VlaaiVisModel} from './model';
import {DEFAULTS, LABELS, PANE_HEADER_SELECTOR, TOOLTIPS} from './constants';

import '../../css/vlaaivis.css';

function clamp01(v: number | null): number {
  return Math.max(0, Math.min(1, v ?? 0));
}

function draggedProperty(o: any): string | null {
  return typeof o?.vlaaivisProperty === 'string' ? o.vlaaivisProperty : null;
}

/// One accordion pane per sector, each holding the MPO rows the profile editor builds from.
export class VlaaiVisEditor {
  readonly root = ui.div([], 'power-grid-vlaaivis');

  private model: VlaaiVisModel;
  private gc: DG.GridColumn;
  private body = ui.div([]);

  private rows: MpoRow[] = [];
  private panes = new Map<Sector, DG.AccordionPane>();
  private collapsed = new WeakSet<Sector>();
  private subs: Subscription[] = [];

  constructor(settings: PieChartSettings, gc: DG.GridColumn) {
    this.gc = gc;
    this.model = new VlaaiVisModel(settings, gc.grid.dataFrame);
    this.root.append(this.buildBoundsForm(), this.body);
    this.subs.push(this.model.onChanged.subscribe((change) => this.onModelChanged(change)));
    this.render();
  }

  refresh(): void {
    this.model.syncColumns();
  }

  detach(): void {
    this.disposeRows();
    for (const sub of this.subs)
      sub.unsubscribe();
    this.subs = [];
  }

  private onModelChanged(change: VlaaiVisChange): void {
    if (change === VlaaiVisChange.Structure)
      this.render();
    this.gc.grid.invalidate();
  }

  private buildBoundsForm(): HTMLElement {
    const bound = (caption: string, key: 'lowerBound' | 'upperBound', tooltip: string) => {
      const input = ui.input.float(caption, {value: this.model.bound(key), min: 0, max: 1, showSlider: false,
        onValueChanged: (v) => this.model.setBound(key, clamp01(v))});
      input.setTooltip(tooltip);
      return input;
    };
    return ui.divV([
      bound(LABELS.LOWER_BOUND, 'lowerBound', TOOLTIPS.LOWER_BOUND),
      bound(LABELS.UPPER_BOUND, 'upperBound', TOOLTIPS.UPPER_BOUND),
    ]);
  }

  private render(): void {
    this.rememberCollapsed();
    this.disposeRows();
    this.panes.clear();
    ui.empty(this.body);

    if (this.model.sectors.length === 0) {
      this.body.append(this.buildEmptyState());
      return;
    }

    const accordion = ui.accordion();
    for (const sector of this.model.sectors)
      this.addDropPane(accordion, sector, () => ui.divV(sector.subsectors.map((p) => this.buildRow(p.name))));

    const unassigned = this.model.unassigned;
    this.addDropPane(accordion, null, () => unassigned.length > 0 ?
      ui.divV(unassigned.map((name) => this.makeDraggable(ui.divText(name), name))) :
      ui.divText(LABELS.DROP_HINT, 'power-grid-vlaaivis-hint'));

    this.body.append(accordion.root,
      ui.divH([ui.icons.add(() => this.addSector(), TOOLTIPS.NEW_SECTOR)], 'power-grid-vlaaivis-add'));
  }

  private rememberCollapsed(): void {
    for (const [sector, pane] of this.panes)
      pane.expanded ? this.collapsed.delete(sector) : this.collapsed.add(sector);
  }

  private disposeRows(): void {
    for (const row of this.rows)
      row.sub.unsubscribe();
    this.rows = [];
  }

  private buildEmptyState(): HTMLElement {
    return ui.divV([
      ui.iconFA('chart-pie'),
      ui.h3('No sectors yet'),
      ui.p('A sector is a colored group of columns. The length of each wedge is that column\'s desirability ' +
        'score, and the sector\'s share of the circle is the sum of its weights.'),
      ui.divH([
        ui.bigButton(`Auto-group first ${DEFAULTS.AUTO_GROUP_COLUMNS}`, () => this.autoGroup()),
        ui.button('New sector', () => this.addSector()),
      ], 'power-grid-vlaaivis-actions'),
    ], 'statistics-mpo-empty-state');
  }

  private addDropPane(accordion: DG.Accordion, sector: Sector | null, content: () => HTMLElement): void {
    const pane = accordion.addPane(sector?.name ?? LABELS.UNASSIGNED, content,
      !sector || !this.collapsed.has(sector));
    this.makeDroppable(pane.root, sector);
    if (sector)
      this.panes.set(sector, pane);

    const header = pane.root.querySelector(PANE_HEADER_SELECTOR) as HTMLElement | null;
    if (!header)
      return;
    ui.tooltip.bind(header, sector ? TOOLTIPS.SECTOR : TOOLTIPS.UNASSIGNED);
    if (sector) {
      header.prepend(this.buildSwatch(sector));
      header.append(this.buildSectorButtons(sector));
    }
  }

  private buildSectorButtons(sector: Sector): HTMLElement {
    const rename = ui.icons.edit(() => this.renameSector(sector), TOOLTIPS.RENAME);
    const remove = ui.icons.delete(() => setTimeout(() => this.model.deleteSector(sector)), TOOLTIPS.DELETE);
    const buttons = ui.divH([rename, remove], 'power-grid-vlaaivis-sector-buttons');
    buttons.addEventListener('click', (e) => e.stopPropagation());
    return buttons;
  }

  private buildSwatch(sector: Sector): HTMLElement {
    const swatch = ui.div([], 'power-grid-vlaaivis-swatch');
    swatch.style.backgroundColor = sector.sectorColor;
    ui.tooltip.bind(swatch, TOOLTIPS.COLOR);

    swatch.addEventListener('click', (e) => {
      e.stopPropagation();
      const original = sector.sectorColor;
      // Paint straight onto the render model so the wedges follow the picker; the model call on OK persists it.
      const preview = (color: string) => {
        sector.sectorColor = color;
        swatch.style.backgroundColor = color;
        this.gc.grid.invalidate();
      };
      ui.showColorPicker(DG.Color.fromHtml(original), (c) => preview(DG.Color.toHtml(c)),
        () => this.model.setSectorColor(sector, sector.sectorColor), () => preview(original));
    });

    return swatch;
  }

  private buildRow(name: string): HTMLElement {
    const row = createMpoRow(this.model.locate(name)!.property, {
      name: () => name,
      column: () => this.model.column(name),
      height: DEFAULTS.PLOT_HEIGHT,
      design: true,
      onChanged: () => this.model.notifyChanged(name),
      onReplaced: (replacement) => this.model.replaceProperty(name, replacement),
    });

    const label = ui.divText(name, 'power-grid-vlaaivis-name');
    ui.tooltip.bind(label, () => `${name} — ${TOOLTIPS.DRAG}`);
    this.makeDraggable(label, name);
    row.propertyCell.append(label);
    row.weightInput.setTooltip(TOOLTIPS.WEIGHT);

    this.rows.push(row);
    return row.root;
  }

  private makeDraggable(element: HTMLElement, name: string): HTMLElement {
    ui.makeDraggable(element, {
      getDragObject: () => ({vlaaivisProperty: name}),
      getDragCaption: () => name,
    });
    return element;
  }

  private makeDroppable(element: HTMLElement, sector: Sector | null): void {
    ui.makeDroppable(element, {
      acceptDrop: (o) => draggedProperty(o) !== null,
      doDrop: (args) => {
        const name = draggedProperty(args.dragObject);
        // The panel is rebuilt on assignment, so it has to wait for the drag session to finish.
        if (name)
          setTimeout(() => this.model.assign(name, sector));
      },
      dropSuggestion: sector ? `Add to ${sector.name}` : 'Remove from sector',
    });
  }

  private autoGroup(): void {
    this.model.autoGroup(DEFAULTS.AUTO_GROUP_COLUMNS);
  }

  private addSector(): void {
    this.nameDialog('New sector', '', (name) => this.model.addSector(name));
  }

  private renameSector(sector: Sector): void {
    this.nameDialog('Rename sector', sector.name, (name) => this.model.renameSector(sector, name));
  }

  private nameDialog(title: string, current: string, apply: (name: string) => void): void {
    const input = ui.input.string('Name', {value: current});
    const error = (value: string): string | null => {
      const name = value.trim();
      if (!name)
        return 'Enter a name';
      return name !== current && this.model.sector(name) ? `Sector "${name}" already exists` : null;
    };
    input.addValidator(error);
    ui.dialog(title).add(input).onOK(() => {
      if (!error(input.value))
        apply(input.value.trim());
    }).show();
  }
}
