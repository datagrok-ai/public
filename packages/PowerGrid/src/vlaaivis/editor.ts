import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {createMpoRow, MpoRow} from '@datagrok-libraries/statistics/src/mpo/editors/mpo-property-row';

import {Subscription} from 'rxjs';

import {PieChartSettings, Sector} from '../sparklines/piechart';
import {VlaaiVisChange, VlaaiVisModel} from './model';
import {DEFAULTS, LABELS, PANE_HEADER_SELECTOR, TOOLTIPS} from './constants';

import '../../css/vlaaivis.css';

function draggedProperty(o: any): string | null {
  return typeof o?.vlaaivisProperty === 'string' ? o.vlaaivisProperty : null;
}

/// One accordion pane per sector, each holding the MPO rows the profile editor builds from.
export class VlaaiVisEditor {
  readonly root = ui.div([], 'power-grid-vlaaivis');
  readonly boundsInputs: DG.InputBase<number | null>[];

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
    if (this.model.sectors.length === 0)
      this.model.autoGroup(DEFAULTS.AUTO_GROUP_COLUMNS);
    this.boundsInputs = this.buildBoundsInputs();
    this.root.append(
      ui.divH([ui.iconFA('info-circle'), ui.divText(LABELS.TIP)], 'power-grid-vlaaivis-tip'),
      this.body);
    this.subs.push(this.model.onChanged.subscribe((change) => this.onModelChanged(change)));
    this.render();
  }

  refresh(): void {
    this.model.syncColumns();
  }

  detach(): void {
    for (const input of this.boundsInputs)
      input.root.remove();
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

  private buildBoundsInputs(): DG.InputBase<number | null>[] {
    const bound = (caption: string, key: 'lowerBound' | 'upperBound', tooltip: string) =>
      ui.input.float(caption, {value: this.model.bound(key), min: 0, max: 1, showSlider: false}).setTooltip(tooltip);
    const lower = bound(LABELS.LOWER_BOUND, 'lowerBound', TOOLTIPS.LOWER_BOUND);
    const upper = bound(LABELS.UPPER_BOUND, 'upperBound', TOOLTIPS.UPPER_BOUND);

    const values = (): [number, number] => [lower.value ?? 0, upper.value ?? 0];
    const inverted = () => values()[0] > values()[1];
    lower.addValidator(() => inverted() ? LABELS.LOWER_ABOVE_UPPER : null);
    upper.addValidator(() => inverted() ? LABELS.UPPER_BELOW_LOWER : null);

    const apply = (other: DG.InputBase<number | null>) => {
      other.validate();
      if (!inverted())
        this.model.setBounds(...values());
    };
    this.subs.push(lower.onChanged.subscribe(() => apply(upper)), upper.onChanged.subscribe(() => apply(lower)));
    return [lower, upper];
  }

  private render(): void {
    this.rememberCollapsed();
    this.disposeRows();
    this.panes.clear();
    ui.empty(this.body);

    const accordion = ui.accordion();
    for (const sector of this.model.sectors)
      this.addDropPane(accordion, sector, () => ui.divV(sector.subsectors.map((p) => this.buildRow(p.name))));
    accordion.root.append(
      ui.link(LABELS.ADD_SECTOR, () => this.addSector(), TOOLTIPS.NEW_SECTOR, 'power-grid-vlaaivis-add'));

    const unassigned = this.model.unassigned;
    this.addDropPane(accordion, null, () => unassigned.length > 0 ?
      ui.divV(unassigned.map((name) => this.makeDraggable(ui.divText(name), name))) :
      ui.divText(LABELS.NO_UNASSIGNED, 'power-grid-vlaaivis-empty'));

    this.body.append(accordion.root);
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
