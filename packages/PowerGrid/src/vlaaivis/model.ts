import * as DG from 'datagrok-api/dg';

import {Subject} from 'rxjs';

import {
  createDefaultCategorical, createDefaultNumerical, migrateDesirability, PropertyDesirability
} from '@datagrok-libraries/statistics/src/mpo/mpo';

import {PieChartSettings, Sector, Subsector} from '../sparklines/piechart';
import {DEFAULTS, VLAAIVIS_METADATA_TAG} from './constants';

type ColumnMetadata = PropertyDesirability & {groupName?: string; sectorColor?: string};

export enum VlaaiVisChange { Structure, Value }

/// Owns the VlaaiVis configuration. `PieChartSettings.sectors` is the model — it is what the renderer
/// reads — and the column tags mirror it so a configuration survives a reopened project.
export class VlaaiVisModel {
  readonly onChanged = new Subject<VlaaiVisChange>();

  constructor(private settings: PieChartSettings, private dataFrame: DG.DataFrame) {
    this.settings.sectors ??= {
      lowerBound: DEFAULTS.LOWER_BOUND, upperBound: DEFAULTS.UPPER_BOUND,
      sectors: this.readColumnTags(), values: null,
    };
  }

  get sectors(): Sector[] { return this.settings.sectors!.sectors; }

  bound(key: 'lowerBound' | 'upperBound'): number { return this.settings.sectors![key]; }

  setBound(key: 'lowerBound' | 'upperBound', value: number): void {
    this.settings.sectors![key] = value;
    this.onChanged.next(VlaaiVisChange.Value);
  }

  column(name: string): DG.Column | null { return this.dataFrame.col(name); }

  sector(name: string): Sector | undefined { return this.sectors.find((s) => s.name === name); }

  locate(name: string): {sector: Sector; property: Subsector} | null {
    for (const sector of this.sectors) {
      const property = sector.subsectors.find((s) => s.name === name);
      if (property)
        return {sector, property};
    }
    return null;
  }

  get unassigned(): string[] {
    const assigned = new Set(this.sectors.flatMap((s) => s.subsectors.map((p) => p.name)));
    return this.settings.columnNames.filter((name) => !assigned.has(name));
  }

  addSector(name: string): Sector {
    const sector = this.newSector(name);
    this.changed(VlaaiVisChange.Structure, []);
    return sector;
  }

  renameSector(sector: Sector, name: string): void {
    sector.name = name;
    this.changed(VlaaiVisChange.Structure, []);
  }

  deleteSector(sector: Sector): void {
    this.settings.sectors!.sectors = this.sectors.filter((s) => s !== sector);
    this.changed(VlaaiVisChange.Structure, []);
  }

  setSectorColor(sector: Sector, color: string): void {
    sector.sectorColor = color;
    this.changed(VlaaiVisChange.Value, sector.subsectors.map((p) => p.name));
  }

  /// Moves a column into a sector, or out of every sector when `target` is null.
  assign(name: string, target: Sector | null): void {
    const found = this.locate(name);
    if ((found?.sector ?? null) === target)
      return;
    const property = found?.property ?? this.createProperty(name);
    if (found)
      found.sector.subsectors = found.sector.subsectors.filter((p) => p !== property);
    if (target)
      target.subsectors.push(property);
    this.changed(VlaaiVisChange.Structure, []);
  }

  autoGroup(count: number): void {
    for (const name of this.unassigned.slice(0, count))
      this.newSector(name).subsectors.push(this.createProperty(name));
    this.changed(VlaaiVisChange.Structure, []);
  }

  /// The desirability editors mutate the property in place; this republishes what they changed.
  notifyChanged(name: string): void {
    this.changed(VlaaiVisChange.Value, [name]);
  }

  replaceProperty(name: string, prop: PropertyDesirability): void {
    const found = this.locate(name);
    if (!found)
      return;
    found.sector.subsectors[found.sector.subsectors.indexOf(found.property)] = {...prop, name};
    this.changed(VlaaiVisChange.Structure, []);
  }

  private newSector(name: string): Sector {
    const sector: Sector = {
      name,
      sectorColor: DG.Color.toHtml(DG.Color.getCategoricalColor(this.sectors.length)),
      subsectors: [],
    };
    this.sectors.push(sector);
    return sector;
  }

  private createProperty(name: string): Subsector {
    const col = this.column(name);
    const desirability = col?.isCategorical ?
      createDefaultCategorical(DEFAULTS.WEIGHT, col.categories.map((c: string) => ({name: c, desirability: 1}))) :
      createDefaultNumerical(DEFAULTS.WEIGHT, col?.min ?? 0, col?.max ?? 1);
    return {...desirability, name};
  }

  private changed(type: VlaaiVisChange, names: string[]): void {
    this.writeColumnTags(type === VlaaiVisChange.Structure ? this.settings.columnNames : names);
    this.onChanged.next(type);
  }

  private writeColumnTags(names: string[]): void {
    for (const name of names) {
      const col = this.column(name);
      if (!col)
        continue;
      const found = this.locate(name);
      const meta: ColumnMetadata | null = found ?
        {...found.property, groupName: found.sector.name, sectorColor: found.sector.sectorColor} : null;
      col.setTag(VLAAIVIS_METADATA_TAG, meta ? JSON.stringify(meta) : '');
    }
  }

  private readColumnTags(): Sector[] {
    const sectors: Sector[] = [];
    for (const name of this.settings.columnNames) {
      const tag = this.column(name)?.getTag(VLAAIVIS_METADATA_TAG);
      if (!tag)
        continue;

      const {groupName, sectorColor, ...property} = migrateDesirability(JSON.parse(tag)) as ColumnMetadata;
      if (!groupName)
        continue;

      let sector = sectors.find((s) => s.name === groupName);
      if (!sector) {
        sector = {name: groupName, subsectors: [],
          sectorColor: sectorColor ?? DG.Color.toHtml(DG.Color.getCategoricalColor(sectors.length))};
        sectors.push(sector);
      }
      sector.subsectors.push({...property, name});
    }
    return sectors;
  }
}
