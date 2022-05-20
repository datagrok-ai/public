import * as DG from "datagrok-api/dg";

export function names(columns: Iterable<DG.Column>): string[] {
    return Array.from(columns).map((c: any) => c.name)
}

export interface SummarySettingsBase {
    columnNames: string[];
}

type getSettingsFunc<Type extends SummarySettingsBase> = (gs: DG.GridColumn) => Type;

export function getSettingsBase<Type extends SummarySettingsBase>(gc: DG.GridColumn): Type {
    return gc.settings ??= {
        columnNames: names(gc.grid.dataFrame.columns.numerical),
    }
}

function getDataColumns<Type extends SummarySettingsBase>(
    gc: DG.GridColumn,
    getSettings: getSettingsFunc<Type>): DG.Column[] {
  return gc.grid.dataFrame.columns.byNames(getSettings(gc).columnNames);
}