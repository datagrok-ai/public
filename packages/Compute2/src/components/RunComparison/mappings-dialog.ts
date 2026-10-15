import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {NameMapping, UnitMapping} from './types';
import {isValidUnitMapping} from './matching';

/** Edits name and unit mappings; OK hands back the complete, valid rows only. */
export function showMappingsDialog(
  nameMappings: NameMapping[],
  unitMappings: UnitMapping[],
  onApply: (names: NameMapping[], units: UnitMapping[]) => void,
): void {
  const nameRows: {from: DG.InputBase<string>, to: DG.InputBase<string>}[] = [];
  const unitRows: {
    from: DG.InputBase<string>, factor: DG.InputBase<number | null>, to: DG.InputBase<string>,
  }[] = [];
  const nameList = ui.divV([], 'c2-mapping-list');
  const unitList = ui.divV([], 'c2-mapping-list');
  const removable = <T, >(rows: T[], row: T, rowDiv: HTMLElement) => ui.icons.delete(() => {
    rows.splice(rows.indexOf(row), 1);
    rowDiv.remove();
  }, 'Remove mapping');
  const addNameRow = (from = '', to = '') => {
    const row = {from: ui.input.string('', {value: from}), to: ui.input.string('', {value: to})};
    (row.from.input as HTMLInputElement).placeholder = 'name';
    (row.to.input as HTMLInputElement).placeholder = 'same as';
    nameRows.push(row);
    const rowDiv = ui.divH([row.from.root, row.to.root], 'c2-mapping-row');
    rowDiv.append(removable(nameRows, row, rowDiv));
    nameList.append(rowDiv);
  };
  const addUnitRow = (from = '', factor?: number, to = '') => {
    const row = {
      from: ui.input.string('', {value: from}),
      factor: ui.input.float('', {value: factor}),
      to: ui.input.string('', {value: to}),
    };
    (row.from.input as HTMLInputElement).placeholder = 'units';
    (row.factor.input as HTMLInputElement).placeholder = 'factor';
    (row.to.input as HTMLInputElement).placeholder = 'units';
    row.from.root.classList.add('c2-mapping-unit');
    row.factor.root.classList.add('c2-mapping-factor');
    row.to.root.classList.add('c2-mapping-unit');
    unitRows.push(row);
    const rowDiv = ui.divH([ui.divText('1'), row.from.root, ui.divText('='), row.factor.root, row.to.root],
      'c2-mapping-row');
    rowDiv.append(removable(unitRows, row, rowDiv));
    unitList.append(rowDiv);
  };
  for (const {from, to} of nameMappings)
    addNameRow(from, to);
  for (const {from, factor, to} of unitMappings)
    addUnitRow(from, factor, to);
  if (nameRows.length === 0)
    addNameRow();
  if (unitRows.length === 0)
    addUnitRow();
  const section = (title: string, hint: string, list: HTMLElement, onAdd: () => void) => ui.divV([
    ui.divH([ui.divText(title, 'c2-mapping-title'), ui.divText(hint, 'c2-mapping-hint')], 'c2-mapping-header'),
    list,
    ui.button('Add', onAdd),
  ], 'c2-mapping-section');
  ui.dialog('Mappings')
    .add(section('Names', 'Mapped names match as the same value, ignoring case and separators.',
      nameList, () => addNameRow()))
    .add(section('Units', 'Values in one unit are multiplied into the other; works both ways.',
      unitList, () => addUnitRow()))
    .onOK(() => onApply(
      nameRows
        .map((row) => ({from: row.from.value.trim(), to: row.to.value.trim()}))
        .filter((row) => row.from && row.to),
      unitRows
        .map((row) => ({from: row.from.value.trim(), factor: row.factor.value ?? NaN, to: row.to.value.trim()}))
        .filter(isValidUnitMapping),
    ))
    .show({center: true, width: 520});
}
