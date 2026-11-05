import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {NameMapping, UnitMapping} from './types';
import {isValidUnitMapping} from './matching';
import {UnmatchedItem} from './selection';

const UNMATCHED_SHOWN = 10;

/** Edits name and unit mappings; OK hands back the complete, valid rows only. */
export function showMappingsDialog(
  nameMappings: NameMapping[],
  unitMappings: UnitMapping[],
  unmatched: UnmatchedItem[],
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
  const header = (title: string, hint: string) =>
    ui.divH([ui.divText(title, 'c2-mapping-title'), ui.divText(hint, 'c2-mapping-hint')], 'c2-mapping-header');
  const section = (title: string, hint: string, list: HTMLElement, onAdd: () => void) => ui.divV([
    header(title, hint),
    list,
    ui.button('Add', onAdd),
  ], 'c2-mapping-section');
  // fills the first empty row, so the starter row is used before new ones are added
  const mapName = (name: string) => {
    const empty = nameRows.find((row) => !row.from.value && !row.to.value);
    if (empty)
      empty.from.value = name;
    else
      addNameRow(name);
  };
  const mapUnits = (units: string) => {
    const empty = unitRows.find((row) => !row.from.value && row.factor.value == null && !row.to.value);
    if (empty)
      empty.from.value = units;
    else
      addUnitRow(units);
  };
  const unmatchedRow = (item: UnmatchedItem) => ui.divH([
    ui.divText(item.units ? `${item.name} (${item.units})` : item.name, 'c2-mapping-unmatched-name'),
    ui.divText(item.runCount === 1 ? '1 run' : `${item.runCount} runs`, 'c2-mapping-hint'),
    ui.link('Map name', () => mapName(item.name)),
    ...item.units ? [ui.link('Map units', () => mapUnits(item.units!))] : [],
  ], 'c2-mapping-row');
  const unmatchedList = ui.divV(unmatched.slice(0, UNMATCHED_SHOWN).map(unmatchedRow), 'c2-mapping-list');
  if (unmatched.length > UNMATCHED_SHOWN) {
    const more = ui.link(`...more (${unmatched.length - UNMATCHED_SHOWN})`, () => {
      more.remove();
      unmatchedList.append(...unmatched.slice(UNMATCHED_SHOWN).map(unmatchedRow));
      unmatchedList.classList.add('c2-mapping-unmatched-scroll');
    });
    unmatchedList.append(more);
  }
  const dialog = ui.dialog('Mappings');
  if (unmatched.length > 0) {
    dialog.add(ui.divV([header('Not in any compare row', 'Map a name or units to let the value join a row.'),
      unmatchedList], 'c2-mapping-section'));
  }
  dialog
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
