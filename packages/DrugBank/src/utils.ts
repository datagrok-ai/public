import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type drugBankSearchTypes = 'similarity' | 'substructure';
export type drugBankSearchResult = DG.DataFrame | null;

export function getTooltip(value: string) {
  const props = {
    'DRUGBANK_ID': value,
  };
  return ui.divV([ui.tableFromMap(props), ui.divText('Click to open in the store.')]);
}
