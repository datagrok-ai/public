import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type drugBankSearchTypes = 'similarity' | 'substructure';
export type drugBankSearchResult = DG.DataFrame | null;

export function getTooltip(value: string) {

  return ui.divText(`Common name: ${value}\nClick to open in DrugBank Online`);
}
