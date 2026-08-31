import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {T_COMPOUND, T_PROGRAM_COMPOUND, T_TAG} from '../domain/constants';

/** Registration codes of the compounds linked to [programId]. */
export async function programCompoundCodes(programId: string): Promise<string[]> {
  const links = await grok.dapi.domains.table(T_PROGRAM_COMPOUND).query({
    filter: DG.cond('program_id', '=', programId), limit: 100});
  if (links.length === 0)
    return [];
  return (await grok.dapi.domains.table(T_COMPOUND).query({
    filter: DG.or(...links.map((l: any) => DG.cond('id', '=', l.compound_id))), limit: 100}))
    .map((c: any) => c.registration_code).sort();
}

/** Chips over the compound registry; empty text suggests [linkedCodes] (the
 * program's own compounds). allowNew is off — a typo must not mint a compound. */
export function compoundsChipInput(linkedCodes: () => string[], value: string[] = []): DG.InputBase {
  const input = ui.input.tags('Molecules', {
    value: value as any,
    allowNew: false,
    getSuggestions: async (text: string) => {
      const t = text?.trim() ?? '';
      if (t === '')
        return linkedCodes();
      const rows = await grok.dapi.domains.table(T_COMPOUND).query({
        filter: DG.cond('registration_code', 'like', `%${t}%`), sort: 'registration_code', limit: 20});
      return rows.map((r: any) => r.registration_code);
    },
  });
  input.setTooltip('Registration codes from the compound registry');
  return input;
}

/** Chips over the tag registry; a new name creates the tag on save. */
export function tagsChipInput(value: string[] = []): DG.InputBase {
  const input = ui.input.tags('Tags', {
    value: value as any,
    allowNew: true,
    createNewItem: (text: string) => text.trim(),
    getSuggestions: async (text: string) => {
      const rows = await grok.dapi.domains.table(T_TAG).query({
        filter: DG.cond('name', 'like', `%${text?.trim() ?? ''}%`), sort: 'name', limit: 20});
      return rows.map((r: any) => r.name);
    },
  });
  input.setTooltip('Existing tags are suggested; a new name creates the tag');
  return input;
}

/** A chip input's value as trimmed non-empty strings. */
export const chips = (value: unknown): string[] =>
  ((value ?? []) as string[]).map((s) => `${s}`.trim()).filter((s) => s.length > 0);
