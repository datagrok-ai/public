import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {map, NORMALIZED, ADDITIONAL_MODS_COL_NAMES, BASE_MODIFICATIONS, USER_GROUP_NAME} from './constants';

export async function isCurrentUserAppAdmin() {
  const userGroup = await grok.dapi.groups.filter(USER_GROUP_NAME).first();
  if (userGroup == undefined)
    return false;
  const membersLogins = userGroup.members.map((e) => e.name.toLowerCase());
  const adminsLogins = userGroup.adminMembers.map((e) => e.name.toLowerCase());
  const allLogins = membersLogins.concat(adminsLogins);
  const loginOfCurrentUser = (await grok.dapi.users.current()).login.toLowerCase();
  return allLogins.includes(loginOfCurrentUser);
}

export function normalizeSequence(sequence: string, synthesizer: string | null, technology: string | null,
  additionalModsDf: DG.DataFrame): string {
  const additionalCodesCol = additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION);
  const baseModifsCol = additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.BASE_MODIFICATION);

  const codes = (technology == null) ?
    getAllCodesOfSynthesizer(synthesizer!).concat(additionalCodesCol.toList()) :
    Object.keys(map[synthesizer!][technology]);

  for (let i = 0; i < additionalModsDf.rowCount; i++) {
    NORMALIZED[additionalCodesCol.getString(i)] = (baseModifsCol.get(i) != BASE_MODIFICATIONS.NO) ?
      baseModifsCol.get(i) : '';
  }

  for (let i = 0; i < codes.length; i++) {
    codes[i] = codes[i].replace('(', '\\(');
    codes[i] = codes[i].replace(')', '\\)');
  }

  const sortedCodes = sortByStringLengthInDescOrder(codes);
  const regExp = new RegExp('(' + sortedCodes.join('|') + ')', 'g');
  return sequence.replace(regExp, function(code) {return NORMALIZED[code];});
}

export function getAllCodesOfSynthesizer(synthesizer: string): string[] {
  let codes: string[] = [];
  for (const technology of Object.keys(map[synthesizer]))
    codes = codes.concat(Object.keys(map[synthesizer][technology]));
  return codes;
}

export function deleteWord(sequence: string, searchTerm: string): string {
  let n = sequence.search(searchTerm);
  while (sequence.search(searchTerm) > -1) {
    n = sequence.search(searchTerm);
    sequence = sequence.substring(0, n) + sequence.substring(n + searchTerm.length, sequence.length);
  }
  return sequence;
}

export function saveAsCsv(table: DG.DataFrame): void {
  const link = document.createElement('a');
  link.setAttribute('href', 'data:text/csv;charset=utf-8,\uFEFF' + encodeURI(table.toCsv()));
  link.setAttribute('download', `Oligo Properties ${new Date()}.csv`);
  link.click();
}

export function sortByStringLengthInDescOrder(array: string[]): string[] {
  return array.sort(function(a, b) {return b.length - a.length;});
}

export function mergeOptions(obj1: {[ind: string]: number}, obj2: {[ind: string]: number}): {[ind: string]: number} {
  const obj3: {[index: string]: number} = {};
  for (const attrname in obj1) {
    if (Object.prototype.hasOwnProperty.call(obj1, attrname))
      obj3[attrname] = obj1[attrname];
  }
  for (const attrname in obj2) {
    if (Object.prototype.hasOwnProperty.call(obj2, attrname))
      obj3[attrname] = obj2[attrname];
  }
  return obj3;
}

export function stringify(items: string[]): string {
  return '["' + items.join('", "') + '"]';
}
