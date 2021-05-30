/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export async function entitySearch(s: string): Promise<any[]> {
  return [`entity results for ${s}`];
}

export async function googleSearch(s: string): Promise<any[]> {
  return [`google results for  ${s}`];
}

export async function functionSearch(s: string): Promise<any[]> {
  s = s.toLowerCase();
  return DG.Func.find()
    .filter(value =>
      value.name.toLowerCase().includes(s) ||
      value.description.toLowerCase().includes(s))
    .map((f) => f.name);
}

