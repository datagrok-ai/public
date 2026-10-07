import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: helper
//input: int a
//output: int b
export function helper(a: number): number {
  return a;
}

//name: Plain App
//tags: app
//meta.browsePath: Misc
export const plainApp = async (): Promise<void> => {
  await grok.functions.call('Demo:toMolfile', {mol: 'C'});
};
