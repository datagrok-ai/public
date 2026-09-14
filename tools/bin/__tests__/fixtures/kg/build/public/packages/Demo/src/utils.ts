import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export async function callThings(): Promise<void> {
  await grok.functions.call('Demo:toHelm', {});
  await grok.functions.call('demo:ToHELM', {});
  const f = DG.Func.byName('Plain:helper');
  const g = DG.Func.find({package: 'Nowhere', name: 'missing'});
  await grok.functions.eval('Demo:CalculatelogD');
  await grok.functions.call('Demo:ParetoFront', {});
  await grok.functions.call('Demo:TradeOff', {});
  console.log(f, g);
}
