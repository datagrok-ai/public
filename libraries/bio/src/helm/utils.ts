import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export function cleanupHelmSymbol(symbol: string): string {
  if (symbol.startsWith('[') && symbol.endsWith(']'))
    return symbol.slice(1, -1);
  return symbol;
}
