import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NodeType} from './types';

export function isLeaf(node: NodeType) {
  return !node.children || node.children.length == 0;
}
