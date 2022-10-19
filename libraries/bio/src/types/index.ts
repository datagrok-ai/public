// import {NodeType} from '@phylocanvas/phylocanvas.gl';

interface NodeType {
  name: string;
  children: NodeType[];
  branch_length: number;
  isLeaf: boolean;
}

export {NodeType as NodeType};