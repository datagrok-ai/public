// import {NodeType} from '@phylocanvas/phylocanvas.gl';

interface NodeType {
  name: string;
  children: NodeType[];
  branch_length: number;
}

export {NodeType as NodeType};