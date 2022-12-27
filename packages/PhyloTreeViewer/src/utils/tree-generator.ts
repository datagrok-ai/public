import {NodeType} from '@datagrok-libraries/bio';

export function generateTree(size: number): NodeType {
  function placeNode(currentNode: NodeType, newNode: NodeType): void {
    if (currentNode.children!.length < 2) {
      currentNode.children!.push(newNode);
    } else {
      const rnd: number = Math.random();
      const tgtNodeI = Math.floor(rnd / (1 / currentNode.children!.length));
      const tgtNode = currentNode.children![tgtNodeI];
      placeNode(tgtNode, newNode);
    }
  }

  const root = {
    name: `node-0`,
    branch_length: Math.random(),
    children: [],
  };

  for (let nodeI = 1; nodeI < size; nodeI++) {
    const newNode = {
      name: `node-${nodeI}`,
      branch_length: Math.random(),
      children: [],
    };
    placeNode(root, newNode);
  }

  return root;
}
