declare module 'phylotree' {
  export interface PhylotreeNode {
    id: string;
    name: string;
    parent: PhylotreeNode;
    children: PhylotreeNode[];
    attribute: string; // height / branch length
    annotation: string;
  }
}

declare module 'phylotree/src/formats/newick' {
  import {PhylotreeNode} from 'phylotree';

  function newickParser(newick: string): { json: PhylotreeNode, error: string | null } ;

  export = newickParser; // default
}

declare module 'phylotree' {
  class phylotree {

    get display(): TreeRender {}

    constructor(nwk: string, options: {});

    modify_selection(selection: string);

    path_to_root(node: PhylotreeNode);

    select_all_descendants(node: PhylotreeNode, flag1: boolean, flag2: boolean);

    getNodes(): PhylotreeNode[];

    render(props: { [propName: string]: any });
  }

  class TreeRender {
    show();

    modifySelection(predicate: () => boolean);
  }
}