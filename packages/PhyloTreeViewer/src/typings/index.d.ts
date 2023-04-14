declare module 'phylotree' {
  export interface PhylotreeNode {
    id: string;
    name: string;
    parent: PhylotreeNode;
    children: PhylotreeNode[];
    attribute: string; // height / branch length
    annotation: string;
    depth?: number;
  }
}

declare module 'phylotree/src/formats/newick' {
  import {PhylotreeNode} from 'phylotree';
  import {NO_NAME_ROOT} from '@datagrok-libraries/bio/src/trees/phylocanvas';

  /** Does not correct empty root name to {@link NO_NAME_ROOT}*/
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

    /** Get the root node */
    getRootNode(): PhylotreeNode

    /** Implementation returns root node, the same {@link getRootNode} */
    getNodes(): PhylotreeNode;

    /** Update with new hiearchy layout */
    update(json): void;

    render(props: { [propName: string]: any }): TreeRender;
  }

  class TreeRender {
    show(): SVGSVGElement;

    modifySelection(predicate: () => boolean);
  }
}

declare module '@phylocanvas/phylocanvas.gl' {
  export class PhylocanvasGL {
    get deck(): Deck;

    get view(): HTMLDivElement;

    get props(): { [propName: string]: any };

    constructor(element: HTMLElement, props: { [propName: string]: any });

    render(): void;

    fitInCanvas() : void;

    resume(): void;

    setProps(props: { [propName: string]: any }): void;

    selectNode: (nodeOrId: any, append: boolean = false) => void;

    destroy(): void;

    getBranchScale(...arguments): number;
  }
}

declare module '@phylocanvas/phylocanvas.gl' {
  export class PhylocanvasGL {
    get deck(): Deck;

    get view(): HTMLDivElement;

    get props(): { [propName: string]: any };

    constructor(element: HTMLElement, props: { [propName: string]: any });

    render(): void;

    resume(): void;

    fitInCanvas() : void;

    setProps(props: { [propName: string]: any }): void;

    selectNode: (nodeOrId: any, append: boolean = false) => void;

    destroy(): void;

    getBranchScale(...arguments): number;
  }
}

declare module '@phylocanvas/phylocanvas.gl' {
  export class PhylocanvasGL {
    get deck(): Deck;

    get view(): HTMLDivElement;

    get props(): { [propName: string]: any };

    constructor(element: HTMLElement, props: { [propName: string]: any });

    fitInCanvas() : void;

    render(): void;

    resume(): void;

    setProps(props: { [propName: string]: any }): void;

    selectNode: (nodeOrId: any, append: boolean = false) => void;

    destroy(): void;

    getBranchScale(...arguments): number;
  }
}
