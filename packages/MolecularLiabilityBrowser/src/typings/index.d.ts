declare module 'logojs-react' {
  function embedProteinLogo(div: HTMLElement, props: { [propName: string]: any });
}

declare module '@phylocanvas/phylocanvas.gl' {
  import {Deck} from '@deck.gl/core/typed';

  export class PhylocanvasGL {
    get deck(): Deck;

    get view(): HTMLDivElement;

    get props(): { [propName: string]: any };

    constructor(element: HTMLElement, props: { [propName: string]: any });

    render(): void;

    resume(): void;

    setProps(props: { [propName: string]: any }): void;

    selectNode: (nodeOrId: any, append: boolean = false) => void;

    destroy(): void;

    getBranchScale(...arguments): number;
  }

  module Utils {
    function treeTraversal(root: PhylocanvasTreeNode);
  }
}