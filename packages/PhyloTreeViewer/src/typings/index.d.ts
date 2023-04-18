declare module '@phylocanvas/phylocanvas.gl' {
  export class PhylocanvasGL {
    get deck(): Deck;

    get view(): HTMLDivElement;

    get props(): { [propName: string]: any };

    constructor(element: HTMLElement, props: { [propName: string]: any });

    render(): void;

    fitInCanvas(): void;

    resume(): void;

    setProps(props: { [propName: string]: any }): void;

    selectNode: (nodeOrId: any, append: boolean = false) => void;

    destroy(): void;

    getBranchScale(...arguments): number;
  }
}
