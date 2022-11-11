declare module 'logojs-react' {
  function embedProteinLogo(div: HTMLElement, props: { [propName: string]: any });
}

// declare module '@phylocanvas/phylocanvas.gl' {
//   import {Deck} from '@deck.gl/core/typed';
//
//   class PhylocanvasGL {
//     get deck(): Deck;
//
//     constructor(element: HTMLElement, props: { [propName: string]: any });
//
//     setProps(props: { [propName: string]: any });
//
//     selectNode: (nodeOrId: any, append: boolean = false) => void;
//
//     destroy(): void;
//   }
//
//   const TreeTypes;
//   const Shapes;
// }