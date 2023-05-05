declare module 'NGL' {
  export type LoaderParameters = {
    ext: string,
    compressed: boolean,
    binary: boolean,
    name: string,
    defaultRepresentation: boolean
  };

  export class ColormakerRegistryClass {
    add(id: String, scheme: Colormaker): undefined

    addScheme(scheme: Function | Colormaker, label: String): String

    addSelectionScheme(dataList: Array, label?: String): String
  }

  export const ColormakerRegistry: ColormakerRegistryClass;

  export type RepresentationParameters = { [k: string]: any }

  export type StructureRepresentationType = (
    'angle' | 'axes' | 'backbone' | 'ball+stick' | 'base' | 'cartoon' | 'contact' | 'dihedral' |
    'distance' | 'helixorient' | 'hyperball' | 'label' | 'licorice' | 'line' | 'surface' |
    'ribbon' | 'rocket' | 'rope' | 'spacefill' | 'trace' | 'tube' | 'unitcell'
    ) // from NGL lib

  export class Component {
    addRepresentation(type: StructureRepresentationType, params?: RepresentationParameters = {}): RepresentationElement;

    autoView(duration?: number): undefined;

    removeAllRepresentations(): undefined;
  }

  export class RepresentationElement extends Component {

  }

  export class Stage {
    get compList(): Component[];

    //TODO: Find out is host arg mandatory
    constructor(host: HTMLElement);

    handleResize(): undefined;

    removeComponent(component: Component): undefined;

    removeAllComponents(type?: string): undefined;

    async loadFile(path: String | File | Blob, params: Partial<LoaderParameters>): Promise<void>;

    dispose(): undefined;
  }
}
