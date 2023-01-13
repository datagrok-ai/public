declare module 'NGL' {
  export type LoaderParameters = {
    ext: string,
    compressed: boolean,
    binary: boolean,
    name: string
  };

  export class ColormakerRegistryClass {
    add(id: String, scheme: Colormaker): undefined

    addScheme(scheme: Function | Colormaker, label: String): String

    addSelectionScheme(dataList: Array, label?: String): String
  }

  export const ColormakerRegistry: ColormakerRegistryClass;

  export type RepresentationParameters = {
    [p: name]: any
  }

  export class Component {
    addRepresentation(type: string, object: object, params?: RepresentationParameters): RepresentationComponent

    autoView();
  }

  export class RepresentationComponent extends Component {

  }

  export class Stage {
    get compList(): Component[];

    constructor(host: HTMLElement);

    removeAllComponents(type?: string): undefined;

    async loadFile(path: String | File | Blob, params: LoaderParameters): Promise<void>;
  }
}