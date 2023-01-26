declare module 'NGL' {
  export type LoaderParameters = {
    ext?: string,
    compressed?: boolean,
    binary?: boolean,
    name?: string,
    defaultRepresentation?: boolean
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
    addRepresentation(type: string, object: object, params?: RepresentationParameters): RepresentationComponent;

    autoView(duration?: number): undefined;

    removeAllRepresentations(): undefined;
  }

  export class RepresentationComponent extends Component {

  }

  export class Stage {
    get compList(): Component[];

    //TODO: Find out is host arg mandatory
    constructor(host: HTMLElement);

    handleResize(): undefined;

    removeComponent(component: Component): undefined;

    removeAllComponents(type?: string): undefined;

    async loadFile(path: String | File | Blob, params: LoaderParameters): Promise<void>;

    dispose(): undefined;
  }
}