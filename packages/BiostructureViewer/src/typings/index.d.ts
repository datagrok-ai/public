declare module 'NGL' {
  import {Signal} from 'signals';
  import {Box3, Scene, WebGLRenderer} from 'three';

  export type LoaderParameters = {
    ext: string,
    compressed: boolean,
    binary: boolean,
    name: string,
    defaultRepresentation: boolean
  };

  export type Colormaker = {};

  export class ColormakerRegistryClass {
    add(id: String, scheme: Colormaker): undefined

    addScheme(scheme: Function | Colormaker, label: String): String

    addSelectionScheme(dataList: [string, string][], label?: String): String
  }

  export const ColormakerRegistry: ColormakerRegistryClass;

  export type RepresentationParameters = { [k: string]: any }

  export type StructureRepresentationType = (
    'angle' | 'axes' | 'backbone' | 'ball+stick' | 'base' | 'cartoon' | 'contact' | 'dihedral' |
    'distance' | 'helixorient' | 'hyperball' | 'label' | 'licorice' | 'line' | 'surface' |
    'ribbon' | 'rocket' | 'rope' | 'spacefill' | 'trace' | 'tube' | 'unitcell'
    ) // from NGL lib

  export class Component {
    addRepresentation(type: StructureRepresentationType, params?: RepresentationParameters): RepresentationElement;

    autoView(duration?: number): undefined;

    removeAllRepresentations(): undefined;
  }

  export class RepresentationElement extends Component {

  }

  export class Stage {
    viewer: Viewer;

    get compList(): Component[];

    //TODO: Find out is host arg mandatory
    constructor(host: HTMLElement);

    handleResize(): undefined;

    removeComponent(component: Component): undefined;

    removeAllComponents(type?: string): undefined;

    loadFile(path: String | File | Blob, params: Partial<LoaderParameters>): Promise<void>;

    dispose(): undefined;
  }

  export class Viewer {
    signals: ViewerSignals;

    container: HTMLElement;
    wrapper: HTMLElement;

    sampleLevel: number;

    stats: Stats;

    width: number;
    height: number;

    scene: Scene;

    renderer: WebGLRenderer;
    boundingBox: Box3;

    rendering: boolean;
    renderPending: boolean;

    render(picking: boolean): void;
    requestRender(): void;

    setSize(width: number, height: number): void;
  }

  export class Stats {
    signals: {
      updated: Signal
    };

    maxDuration: number;
    minDuration: number;
    avgDuration: number;
    lastDuration: number;

    prevFpsTime: number;
    lastFps: number;
    lastFrames: number;
    frames: number;
    count: number;

    startTime: number;
    currentTime: number;
  }

  export interface ViewerSignals {
    ticked: Signal,
    rendered: Signal
  }
}
