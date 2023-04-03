import {LoadMolstarInterface} from "./MolstarPlugin";
import {PluginContext} from "molstar/lib/mol-plugin/context";
import {StructureRepresentationRegistry} from "molstar/lib/mol-repr/structure/registry";
import {ColorTheme} from "molstar/lib/mol-theme/color";

export type SaguaroPluginModelMapType = Map<string,{entryId: string; chains:Array<{label:string, auth:string; entityId: string; title: string; type:"polymer"|"water"|"branched"|"non-polymer"|"macrolide";}>;}>;

export interface SaguaroPluginInterface extends SaguaroPluginPublicInterface{
    init: (elementId: string, props?: any) => void;
    load: (args: LoadMolstarInterface) => void;
    pluginCall: (f:(plugin: PluginContext)=>void) => void;
    clear: () => void;
    setSelectCallback: (g:(flag?:boolean)=>void)=>void;
    setModelChangeCallback: (f:(modelMap:SaguaroPluginModelMapType)=>void)=>void;
    setHoverCallback:(g:(flag?:boolean)=>void)=>void;
    setRepresentationChangeCallback:(g:(flag?:boolean)=>void)=>void;
    unsetCallbacks:()=>void;
}

export interface SaguaroPluginPublicInterface {
    select(modelId:string, labelAsymId: string, x: number, y: number, mode: 'select'|'hover', operation:'set'|'add'): void;
    select(selection: Array<{modelId:string; labelAsymId: string; position: number;}>, mode: 'select'|'hover', operation:'add'|'set'): void;
    select(selection: Array<{modelId:string; labelAsymId: string; begin: number; end:number;}>, mode: 'select'|'hover', operation:'add'|'set'): void;
    clearSelection: (mode:'select'|'hover', option?:{modelId:string; labelAsymId:string;}) => void;
    createComponent(componentId: string, modelId:string, labelAsymId: string, begin: number, end : number, representationType: StructureRepresentationRegistry.BuiltIn): Promise<void>;
    createComponent(componentId: string, modelId:string, labelAsymId: string, representationType: StructureRepresentationRegistry.BuiltIn): Promise<void>;
    createComponent(componentId: string, modelId:string, residues: Array<{labelAsymId: string; position: number;}>, representationType: StructureRepresentationRegistry.BuiltIn): Promise<void>;
    createComponent(componentId: string, modelId:string, residues: Array<{labelAsymId: string; begin: number; end: number;}>, representationType: StructureRepresentationRegistry.BuiltIn): Promise<void>;
    colorComponent(componentId: string, color: ColorTheme.BuiltIn): Promise<void>;
    removeComponent: (componentId?: string) => void;
    isComponent: (componentId: string) => boolean;
    displayComponent(componentLabel: string, visibilityFlag: boolean): void;
    displayComponent(componentLabel: string): boolean;
    getComponentSet: () => Set<string>;
    setFocus(modelId: string, labelAsymId: string, begin: number, end: number): void;
    clearFocus(): void;
    cameraFocus(modelId: string, labelAsymId: string, positions:Array<number>): void;
    cameraFocus(modelId: string, labelAsymId: string, begin: number, end: number): void;
    resetCamera: ()=>void;
}
