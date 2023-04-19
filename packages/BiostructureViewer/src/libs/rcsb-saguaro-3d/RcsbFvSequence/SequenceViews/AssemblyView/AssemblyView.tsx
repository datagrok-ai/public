import {asyncScheduler} from "rxjs";
import * as React from "react";

import {RcsbFvDOMConstants} from "../../../RcsbFvConstants/RcsbFvConstants";
import {
    buildInstanceSequenceFv,
    unmount
} from "@rcsb/rcsb-saguaro-app";
import {AbstractView, AbstractViewInterface} from "../AbstractView";
import {InstanceSequenceOnchangeInterface} from "@rcsb/rcsb-saguaro-app/build/dist/RcsbFvWeb/RcsbFvBuilder/RcsbFvInstanceBuilder";
import {RcsbFvBoardConfigInterface, RcsbFvTrackDataElementInterface} from "@rcsb/rcsb-saguaro";
import {ChainSelectionInterface} from "../../../RcsbFvSelection/RcsbFvSelectorManager";
import {
    SaguaroPluginInterface,
    SaguaroPluginModelMapType
} from "../../../RcsbFvStructure/StructurePlugins/SaguaroPluginInterface";
import {OptionPropsInterface} from "@rcsb/rcsb-saguaro-app/build/dist/RcsbFvWeb/WebTools/SelectButton";

import {OptionProps} from "react-select/src/components/Option";
import {components} from 'react-select';
import {ChainDisplay} from "./ChainDisplay";

import {
    RcsbFvAdditionalConfig,
    RcsbFvModulePublicInterface
} from "@rcsb/rcsb-saguaro-app/build/dist/RcsbFvWeb/RcsbFvModule/RcsbFvModuleInterface";

export interface AssemblyViewInterface {
    entryId: string;
    additionalConfig?: RcsbFvAdditionalConfig;
}

export class AssemblyView extends AbstractView<AssemblyViewInterface & AbstractViewInterface, {}>{

    private currentLabelAsymId: string;
    private currentEntryId: string;
    private currentModelId: string;
    private currentModelNumber: string;
    private createComponentThreshold: number = 3;
    private innerSelectionFlag: boolean = false;
    private currentSelectedComponentId: string;
    private currentModelMap:SaguaroPluginModelMapType;
    private boardConfig: Partial<RcsbFvBoardConfigInterface>;
    private rcsbFvModule: RcsbFvModulePublicInterface | null;
    //private readonly componentSet = new Map<string, {current: Set<string>, previous: Set<string>}>();

    additionalContent(): JSX.Element {
        return (
            <div style={{marginTop:10}}>
                <div>
                    <div id={RcsbFvDOMConstants.SELECT_INSTANCE_PFV_ID} style={{display:"inline-block"}}/>
                    <div style={{display:"inline-block", marginLeft:25}}>
                        <a href={"/docs/sequence-viewers/protein-feature-view"} target={"_blank"}>Help</a>
                    </div>
                </div>
                <div style={{position:"absolute", top:5, right:5}} >
                    <a style={{textDecoration:"none", color:"#337ab7", cursor:"pointer", marginRight:15}} target={"_blank"} href={"/docs/sequence-viewers/3d-protein-feature-view"}>
                        Help
                    </a>
                    <a style={{textDecoration:"none", color: "#337ab7", cursor:"pointer"}} onClick={()=>{this.props.unmount(true)}}>
                        Back
                    </a>
                </div>
            </div>
        );
    }

    componentDidMount (): void {
        super.componentDidMount();
        const width: number | undefined = document.getElementById(this.componentDivId)?.getBoundingClientRect().width;
        if(width == null)
            return;
        const trackWidth: number = width - 190 - 55;
        this.boardConfig = {
            ...this.props.additionalConfig?.boardConfig,
            trackWidth: trackWidth,
            highlightHoverPosition:true,
            highlightHoverElement:true,
            elementClickCallBack:(d?:RcsbFvTrackDataElementInterface, e?: MouseEvent)=>{
                this.elementClickCallback(d!);
                if(typeof this.props.additionalConfig?.boardConfig?.elementClickCallBack === "function")
                    this.props.additionalConfig?.boardConfig.elementClickCallBack(d);
            },
            selectionChangeCallBack:(selection: Array<RcsbFvTrackDataElementInterface>)=>{
                this.selectionChangeCallback(selection);
                if(typeof this.props.additionalConfig?.boardConfig?.selectionChangeCallBack === "function")
                    this.props.additionalConfig?.boardConfig.selectionChangeCallBack(selection);
            },
            highlightHoverCallback:(selection: RcsbFvTrackDataElementInterface[])=>{
                this.highlightHoverCallback(selection);
                if(typeof this.props.additionalConfig?.boardConfig?.highlightHoverCallback === "function")
                    this.props.additionalConfig?.boardConfig.highlightHoverCallback(selection);
            },
        };
    }

    componentWillUnmount() {
        super.componentWillUnmount();
        unmount(this.rcsbFvDivId);
    }

    async structureSelectionCallback(): Promise<void> {
        await this.pluginSelectCallback('select');
    }

    async structureHoverCallback(): Promise<void> {
        await this.pluginSelectCallback('hover');
    }

    representationChangeCallback(): void{
        //TODO
    }

    async updateDimensions(): Promise<void> {
        const width: number = window.document.getElementById(this.componentDivId)?.getBoundingClientRect().width ?? 0;
        const trackWidth: number = width - 190 - 55;
        this.boardConfig.trackWidth = trackWidth;
        await this.rcsbFvModule?.getFv().updateBoardConfig({boardConfigData:{trackWidth:trackWidth}})
        await this.structureSelectionCallback();
        return void 0;
    }

    private resetPluginView(): void {
        this.props.plugin.clearFocus();
        this.props.plugin.resetCamera();
    }

    private async pluginSelectCallback(mode:'select'|'hover'): Promise<void> {
        if(this.rcsbFvModule == null)
            return;
        this.innerSelectionFlag = true;
        if(mode === 'select' && this.currentSelectedComponentId != null){
            this.props.plugin.removeComponent(this.currentSelectedComponentId);
        }
        const allSel: Array<ChainSelectionInterface> | undefined = this.props.selectorManager.getSelection(mode);
        if(allSel == null || allSel.length ===0) {
            this.rcsbFvModule?.getFv().clearSelection(mode);
            if(mode === 'select')
                this.resetPluginView();
        }else if(mode === 'select' && this.props.selectorManager.getLastSelection('select')?.labelAsymId != null && this.props.selectorManager.getLastSelection('select')?.labelAsymId != this.currentLabelAsymId){
            const authId: string | undefined = this.currentModelMap
                .get(this.currentModelId)?.chains
                .filter(ch=>(ch.label===this.props.selectorManager.getLastSelection('select')?.labelAsymId))[0]?.auth;
            await this.modelChangeCallback(this.currentModelMap, authId);
        }else{
            const sel: ChainSelectionInterface | undefined = this.props.selectorManager.getSelectionWithCondition(this.currentModelId, this.currentLabelAsymId, mode);
            if (sel == null) {
                this.rcsbFvModule?.getFv().clearSelection(mode);
                if(mode === 'select')
                    this.resetPluginView();
            } else {
                this.rcsbFvModule?.getFv().setSelection({elements: sel.regions, mode: mode});
            }
        }
        this.innerSelectionFlag = false;
    }

    async modelChangeCallback(modelMap:SaguaroPluginModelMapType, defaultAuthId?: string): Promise<void> {
        this.currentModelMap = modelMap;
        this.props.plugin.clearFocus();
        const onChangeCallback: Map<string, (x: InstanceSequenceOnchangeInterface)=>void> = new Map<string, (x: InstanceSequenceOnchangeInterface) => {}>();
        const filterInstances: Map<string, Set<string>> = new Map<string, Set<string>>();
        modelMap.forEach((v,k)=>{
            onChangeCallback.set(v.entryId,(x)=>{
                this.currentEntryId = v.entryId;
                this.currentLabelAsymId = x.asymId;
                this.currentModelId = k;
                asyncScheduler.schedule(()=>{
                    this.props.selectorManager.setLastSelection('select', null);
                    this.structureSelectionCallback();
                },1000);
            });
            filterInstances.set(v.entryId,new Set<string>(v.chains.map(d=>d.label)));
        });
        this.unmountRcsbFv();
        const entryId: string = Array.from(modelMap.values()).map(d=>d.entryId)[0];
        if(entryId != null) {
            this.rcsbFvModule = await buildInstanceSequenceFv(
                this.rcsbFvDivId,
                RcsbFvDOMConstants.SELECT_INSTANCE_PFV_ID,
                entryId,
                {
                    defaultValue: defaultAuthId,
                    onChangeCallback: onChangeCallback.get(entryId),
                    filterInstances: filterInstances.get(entryId),
                    selectButtonOptionProps: (props: OptionProps<OptionPropsInterface>) => (components.Option &&
                        <div style={{display: 'flex'}}>
                            <ChainDisplay plugin={this.props.plugin} label={props.data.label}/>
                            <components.Option {...props}/>
                        </div>)
                },
                {
                    ...this.props.additionalConfig,
                    boardConfig: this.boardConfig
                }
            );
        }
        if(!defaultAuthId)
            await createComponents(this.props.plugin, modelMap);
    }

    private unmountRcsbFv(): void {
        this.rcsbFvModule = null;
        unmount(this.rcsbFvDivId);
    }

    private highlightHoverCallback(selection: RcsbFvTrackDataElementInterface[]): void {
        if(selection != null && selection.length > 0) {
            if(selection[0].isEmpty){
                const selectionList = [{modelId: this.currentModelId, labelAsymId: this.currentLabelAsymId, position: selection[0].begin}];
                if(selection[0].end != null) selectionList.push({modelId: this.currentModelId, labelAsymId: this.currentLabelAsymId, position: selection[0].end})
                this.props.plugin.select(
                    selectionList,
                    'hover',
                    'set'
                );
            }else {
                this.props.plugin.select(processMultipleGaps(this.currentModelId, this.currentLabelAsymId, selection), 'hover', 'set');
            }
        }else{
            this.props.plugin.clearSelection('hover');
        }
    }

    private selectionChangeCallback(selection: Array<RcsbFvTrackDataElementInterface>): void {
        if(this.innerSelectionFlag)
            return;
        this.props.plugin.clearSelection('select', {modelId: this.currentModelId, labelAsymId: this.currentLabelAsymId});
        this.props.selectorManager.clearSelection('select', {labelAsymId: this.currentLabelAsymId});
        if(selection == null || selection.length === 0) {
            this.resetPluginView();
        }else{
            this.select(selection);
        }
    }

    private select(selection: Array<RcsbFvTrackDataElementInterface>): void{
        selection.forEach(e=>{
            const x = e.begin;
            const y = e.end ?? e.begin;
            if(e.isEmpty){
                this.props.plugin.select(
                    [{modelId: this.currentModelId, labelAsymId: this.currentLabelAsymId, position: x},{modelId: this.currentModelId, labelAsymId: this.currentLabelAsymId, position: y}], 'select',
                    'add'
                );
                this.props.selectorManager.addSelectionFromRegion(this.currentModelId, this.currentLabelAsymId, {begin:x, end:y, isEmpty: true, source: 'sequence'}, 'select');
            }else{
                this.props.plugin.select(processGaps(this.currentModelId, this.currentLabelAsymId, e), 'select', 'add');
                this.props.selectorManager.addSelectionFromRegion(this.currentModelId, this.currentLabelAsymId, {begin:x, end:y, source: 'sequence'}, 'select');
            }
        });
    }

    private elementClickCallback(e:RcsbFvTrackDataElementInterface): void {
        this.props.plugin.clearFocus();
        if(this.currentSelectedComponentId != null)
            this.props.plugin.removeComponent(this.currentSelectedComponentId);
        if(e == null)
            return;
        const x = e.begin;
        const y = e.end ?? e.begin;
        if(e.isEmpty){
            this.props.plugin.cameraFocus(this.currentModelId, this.currentLabelAsymId, [x,y]);
            this.currentSelectedComponentId = this.currentLabelAsymId +":"+ ((x === y) ? x.toString() : x.toString()+","+y.toString());
            asyncScheduler.schedule(async ()=>{
                await this.props.plugin.createComponent(
                    this.currentSelectedComponentId,
                    this.currentModelId,
                    [{labelAsymId: this.currentLabelAsymId, position: x}, {labelAsymId: this.currentLabelAsymId, position: y}],
                    'ball-and-stick'
                )
                if(x === y)
                    asyncScheduler.schedule(()=>{
                        this.props.plugin.setFocus(this.currentModelId, this.currentLabelAsymId, x, y);
                    },200);
            },100);

        }else{
            this.props.plugin.cameraFocus(this.currentModelId, this.currentLabelAsymId, x, y);
            if((y-x)<this.createComponentThreshold){
                this.currentSelectedComponentId = this.currentLabelAsymId +":"+ (x === y ? x.toString() : x.toString()+"-"+y.toString());
                asyncScheduler.schedule(async ()=>{
                    await this.props.plugin.createComponent(
                        this.currentSelectedComponentId,
                        this.currentModelId,
                        processGaps(this.currentModelId, this.currentLabelAsymId, e),
                        'ball-and-stick'
                    )
                    if(x === y)
                        asyncScheduler.schedule(()=>{
                            this.props.plugin.setFocus(this.currentModelId, this.currentLabelAsymId, x, y);
                        },200);
                },100);
            }
        }
    }

}

function processGaps(modelId: string, labelAsymId: string, e: RcsbFvTrackDataElementInterface): Array<{modelId: string; labelAsymId: string; begin: number; end: number;}>{
    const regions: Array<{modelId: string; labelAsymId: string; begin: number; end: number;}> = new Array<{modelId: string; labelAsymId: string; begin: number; end: number}>();
    let lastIndex: number = e.begin;
    e.gaps?.forEach((g)=>{
        regions.push({
            modelId: modelId,
            labelAsymId: labelAsymId,
            begin: lastIndex,
            end: g.begin
        });
        lastIndex = g.end;
    });
    regions.push({
        modelId: modelId,
        labelAsymId: labelAsymId,
        begin: lastIndex,
        end: e.end ?? e.begin
    });
    return regions;
}

function processMultipleGaps(modelId: string, labelAsymId: string, list: Array<RcsbFvTrackDataElementInterface>): Array<{modelId: string; labelAsymId: string; begin: number; end: number;}>{
    let regions: Array<{modelId: string; labelAsymId: string; begin: number; end: number;}> = new Array<{modelId: string; labelAsymId: string; begin: number; end: number}>();
    list.forEach(e=>{
        regions = regions.concat(processGaps(modelId, labelAsymId, e));
    });
    return regions;
}

async function createComponents(plugin: SaguaroPluginInterface, modelMap:SaguaroPluginModelMapType): Promise<void> {
    plugin.displayComponent("Water", false);
    await plugin.colorComponent("Polymer", 'chain-id');
    const chains: Array<{modelId: string; auth: string; label: string;}> = new Array<{modelId: string; auth: string; label: string;}>();
    modelMap.forEach((entry, modelId)=>{
        entry.chains.forEach(ch=>{
            if(ch.type === "polymer") {
                chains.push({modelId: modelId, auth: ch.auth, label: ch.label});
            }
        });
    });
    plugin.removeComponent();
    plugin.clearFocus();
    for(const ch of chains) {
        const label: string = ch.auth === ch.label ? ch.label : `${ch.label} [auth ${ch.auth}]`;
        await plugin.createComponent(label, ch.modelId, ch.label, 'cartoon');
        await plugin.colorComponent(label, 'chain-id');
    }
    await plugin.removeComponent("Polymer");
}
