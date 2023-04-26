import {asyncScheduler} from "rxjs";

import {AbstractView, AbstractViewInterface} from "./AbstractView";
import {
    RcsbFvBoardConfigInterface,
    RcsbFvRowConfigInterface,
    RcsbFv,
    RcsbFvTrackDataElementInterface
} from "@rcsb/rcsb-saguaro";
import * as React from "react";
import {RcsbFvSelectorManager} from "../../RcsbFvSelection/RcsbFvSelectorManager";
import {
    SaguaroPluginModelMapType,
    SaguaroPluginPublicInterface
} from "../../RcsbFvStructure/StructurePlugins/SaguaroPluginInterface";

export type CustomViewStateInterface = Omit<CustomViewInterface, "modelChangeCallback">;

export interface CustomViewInterface {
    blockConfig: FeatureBlockInterface | Array<FeatureBlockInterface>;
    blockSelectorElement?: (blockSelector: BlockSelectorManager) => JSX.Element;
    modelChangeCallback?: (modelMap: SaguaroPluginModelMapType) => CustomViewStateInterface;
    blockChangeCallback?: (plugin: SaguaroPluginPublicInterface, pfvList: Array<RcsbFv>, selection: RcsbFvSelectorManager) => void;
}

export interface FeatureBlockInterface {
    blockId:string;
    blockTitle?: string;
    blockShortName?: string;
    featureViewConfig: Array<FeatureViewInterface> | FeatureViewInterface;
}

export interface FeatureViewInterface {
    boardId?:string;
    boardConfig: RcsbFvBoardConfigInterface;
    rowConfig: Array<RcsbFvRowConfigInterface>;
    sequenceSelectionChangeCallback: (plugin: SaguaroPluginPublicInterface, selectorManager: RcsbFvSelectorManager, sequenceRegion: Array<RcsbFvTrackDataElementInterface>) => void;
    sequenceElementClickCallback: (plugin: SaguaroPluginPublicInterface, selectorManager: RcsbFvSelectorManager, d: RcsbFvTrackDataElementInterface) => void;
    sequenceHoverCallback: (plugin: SaguaroPluginPublicInterface, selectorManager: RcsbFvSelectorManager, hoverRegion: Array<RcsbFvTrackDataElementInterface>) => void;
    structureSelectionCallback: (plugin: SaguaroPluginPublicInterface, pfv: RcsbFv, selectorManager: RcsbFvSelectorManager) => void;
    structureHoverCallback: (plugin: SaguaroPluginPublicInterface, pfv: RcsbFv, selectorManager: RcsbFvSelectorManager) => void;
}

export class BlockSelectorManager {
    private blockId: string;
    private previousBlockId: string;
    private readonly blockChangeCallback: ()=>void = ()=>{};
    constructor(f:()=>void){
        this.blockChangeCallback = f;
    }
    setActiveBlock(blockId:string): void{
        this.previousBlockId = this.blockId;
        this.blockId = blockId;
        this.blockChangeCallback();
    }
    getActiveBlock(): string{
        return this.blockId;
    }
    getPreviousBlock(): string{
        return this.previousBlockId;
    }
}

export class CustomView extends AbstractView<CustomViewInterface & AbstractViewInterface, CustomViewStateInterface> {

    private blockViewSelector: BlockSelectorManager = new BlockSelectorManager( this.blockChange.bind(this) );
    private boardMap: Map<string, FeatureViewInterface> = new Map<string, FeatureViewInterface>();
    private blockMap: Map<string, Array<string>> = new Map<string, Array<string>>();
    private rcsbFvMap: Map<string, RcsbFv> = new Map<string, RcsbFv>();
    private firstModelLoad: boolean = true;
    private innerSelectionFlag: boolean = false;

    readonly state: CustomViewStateInterface = {
        blockConfig: this.props.blockConfig,
        blockSelectorElement: this.props.blockSelectorElement,
        blockChangeCallback: this.props.blockChangeCallback
    };

    constructor(props: CustomViewInterface & AbstractViewInterface) {
        super(props);
        this.mapBlocks(props.blockConfig);
    }

    componentDidMount(): void {
        super.componentDidMount();
        this.blockViewSelector.setActiveBlock( (this.state.blockConfig instanceof Array ? this.state.blockConfig : [this.state.blockConfig])[0].blockId! );
    }

    componentWillUnmount() {
        super.componentWillUnmount();
        this.rcsbFvMap.forEach((pfv,id)=>{
            pfv.unmount();
        });
    }

    private mapBlocks(config: FeatureBlockInterface | Array<FeatureBlockInterface>){
        this.blockMap.clear();
        this.boardMap.clear();
        ( config instanceof Array ? config : [config]).forEach(block=>{
            if(block.blockId == null)block.blockId = "block_"+Math.random().toString(36).substr(2);
            if(!this.blockMap.has(block.blockId))this.blockMap.set(block.blockId, new Array<string>());
            (block.featureViewConfig instanceof Array ? block.featureViewConfig : [block.featureViewConfig]).forEach(board=>{
                if(board.boardId == null)board.boardId = "board_"+Math.random().toString(36).substr(2);
                this.blockMap.get(block.blockId!)?.push(board.boardId);
                this.boardMap.set(board.boardId, board);
            });
        });
    }

    private blockChange(): void{
        this.unmountBlockFv();
        this.buildBlockFv();
        asyncScheduler.schedule(()=>{
            if(typeof this.state.blockChangeCallback === "function")
                this.state.blockChangeCallback(this.props.plugin, Array.from(this.blockMap.get(this.blockViewSelector.getActiveBlock())!.values()).map(boardId=>(this.rcsbFvMap.get(boardId)!)), this.props.selectorManager);
            else
                this.structureSelectionCallback();
        },1000);
    }

    private unmountBlockFv(){
        this.blockMap.get(this.blockViewSelector.getPreviousBlock())?.forEach(boardId=>{
            if(this.rcsbFvMap.get(boardId) == null)
                return;
            this.rcsbFvMap.get(boardId)!.unmount();
            document.getElementById("boardDiv_"+boardId)?.remove()
        });
        this.rcsbFvMap.clear();
        this.props.plugin.unsetCallbacks();
    }

    private buildBlockFv(){
        this.blockMap.get(this.blockViewSelector.getActiveBlock())?.forEach(boardId=>{
            if(this.boardMap.get(boardId) == null)
                return;
            const div: HTMLDivElement = document.createElement<"div">("div");
            div.setAttribute("id", "boardDiv_"+boardId);
            document.getElementById(this.componentDivId)?.append(div);
            const width: number = window.document.getElementById(this.componentDivId)?.getBoundingClientRect().width ?? 0;
            const trackWidth: number = width - (this.boardMap.get(boardId)!.boardConfig?.rowTitleWidth ?? 190) - 55;
            const rcsbFv: RcsbFv = new RcsbFv({
                elementId: "boardDiv_"+boardId,
                boardConfigData:{
                    highlightHoverPosition:true,
                    highlightHoverElement:true,
                    ...this.boardMap.get(boardId)!.boardConfig,
                    trackWidth:this.boardMap.get(boardId)!.boardConfig?.trackWidth ? this.boardMap.get(boardId)!.boardConfig?.trackWidth!-4 : trackWidth,
                    selectionChangeCallBack:(selection: RcsbFvTrackDataElementInterface[])=>{
                        if(this.innerSelectionFlag)
                            return;
                        this.boardMap.get(boardId)!.sequenceSelectionChangeCallback(this.props.plugin, this.props.selectorManager, selection);
                    },
                    highlightHoverCallback:(elements:Array<RcsbFvTrackDataElementInterface>)=>{
                        this.boardMap.get(boardId)!.sequenceHoverCallback(this.props.plugin, this.props.selectorManager, elements);
                    },
                    elementClickCallBack: (d?: RcsbFvTrackDataElementInterface, e?: MouseEvent)=>{
                        this.boardMap.get(boardId)!.sequenceElementClickCallback(this.props.plugin, this.props.selectorManager, d!);
                    }
                },
                rowConfigData: this.boardMap.get(boardId)!.rowConfig
            });
            this.rcsbFvMap.set(boardId, rcsbFv);
        });
        this.props.plugin.setSelectCallback(()=>{
           this.structureSelectionCallback();
        });
    }

    structureSelectionCallback(): void {
        this.innerSelectionFlag = true;
        this.blockMap.get(this.blockViewSelector.getActiveBlock())?.forEach(boardId=>{
            const pfv: RcsbFv | undefined = this.rcsbFvMap.get(boardId);
            if(pfv == null)
                return;
            this.boardMap.get(boardId)?.structureSelectionCallback(this.props.plugin, pfv, this.props.selectorManager);
        });
        this.innerSelectionFlag = false;
    }

    structureHoverCallback(): void{
        this.blockMap.get(this.blockViewSelector.getActiveBlock())?.forEach(boardId=>{
            const pfv: RcsbFv | undefined = this.rcsbFvMap.get(boardId);
            if(pfv == null)
                return;
            this.boardMap.get(boardId)?.structureHoverCallback(this.props.plugin, pfv, this.props.selectorManager);
        });
    }

    representationChangeCallback(): void{
        //TODO
    }

    additionalContent(): JSX.Element {
        if(this.state.blockSelectorElement == null)
            return <></>;
        return this.state.blockSelectorElement(this.blockViewSelector);
    }

    modelChangeCallback(modelMap:SaguaroPluginModelMapType): void {
        if(this.firstModelLoad){
            this.firstModelLoad = false;
            return;
        }
        if(typeof this.props.modelChangeCallback === "function") {
            let newConfig: CustomViewStateInterface = this.props.modelChangeCallback(modelMap);
            if(newConfig != null ){
                newConfig = newConfig as CustomViewStateInterface;
                if(newConfig.blockConfig != null && newConfig.blockSelectorElement != null){
                    this.mapBlocks(newConfig.blockConfig);
                    this.setState({blockConfig: newConfig.blockConfig, blockSelectorElement: newConfig.blockSelectorElement})
                }else if(newConfig.blockConfig == null && newConfig.blockSelectorElement != null){
                    this.setState({blockSelectorElement: newConfig.blockSelectorElement})
                }else if(newConfig.blockConfig != null && newConfig.blockSelectorElement == null){
                    this.mapBlocks(newConfig.blockConfig);
                    this.setState({blockConfig: newConfig.blockConfig})
                }
            }
        }
    }

    updateDimensions(): void {
        const div: HTMLElement | undefined | null = document.getElementById(this.componentDivId)?.parentElement;
        const width: number = window.document.getElementById(this.componentDivId)?.getBoundingClientRect().width ?? 0;
        if(div == null || (div.style.width && !div.style.width.includes("%")) )
            return;
        this.rcsbFvMap.forEach((rcsbFv, boardId)=>{
            const trackWidth: number = width - (this.boardMap.get(boardId)!.boardConfig?.rowTitleWidth ?? 190) - 55;
            rcsbFv.updateBoardConfig({boardConfigData:{trackWidth:this.boardMap.get(boardId)!.boardConfig?.trackWidth ? this.boardMap.get(boardId)!.boardConfig?.trackWidth!-4 : trackWidth}});
        });
    }

}
