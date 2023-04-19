import * as React from "react";
import * as classes from '../styles/RcsbFvStyle.module.scss';

import {MolstarPlugin} from '../RcsbFvStructure/StructurePlugins/MolstarPlugin';
import {SaguaroPluginInterface} from '../RcsbFvStructure/StructurePlugins/SaguaroPluginInterface';

import '../styles/RcsbFvMolstarStyle.module.scss';
import {RcsbFvSequence, RcsbFvSequenceInterface} from "../RcsbFvSequence/RcsbFvSequence";
import {RcsbFvStructure, RcsbFvStructureInterface} from "../RcsbFvStructure/RcsbFvStructure";
import {
    EventType,
    RcsbFvContextManager,
    RcsbFvContextManagerInterface,
    UpdateConfigInterface
} from "../RcsbFvContextManager/RcsbFvContextManager";
import {Subscription} from "rxjs";
import {PluginContext} from "molstar/lib/mol-plugin/context";
import {RcsbFvSelectorManager} from "../RcsbFvSelection/RcsbFvSelectorManager";
import {CSSProperties, MouseEvent} from "react";

export interface RcsbFv3DCssConfig {
    overwriteCss?: boolean;
    rootPanel?: CSSProperties;
    structurePanel?: CSSProperties;
    sequencePanel?: CSSProperties;
}

export interface RcsbFv3DComponentInterface {
    structurePanelConfig:RcsbFvStructureInterface;
    sequencePanelConfig: RcsbFvSequenceInterface;
    id: string;
    ctxManager: RcsbFvContextManager;
    cssConfig?:RcsbFv3DCssConfig;
    unmount:(flag:boolean)=>void;
    fullScreen: boolean;
}

interface RcsbFv3DComponentState {
    structurePanelConfig:RcsbFvStructureInterface;
    sequencePanelConfig:RcsbFvSequenceInterface;
    pfvScreenFraction: number;
}

export class RcsbFv3DComponent extends React.Component <RcsbFv3DComponentInterface, RcsbFv3DComponentState> {

    private readonly plugin: SaguaroPluginInterface;
    private readonly selectorManager: RcsbFvSelectorManager = new RcsbFvSelectorManager();
    private subscription: Subscription;
    private readonly ROOT_DIV_ID: string = "rootPanelDiv";

    readonly state: RcsbFv3DComponentState = {
        structurePanelConfig: this.props.structurePanelConfig,
        sequencePanelConfig: this.props.sequencePanelConfig,
        pfvScreenFraction: 0.55
    }

    constructor(props: RcsbFv3DComponentInterface) {
        super(props);
        this.plugin = new MolstarPlugin(this.selectorManager);
    }

    render(): JSX.Element {
        return (
            <div className={this.props.fullScreen ? classes.fullScreen : classes.fullHeight} >
                <div
                    id={this.ROOT_DIV_ID}
                    style={RcsbFv3DComponent.mainDivCssConfig(this.props.cssConfig?.rootPanel)}
                    className={this.useDefaultCss() ? classes.rcsbFvMain : ""}
                    onMouseMove={(evt: MouseEvent<HTMLDivElement>)=>{this.mouseMove(evt)}}
                    onMouseUp={ (e)=>{this.splitPanelMouseUp()} }
                >
                    <div style={this.structureCssConfig(this.props.cssConfig?.structurePanel)} >
                        <RcsbFvStructure
                            {...this.state.structurePanelConfig}
                            componentId={this.props.id}
                            plugin={this.plugin}
                            selectorManager={this.selectorManager}
                        />
                    </div>
                    <div style={this.sequenceCssConfig(this.props.cssConfig?.sequencePanel)}  >
                        <RcsbFvSequence
                            type={this.state.sequencePanelConfig.type}
                            config={this.state.sequencePanelConfig.config}
                            componentId={this.props.id}
                            plugin={this.plugin}
                            selectorManager={this.selectorManager}
                            title={this.state.sequencePanelConfig.title}
                            subtitle={this.state.sequencePanelConfig.subtitle}
                            unmount={this.props.unmount}
                        />
                    </div>
                    {
                        this.panelDelimiter()
                    }
                </div>
            </div>
        );
    }

    componentDidMount() {
        this.subscription = this.subscribe();
    }

    componentWillUnmount() {
        this.unsubscribe();
    }

    private useDefaultCss(): boolean {
       return this.state.sequencePanelConfig.type === "assembly"  || !this.props.cssConfig?.overwriteCss;
    }

    private panelDelimiter(): JSX.Element {
        return  this.useDefaultCss() ? <div
            onMouseDown={() => {
                this.splitPanelMouseDown()
            }}
            className={classes.rcsbFvSplitPanel}
            style={{right: Math.round((1 - this.state.pfvScreenFraction) * 100) + "%"}}
        /> : <></>;
    }

    private structureCssConfig(css: CSSProperties | undefined): CSSProperties{
        const widthFr: number = Math.round((1-this.state.pfvScreenFraction)*100);
        const cssWidth: string = widthFr.toString()+"%";
        const cssHeight: string = "100%";
        return {...(this.useDefaultCss() ? {width:cssWidth, height:cssHeight, zIndex:100} : {}), ...css };
    }

    private sequenceCssConfig(css: CSSProperties | undefined): CSSProperties{
        const widthFr: number = Math.round((this.state.pfvScreenFraction)*100);
        const cssWidth: string = widthFr.toString()+"%";
        const cssHeight: string = "100%";
        return {...(this.useDefaultCss() ? {width:cssWidth, height:cssHeight, overflowY:"auto", overflowX:"hidden", paddingBottom:5} : {}), ...css };
    }

    private static mainDivCssConfig(css: CSSProperties | undefined): CSSProperties{
        return {...{

        }, ...css}
    }

    private subscribe(): Subscription{
        return this.props.ctxManager.subscribe((obj:RcsbFvContextManagerInterface)=>{
            if(obj.eventType == EventType.UPDATE_CONFIG){
                this.updateConfig(obj.eventData as UpdateConfigInterface)
            }else if(obj.eventType == EventType.PLUGIN_CALL){
                this.plugin.pluginCall(obj.eventData as ((f:PluginContext)=>void));
            }
        });
    }

    /**Unsubscribe className to rxjs events. Useful if many panels are created an destroyed.*/
    private unsubscribe(): void{
        this.subscription.unsubscribe();
    }

    private updateConfig(config:UpdateConfigInterface){
        const structureConfig: RcsbFvStructureInterface | undefined = config.structurePanelConfig;
        const sequenceConfig: RcsbFvSequenceInterface | undefined = config.sequencePanelConfig;
        if(structureConfig != null && sequenceConfig != null){
            this.setState({structurePanelConfig:structureConfig, sequencePanelConfig:sequenceConfig});
        }else if(structureConfig != null){
            this.setState({structurePanelConfig:structureConfig});
        }else if(sequenceConfig != null){
            this.setState({sequencePanelConfig: sequenceConfig});
        }
    }

    private splitPanelMouseDown(): void {
        const element: HTMLElement | null = document.getElementById(this.ROOT_DIV_ID);
        if(!element)return;
        element.style.cursor = "ew-resize";
        document.body.classList.add(classes.disableTextSelection);
        this.resize = (evt: MouseEvent<HTMLDivElement>)=>{
            const rect: DOMRect | undefined = element.getBoundingClientRect();
            const x: number = evt.clientX - rect.left;
            this.setState({pfvScreenFraction:x/rect.width});
        };
    }

    private splitPanelMouseUp(): void {
        if(typeof this.resize === "function") {
            const element: HTMLElement | null = document.getElementById(this.ROOT_DIV_ID);
            if (!element) return;
            element.style.cursor = "auto";
            document.body.classList.remove(classes.disableTextSelection);
            window.dispatchEvent(new Event('resize'));
            this.resize = null;
        }
    }

    private mouseMove(evt: MouseEvent<HTMLDivElement>): void{
        if(typeof this.resize === "function")
            this.resize(evt);
    }

    private resize: null | ((evt: MouseEvent<HTMLDivElement>)=>void) = null;

}
