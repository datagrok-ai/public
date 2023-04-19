import * as React from "react";
import {AssemblyView, AssemblyViewInterface} from "./SequenceViews/AssemblyView/AssemblyView";
import {CustomView, CustomViewInterface} from "./SequenceViews/CustomView";
import {SaguaroPluginInterface} from "../RcsbFvStructure/StructurePlugins/SaguaroPluginInterface";
import {PluginContext} from "molstar/lib/mol-plugin/context";
import {RcsbFv, RcsbFvTrackDataElementInterface} from "@rcsb/rcsb-saguaro";
import {RcsbFvSelectorManager} from "../RcsbFvSelection/RcsbFvSelectorManager";

export interface RcsbFvSequenceInterface{
    type: "custom" | "assembly";
    config: AssemblyViewInterface | CustomViewInterface;
    title?: string;
    subtitle?: string;
}

interface CallbackConfig {
    structureCallback?: (plugin: PluginContext, ann: RcsbFvTrackDataElementInterface)=>void;
    sequenceCallback?: (rcsbFv: RcsbFv)=>void;
}

export class RcsbFvSequence extends React.Component <RcsbFvSequenceInterface & CallbackConfig & {unmount:(flag:boolean)=>void, plugin: SaguaroPluginInterface, selectorManager:RcsbFvSelectorManager, componentId:string}, RcsbFvSequenceInterface > {

    render() {
        if(this.props.type == "custom"){
            const config: CustomViewInterface = this.props.config as CustomViewInterface;
            return (<CustomView
                {...config}
                componentId={this.props.componentId}
                plugin={this.props.plugin}
                selectorManager={this.props.selectorManager}
                title={this.props.title}
                subtitle={this.props.subtitle}
                unmount={this.props.unmount}
            />)
        }else if(this.props.type == "assembly"){
            const config: AssemblyViewInterface = this.props.config as AssemblyViewInterface;
            return (<AssemblyView
                {...config}
                componentId={this.props.componentId}
                plugin={this.props.plugin}
                selectorManager={this.props.selectorManager}
                title={this.props.title}
                subtitle={this.props.subtitle}
                unmount={this.props.unmount}
            />)
        }
    }

}
