import * as React from "react";
import {SaguaroPluginInterface} from "./StructurePlugins/SaguaroPluginInterface";
import {RcsbFvDOMConstants} from "../RcsbFvConstants/RcsbFvConstants";
import {ViewerProps} from "@rcsb/rcsb-molstar/build/src/viewer";
import {LoadMolstarInterface} from "./StructurePlugins/MolstarPlugin";
import {RcsbFvSelectorManager} from "../RcsbFvSelection/RcsbFvSelectorManager";

export interface RcsbFvStructureInterface {
    loadConfig: LoadMolstarInterface;
    pluginConfig?: Partial<ViewerProps>;
}

export class RcsbFvStructure extends React.Component <RcsbFvStructureInterface & {plugin: SaguaroPluginInterface, componentId: string, selectorManager: RcsbFvSelectorManager}, RcsbFvStructureInterface > {

    render():JSX.Element {
        return (
            <div id={this.props.componentId+"_"+RcsbFvDOMConstants.MOLSTAR_DIV} >
                <div id={this.props.componentId+"_"+RcsbFvDOMConstants.MOLSTAR_APP_ID} style={{position:"absolute"}}/>
            </div>
        );
    }

    componentDidMount() {
        this.updateDimensions();
        this.props.plugin.init(this.props.componentId+"_"+RcsbFvDOMConstants.MOLSTAR_APP_ID, this.props.pluginConfig);
        this.props.plugin.load(this.props.loadConfig);
        window.addEventListener('resize', this.updateDimensions.bind(this));
    }

    private updateDimensions(): void {
        const div: HTMLElement | undefined | null = document.getElementById(this.props.componentId+"_"+RcsbFvDOMConstants.MOLSTAR_DIV)?.parentElement;
        if(div == null)
            return;
        const rect: DOMRect = div.getBoundingClientRect()
        RcsbFvStructure.setSize(document.getElementById(this.props.componentId+"_"+RcsbFvDOMConstants.MOLSTAR_DIV), rect);
        RcsbFvStructure.setSize(document.getElementById(this.props.componentId+"_"+RcsbFvDOMConstants.MOLSTAR_APP_ID), rect);
    }

    private static setSize(element: HTMLElement | null, rect: DOMRect | undefined): void{
        if(element == null)
            return;
        if(rect == null)
            return;
        element.style.width = rect.width+"px";
        element.style.height = rect.height+"px";
    }

}
