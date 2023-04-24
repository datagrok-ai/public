import * as React from "react";
import {SaguaroPluginInterface} from "../../../RcsbFvStructure/StructurePlugins/SaguaroPluginInterface";

interface ChainDisplayInterface {
    plugin: SaguaroPluginInterface;
    label: string;
}

interface ChainDisplayState {
    display: 'visible' | 'hidden';
}

export class ChainDisplay extends React.Component<ChainDisplayInterface, ChainDisplayState>{

    readonly state: ChainDisplayState = {
        display: this.props.plugin.displayComponent(this.props.label) ? 'visible' : 'hidden'
    };

    private changeDisplay(): void{
        if(this.state.display === 'visible') {
            this.props.plugin.displayComponent(this.props.label, false);
            this.setState({display: 'hidden'});
        }else{
            this.props.plugin.displayComponent(this.props.label, true);
            this.setState({display: 'visible'});
        }
    }

    render(): JSX.Element{
        return(
                <input style={{marginLeft:5, marginRight:5}} type={'checkbox'} checked={this.state.display === 'visible'} onChange={this.changeDisplay.bind(this)}/>
        );
    }

}
