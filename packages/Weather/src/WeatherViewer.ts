import * as DG from "datagrok-api/dg";
import * as ui from 'datagrok-api/ui';

export class WeatherViewer extends DG.JsViewer {

    country: string = this.string('country', 'France');

    constructor() {
        super();

        this.root.appendChild(ui.button('Spain', (_: any) => this.props.country = 'Spain'));
    }
}