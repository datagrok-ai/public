import * as DG from "datagrok-api/dg";
import * as ui from 'datagrok-api/ui';

export class WeatherView extends DG.ViewBase {
    constructor() {
        super({});
        this.root.appendChild(ui.divText('Hi there!'));
    }
}