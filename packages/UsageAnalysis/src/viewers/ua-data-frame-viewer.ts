import {UaViewer} from "./ua-viewer";
import {UaFilter} from "../filter2";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";

export class UaDataFrameViewer extends UaViewer {
    queryName: string;
    viewerFunction: Function;
    setStyle: Function = null as any;

    static splineStyle: Object = {
        "aggrType": "count",
        "innerChartMarginTop": 0,
        "innerChartMarginBottom": 0,
        "outerChartMarginTop": 5,
        "outerChartMarginBottom": 0,
        "yGlobalScale": false,
        "showTopPanel": false,
        "showMouseOverRowLine": false,
        "showXSelector": false,
        "showYSelectors": false,
        "showAggrSelectors": false,
        "showSplitSelector": false,
        "showYAxis": false,
        "showMarkers": "Never",
        "Title":"Users"
    };


    public constructor(name: string, queryName: string, viewerFunction: Function, setStyle?: Function) {
        super(name);
        this.queryName = queryName;
        this.viewerFunction = viewerFunction;
        if (setStyle)
            this.setStyle = setStyle;
        this.init();
    }


    // addCardUsingDataframe(cardName: string, dataFrame: DG.DataFrame, viewer:any, supportUsers = true) {
    //     let host = ui.block([],'d4-item-card card');
    //     host.appendChild(ui.h1(cardName));
    //
    //     // if (cardName === 'Errors')
    //     //     grok.data.detectSemanticTypes(dataFrame);
    //     host.appendChild(viewer(dataFrame));
    //     return host;
    // }

    addCardWithFilters(cardName: string, queryName: string, filter: UaFilter, viewer:any) {
        let host = ui.block([]);
        if (this.setStyle)
            this.setStyle(host);
        host.appendChild(ui.h1(cardName));
        let loader = ui.loader();
        host.appendChild(loader);

        grok.data.query('UsageAnalysis:' + queryName, filter).then((dataFrame) => {
            // if (cardName === 'Errors')
            //     grok.data.detectSemanticTypes(dataFrame);
            host.appendChild(viewer(dataFrame));
            host.removeChild(loader);
        });
        return host;
    }

    reload(filter: UaFilter): void {
    }

    init(): void {
        this.root.innerHTML = '';
        this.root.append(
            this.addCardWithFilters(
                this.name,
                this.queryName,
                new UaFilter(),
                this.viewerFunction
            )
        );
    }

}