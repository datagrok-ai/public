import {VIEWER} from './const';
import {Column, DataFrame} from "./dataframe";
import {StreamSubscription} from "./events";
import {DateTime, Property} from "./entities";
import {Subscription} from "rxjs";

/**
 * Represents a {@link https://datagrok.ai/help/visualize/viewers | viewer}.
 See also {@link https://datagrok.ai/help/develop/js-api#pre-defined-viewers}

 Use Viewer to control the viewers. To develop a custom viewer, {@see JsViewer}.
 * @example
 * let view = grok.shell.addTableView(grok.data.testData('demog', 5000));
 view.scatterPlot({
 *      x: 'height',
 *      y: 'weight',
 *      size: 'age',
 *      color: 'race',
 *    });
 */
export class Viewer {
    /**
     * Visual root.
     */
    root: HTMLElement;

    /**
     * Creates a new viewer of the specified type.
     */
    static fromType(viewerType: VIEWER | string, table: DataFrame, options?: object | null): Viewer;

    static grid(t: DataFrame, options?: object | null): Viewer;

    static histogram(t: DataFrame, options?: object | null): Viewer;

    static barchart(t: DataFrame, options?: object | null): Viewer;

    static boxPlot(t: DataFrame, options?: object | null): Viewer;

    static filters(t: DataFrame, options?: object | null): Viewer;

    static scatterPlot(t: DataFrame, options?: object | null): Viewer;

    /**
     * Sets viewer options. See also {@link getOptions}
     Sample: https://public.datagrok.ai/js/samples/ui/viewers/scatter-plot
     */
    setOptions(map: object): void;

    /**
     * Gets viewer options. See also {@link getOptions}
     Sample: https://public.datagrok.ai/js/samples/ui/viewers/scatter-plot
     */
    getOptions(): object;

    /**
     * Closes and detaches the viewer.
     */
    close(): void;
}


export class JsViewerProps {

}


/**
 * Subclass JsViewer to implement a DataFrame-bound Datagrok viewer in JavaScript.
 See an example on github: {@link https://github.com/datagrok-api/public/tree/master/packages/Leaflet}
 */
export class JsViewer {
    root: HTMLElement;
    properties: Property[];
    dataFrame: DataFrame;
    subs: Subscription[];
    props: any;

    /**
     * Gets called when a table is attached to the viewer.
     */
    onTableAttached(dataFromHandle: any): void;

    /**
     * Gets called when viewer's property is changed.
     * @param property - or null, if multiple properties were changed.
     */
    onPropertyChanged(property: Property): void;

    /**
     * Gets called when this viewer is detached.
     */
    detach(): void;

    /**
     * Gets property ba name (case-sensitive).
     */
    getProperty(name: string): Property;

    /**
     * cleanup() will get called when the viewer is disposed
     */
    registerCleanup(cleanup: (...params: any[]) => any): void;

    /**
     * Returns the column bound to the specified data property.
     Note that "ColumnName" suffix (this determines whether this is a data property) should be omitted.
     */
    column(dataPropertyName: string): Column;

    /**
     * Registers an integer property with the specified name and defaultValue
     */
    int(propertyName: string, defaultValue?: number | null): number;

    /**
     * Registers a floating point property with the specified name and defaultValue
     */
    float(propertyName: string, defaultValue?: number | null): number;

    /**
     * Registers a string property with the specified name and defaultValue
     */
    string(propertyName: string, defaultValue?: string | null): string;

    /**
     * Registers a string list property with the specified name and defaultValue
     */
    stringList(propertyName: string, defaultValue?: string[] | null): string[];

    /**
     * Registers a boolean property with the specified name and defaultValue
     */
    bool(propertyName: string, defaultValue?: boolean | null): boolean;

    /**
     * Registers a datetime property with the specified name and defaultValue
     */
    dateTime(propertyName: string, defaultValue?: DateTime | null): DateTime;
}