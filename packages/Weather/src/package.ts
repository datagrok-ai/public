/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import {WeatherView} from "./WeatherView";

export let _package = new DG.Package();

//name: test
//input: string s
export function test(s: string) {
}

//name: Weather
//tags: app
export function weatherApp() {
    grok.shell.addView(<any>new WeatherView());
}