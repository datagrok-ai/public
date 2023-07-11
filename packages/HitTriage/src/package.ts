/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from './app/hit-triage-app';

export const _package = new DG.Package();

//tags: app
//name: Hit Triage
export async function hitTriageApp() {
  new HitTriageApp();
}

//name: demosmth
//tags: HitTriageFunction
export async function demosmth() {
  return 1;
};

//name: Demo Molecules 100
//tags: HitTriageDataSource
//output: dataframe result
export async function demoFileIngest(): Promise<DG.DataFrame> {
  return grok.data.demo.molecules(100);
}

//name: Demo Molecules 5000
//tags: HitTriageDataSource
//output: dataframe result
export async function demoFileIngest1(): Promise<DG.DataFrame> {
  return grok.data.demo.molecules(5000);
}

//name: Demo File Submit
//tags: HitTriageSubmitFunction
//input: dataframe df [Dataframe]
//input: string molecules [Molecules column name]
export async function demoFileSubmit(df: DG.DataFrame, molecules: string): Promise<void> {
  grok.shell.info(df.name);
  grok.shell.info(molecules);
}

//name: Sample File Submit
//tags: HitTriageSubmitFunction
//input: dataframe df [dataframe]
//input: string molecules [molecules column name]
export async function demoFileSubmit1(df: DG.DataFrame, molecules: string): Promise<void> {
  grok.shell.info('Another function');
}
