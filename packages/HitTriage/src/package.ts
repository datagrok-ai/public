/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from './app/hit-triage-app';
import {HitDesignApp} from './app/hit-design-app';

export const _package = new DG.Package();

//tags: app
//name: Hit Triage
export async function hitTriageApp() {
  new HitTriageApp();
}

//tags: app
//name: Hit Design
export async function hitDesignApp() {
  new HitDesignApp();
}
//name: Demo Molecules 100
//tags: HitTriageDataSource
//output: dataframe result
export async function demoFileIngest(): Promise<DG.DataFrame> {
  const df = grok.data.demo.molecules(100);
  df.name = '100 Molecules';
  return df;
}

//name: Demo Molecules 5000
//tags: HitTriageDataSource
//output: dataframe result
export async function demoFileIngest1(): Promise<DG.DataFrame> {
  const df = grok.data.demo.molecules(5000);
  df.name = '5000 Molecules';
  return df;
}

//name: Demo Molecules variable
//input: int numberOfMolecules [Molecules count]
//tags: HitTriageDataSource
//output: dataframe result
export async function demoFileIngest2(numberOfMolecules: number): Promise<DG.DataFrame> {
  const df = grok.data.demo.molecules(numberOfMolecules);
  df.name = 'Variable Molecules number';
  return df;
}

//name: Demo File Submit
//tags: HitTriageSubmitFunction
//input: dataframe df [Dataframe]
//input: string molecules [Molecules column name]
export async function demoFileSubmit(df: DG.DataFrame, molecules: string): Promise<void> {
  grok.shell.info(df.rowCount);
  grok.shell.info(molecules);
}

//name: Sample File Submit
//tags: HitTriageSubmitFunction
//input: dataframe df [dataframe]
//input: string molecules [molecules column name]
export async function demoFileSubmit1(df: DG.DataFrame, molecules: string): Promise<void> {
  grok.shell.info('Another function');
}
