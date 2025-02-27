/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as commonSpec from 'common-spectrum';
import { buildMassSpecViewer } from './SpectrumCalculations/mass-spectrum';

import '../css/spectroscopy.css';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: readJDXFile
//input: file file
export async function readJDXFile(fileData: DG.FileInfo) {
  const fileString = await fileData.readAsString() ?? '';
  const spectrumTypeRergex = new RegExp("##DATA\\s*TYPE\\s*=\\s*([A-Za-z\\/]*) SPECTRUM");
  const spectrumType = ((fileString).match(spectrumTypeRergex) ?? [])[1];

  if (spectrumType) {
    if (spectrumType === 'MASS') {
      const molFormRegex = new RegExp("##MOLFORM\s*=\s*([^\n]*)");
      const molForm = (fileString.match(molFormRegex) ?? [])[1];
      buildMassSpecViewer(fileData.name, fileString, molForm);
    }
    if (spectrumType === 'UV/VIS')
      buildUVSpecViewer(fileData.name, fileString);
    if (spectrumType === 'INFRARED')
      buildIRSpecViewer(fileData.name, fileString);
    if (spectrumType === 'NMR')
    {
      //TODO call from nmr
    }
  }
}

function buildIRSpecViewer(name: string, data: string) {
  const irSpec = commonSpec.fromJcamp(data);
  const xyPoints = irSpec.getXY();
  
  //@ts-ignore
  let xColumn = DG.Column.fromFloat32Array('x', (xyPoints.x ?? []) as Float32Array);
  //@ts-ignore
  let yColumn = DG.Column.fromFloat32Array('y', (xyPoints.y ?? []) as Float32Array);
  let df = DG.DataFrame.fromColumns([xColumn, yColumn]);
  let lineChart = DG.Viewer.lineChart(df);
  let view = DG.View.fromRoot(lineChart.root);
  view.name = name;
  lineChart.root.style.width = '100%';
  lineChart.root.style.height = '100%'; 

  irSpec.getNormalizedSpectrum();
  irSpec.getNormalizedSpectra();
  grok.shell.addView(view);
}

function buildUVSpecViewer(name: string, data: string) {
  const irSpec = commonSpec.fromJcamp(data);
  const xyPoints = irSpec.getXY();
  
  //@ts-ignore
  let xColumn = DG.Column.fromFloat32Array('x', (xyPoints.x ?? []) as Float32Array);
  //@ts-ignore
  let yColumn = DG.Column.fromFloat32Array('y', (xyPoints.y ?? []) as Float32Array);
  let df = DG.DataFrame.fromColumns([xColumn, yColumn]);
  let lineChart = DG.Viewer.lineChart(df);
  let view = DG.View.fromRoot(lineChart.root);
  view.name = name;
  lineChart.root.style.width = '100%';
  lineChart.root.style.height = '100%'; 

  irSpec.getNormalizedSpectrum();
  irSpec.getNormalizedSpectra();
  grok.shell.addView(view);
}