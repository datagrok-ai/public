import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {parseFcs, fcsToDataFrame} from './fcs/fcs-parser';
import {FcsParsed} from './fcs/fcs-types';


export const _package = new DG.Package();

export * from './package.g';


function findChannel(parsed: FcsParsed, candidates: string[]): string | null {
  for (const candidate of candidates) {
    const upper = candidate.toUpperCase();
    for (const p of parsed.parameters) {
      if (p.shortName.toUpperCase() === upper)
        return p.shortName;
    }
  }
  return null;
}

function buildMetadataPanel(parsed: FcsParsed): HTMLElement {
  const kw = parsed.keywords;
  const items: HTMLElement[] = [];

  items.push(ui.h2('FCS File Info'));
  const infoRows: string[][] = [];
  const addRow = (label: string, value: string | undefined) => {
    if (value)
      infoRows.push([label, value]);
  };
  addRow('Version', parsed.header.version);
  addRow('Instrument', kw.get('$CYT'));
  addRow('Date', kw.get('$DATE'));
  addRow('Events', String(parsed.eventCount));
  addRow('Parameters', String(parsed.paramCount));
  addRow('Data type', parsed.dataType === 'F' ? 'Float32' : parsed.dataType === 'D' ? 'Float64' : 'Integer');
  addRow('Byte order', parsed.littleEndian ? 'Little-endian' : 'Big-endian');
  addRow('Well ID', kw.get('$WELLID'));
  addRow('Plate ID', kw.get('$PLATEID'));
  addRow('Source', kw.get('$SRC'));

  const infoTable = ui.tableFromMap(Object.fromEntries(infoRows));
  items.push(infoTable);

  items.push(ui.h2('Parameters'));
  const paramRows = parsed.parameters.map((p) => ({
    '#': p.index,
    'Name': p.shortName,
    'Marker': p.longName || '',
    'Bits': p.bits,
    'Range': p.range,
  }));
  const paramDf = DG.DataFrame.fromObjects(paramRows)!;
  const grid = DG.Viewer.grid(paramDf);
  grid.root.style.maxHeight = '200px';
  items.push(grid.root);

  if (parsed.spillover)
    items.push(ui.divText(`Spillover matrix: ${parsed.spillover.n} parameters`));

  return ui.divV(items);
}


export class FlowCytometryPackageFunctions {
  @grok.decorators.fileHandler({ext: 'fcs', description: 'Opens FCS (Flow Cytometry Standard) file'})
  static async importFcs(
    @grok.decorators.param({type: 'list'}) bytes: Uint8Array
  ): Promise<DG.DataFrame[]> {
    const buffer = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) as ArrayBuffer;
    const parsed = parseFcs(buffer);
    return [fcsToDataFrame(parsed)];
  }

  @grok.decorators.fileViewer({fileViewer: 'fcs'})
  static async previewFcs(file: DG.FileInfo): Promise<DG.View> {
    const view = DG.View.create();
    view.name = file.name;

    const bytes = await file.readAsBytes();
    const buffer = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) as ArrayBuffer;
    const parsed = parseFcs(buffer);
    const df = fcsToDataFrame(parsed);

    const metaPanel = buildMetadataPanel(parsed);
    metaPanel.style.padding = '8px';
    view.append(metaPanel);

    const xName = findChannel(parsed, ['FSC-A', 'FSC-H', 'FSC.A', 'FSC.H']) ?? parsed.parameters[0]?.shortName;
    const yName = findChannel(parsed, ['SSC-A', 'SSC-H', 'SSC.A', 'SSC.H']) ?? parsed.parameters[1]?.shortName;
    if (xName && yName) {
      const scatter = DG.Viewer.scatterPlot(df, {x: xName, y: yName});
      scatter.root.style.width = '100%';
      scatter.root.style.height = '400px';
      view.append(scatter.root);
    }

    return view;
  }
}
