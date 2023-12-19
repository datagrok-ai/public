import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {v4 as uuidv4} from 'uuid';

import {FileInfo as msFileInfo} from 'molstar/lib/mol-util/file-info';
import {Color as msColor} from 'molstar/lib/mol-util/color';
import {StateObjectSelector} from 'molstar/lib/mol-state';
import {PluginStateObject} from 'molstar/lib/mol-plugin-state/objects';
import {StateTransforms} from 'molstar/lib/mol-plugin-state/transforms';
import {createVolumeRepresentationParams} from 'molstar/lib/mol-plugin-state/helpers/volume-representation-params';
import {PluginUIContext} from 'molstar/lib/mol-plugin-ui/context';
import {VolumeIsovalueInfo} from 'molstar/lib/apps/viewer';

import {BiostructureData} from '@datagrok-libraries/bio/src/pdb/types';

import {molecule3dFileExtensions} from './consts';

enum ParsedParts {
  format = 'format', // No specific treatment
  trajectory = 'trajectory',
  structure = 'structure',
  topology = 'topology',
  shape = 'shape',
  volume = 'volume',
  volumes = 'volumes',
  // density = 'density',
}

export async function parseAndVisualsData(plugin: PluginUIContext, dataEff: BiostructureData): Promise<void> {
  if (!dataEff) return;

  const {binary} = molecule3dFileExtensions[dataEff.ext];
  const data: PluginStateObject.Data.Binary | PluginStateObject.Data.String = binary ?
    new PluginStateObject.Data.Binary(dataEff.data as Uint8Array) :
    new PluginStateObject.Data.String(dataEff.data as string);

  //await plugin.dataTransaction(async () => {
  const dataProvider = plugin.dataFormats.auto({ext: dataEff.ext} as msFileInfo, data)!;
  if (!dataProvider)
    throw new Error(`Can not find data provider for file ext '${dataEff.ext}'.`);

  const entryId = uuidv4();
  let dataVal: BiostructureData;
  if (dataEff.ext === 'pdb' && dataEff.binary) {
    dataVal = {
      binary: false, ext: dataEff.ext, options: dataEff.options,
      data: (new TextDecoder()).decode(dataEff.data as Uint8Array)
    };
  } else
    dataVal = dataEff;

  const _data = await plugin.builders.data.rawData(
    {data: dataVal.data, label: dataVal.options?.dataLabel});
  const parsed = await dataProvider.parse(plugin, _data, {entryId: entryId});

  if (parsed.hasOwnProperty(ParsedParts.trajectory) || parsed.hasOwnProperty(ParsedParts.structure) ||
    parsed.hasOwnProperty(ParsedParts.topology) || parsed.hasOwnProperty(ParsedParts.shape) ||
    parsed.hasOwnProperty(ParsedParts.format)
  ) {
    if (dataProvider.visuals) {
      const visualsPromise = dataProvider.visuals(plugin, parsed);
      if (visualsPromise) await visualsPromise;
    }
  }

  if (parsed.hasOwnProperty(ParsedParts.volume) || parsed.hasOwnProperty(ParsedParts.volumes)) {
    const firstVolume = (parsed.volume || parsed.volumes[0]) as StateObjectSelector<PluginStateObject.Volume.Data>;
    const repr = plugin.build();

    // Defaults
    const isovalues: VolumeIsovalueInfo[] = [{
      type: 'relative',
      value: 1,
      color: msColor(0x3377aa),
    }];

    for (const iso of isovalues) {
      repr
        .to(parsed.volumes?.[iso.volumeIndex ?? 0] ?? parsed.volume)
        .apply(
          StateTransforms.Representation.VolumeRepresentation3D,
          createVolumeRepresentationParams(plugin, firstVolume.data!, {
            type: 'isosurface',
            typeParams: {
              alpha: iso.alpha ?? 1,
              isoValue: iso.type === 'absolute' ? {kind: 'absolute', absoluteValue: iso.value} : {
                kind: 'relative',
                relativeValue: iso.value,
              },
            },
            color: 'uniform',
            colorParams: {value: iso.color},
          }));
    }

    await repr.commit();
  }

  const unhandledParsedPartList = Object.getOwnPropertyNames(parsed)
    .filter((propName) => !(Object.values(ParsedParts) as string[]).includes(propName));
  if (unhandledParsedPartList.length > 0)
    throw new Error(`Unhandled parsed parts '${JSON.stringify(unhandledParsedPartList)}'.`);
}
