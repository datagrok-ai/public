import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../css/media-viewer.css';

export * from './package.g';
export const _package = new DG.Package();

const AUDIO_MIME: Record<string, string> = {
  mp3: 'audio/mpeg',
  wav: 'audio/wav',
  flac: 'audio/flac',
  ogg: 'audio/ogg',
  oga: 'audio/ogg',
  opus: 'audio/ogg',
  m4a: 'audio/mp4',
  aac: 'audio/aac',
  weba: 'audio/webm',
};

const VIDEO_MIME: Record<string, string> = {
  mp4: 'video/mp4',
  m4v: 'video/mp4',
  webm: 'video/webm',
  ogv: 'video/ogg',
  mov: 'video/quicktime',
};

function mediaView(file: DG.FileInfo, kind: 'audio' | 'video',
  mimes: Record<string, string>): DG.View {
  const view = DG.View.create();
  view.name = file.name;

  const host = ui.div([], `media-viewer-host media-viewer-${kind}`);
  view.append(host);

  let objectUrl: string | null = null;

  file.readAsBytes().then((bytes) => {
    const ext = (file.extension ?? '').toLowerCase();
    const type = mimes[ext] ?? '';
    const blob = new Blob([bytes as BlobPart], {type});
    objectUrl = URL.createObjectURL(blob);

    const el = document.createElement(kind) as HTMLMediaElement;
    el.controls = true;
    el.preload = 'metadata';
    el.src = objectUrl;
    el.className = `media-viewer-${kind}-element`;
    host.appendChild(el);
  }).catch((err) => {
    const msg = err instanceof Error ? err.message : String(err);
    grok.shell.error(`Media preview error: ${msg}`);
    host.appendChild(ui.divText('Failed to load media file.', 'media-viewer-error'));
  });

  const sub = grok.events.onViewRemoved.subscribe((removed: DG.ViewBase) => {
    if (removed !== view) return;
    if (objectUrl !== null) {
      URL.revokeObjectURL(objectUrl);
      objectUrl = null;
    }
    sub.unsubscribe();
  });

  return view;
}

//name: audioViewer
//description: Previews audio files using the native browser player
//meta.role: fileViewer
//meta.fileViewer: mp3, wav, flac, ogg, oga, opus, m4a, aac, weba
//input: file file
//output: view v
export function audioViewer(file: DG.FileInfo): DG.View {
  return mediaView(file, 'audio', AUDIO_MIME);
}

//name: videoViewer
//description: Previews video files using the native browser player
//meta.role: fileViewer
//meta.fileViewer: mp4, m4v, webm, ogv, mov
//input: file file
//output: view v
export function videoViewer(file: DG.FileInfo): DG.View {
  return mediaView(file, 'video', VIDEO_MIME);
}
