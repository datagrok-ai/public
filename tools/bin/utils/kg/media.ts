/// Media files: what a page shows rather than cites. The one place that says which extensions are
/// media, so `check` (a missing image), the home extractor (no ownership claim, no source-file) and
/// the media extractor read the same rule.
import * as path from 'path';

export type MediaFormat = 'png' | 'jpg' | 'gif' | 'svg' | 'webp' | 'mp4' | 'webm' | 'pdf' | 'youtube' | 'other';

const FORMATS: Record<string, MediaFormat> = {
  '.png': 'png', '.jpg': 'jpg', '.jpeg': 'jpg', '.gif': 'gif', '.svg': 'svg', '.webp': 'webp', '.mp4': 'mp4', '.webm': 'webm',
  '.pdf': 'pdf', '.ico': 'other', '.mov': 'other', '.avif': 'other', '.bmp': 'other',
};

export const MEDIA_EXT = Object.keys(FORMATS);

export function isMedia(file: string): boolean {
  return path.posix.extname(file).toLowerCase() in FORMATS;
}

export function formatOf(file: string): MediaFormat {
  return FORMATS[path.posix.extname(file).toLowerCase()] ?? 'other';
}
