"""Cuts frames from the start or the end of a GIF without re-encoding the frames it keeps.

Usage: python gif-trim.py <in.gif> <out.gif> [--start N] [--end-at SECONDS] [--hold SECONDS]
  --start N         drop the first N frames; frame N becomes the first one
  --end-at SECONDS  drop the frames that start at or after this time
  --hold SECONDS    delay of the last kept frame (default: unchanged, 2.0 with --end-at)

Find the frames with gif-check.py or a contact sheet first. Requires Pillow and numpy.
"""
import argparse
import io
import numpy as np
from PIL import Image, ImageSequence


def blocks(b):
    """Splits a GIF into its header (with the global colour table) and frames: [GCE bytes, image bytes, delay s]."""
    p = 13
    if b[10] & 0x80:
        p += 3 * (2 << (b[10] & 7))
    header, frames, extensions, gce = b[:p], [], [], b''

    def skip_sub(q):
        while b[q]:
            q += b[q] + 1
        return q + 1

    while b[p] != 0x3B:
        if b[p] == 0x21 and b[p + 1] == 0xF9:
            gce = b[p:p + 8]
            p += 8
        elif b[p] == 0x21:
            e = skip_sub(p + 2)
            extensions.append(b[p:e])
            p = e
        elif b[p] == 0x2C:
            q = p + 10
            if b[p + 9] & 0x80:
                q += 3 * (2 << (b[p + 9] & 7))
            e = skip_sub(q + 1)
            delay = (gce[4] | gce[5] << 8) / 100 if gce else 0
            frames.append([gce, b[p:e], delay])
            gce = b''
            p = e
        else:
            raise ValueError(f'unexpected block {b[p]:#x} at {p}')
    return header, extensions, frames


def full_frame(src, k):
    """Frame k as it is shown (frames after the first usually redraw only what changed), as one full frame."""
    im = Image.open(src)
    im.seek(k)
    rgb = im.convert('RGB')
    buf = io.BytesIO()
    rgb.quantize(256, dither=Image.Dither.NONE).save(buf, 'GIF', interlace=False)
    t = buf.getvalue()
    bits = t[10] & 7
    table = t[13:13 + 3 * (2 << bits)]
    p = 13 + len(table)
    while t[p] != 0x2C:
        p += 2
        while t[p]:
            p += t[p] + 1
        p += 1
    e = p + 11
    while t[e]:
        e += t[e] + 1
    w, h = rgb.size
    image = bytes([0x2C, 0, 0, 0, 0, w & 255, w >> 8, h & 255, h >> 8, 0x80 | bits]) + table + t[p + 10:e + 1]
    return image, rgb


def gce_with(delay, disposal=1):
    d = int(round(delay * 100))
    return bytes([0x21, 0xF9, 4, disposal << 2, d & 255, d >> 8, 0, 0])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('src')
    ap.add_argument('dst')
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--end-at', type=float)
    ap.add_argument('--hold', type=float)
    a = ap.parse_args()

    b = open(a.src, 'rb').read()
    header, extensions, frames = blocks(b)
    total = len(frames)
    t, kept = 0.0, []
    for i, f in enumerate(frames):
        if a.end_at is not None and t >= a.end_at:
            break
        if i >= a.start:
            kept.append(f)
        t += f[2]
    if not kept:
        raise SystemExit('nothing left to keep')

    if a.start:
        image, expected = full_frame(a.src, a.start)
        kept[0] = [gce_with(kept[0][2]), image, kept[0][2]]
    hold = a.hold if a.hold is not None else (2.0 if a.end_at is not None else None)
    if hold is not None:
        g = bytearray(kept[-1][0] or gce_with(0))
        d = int(round(hold * 100))
        g[4], g[5] = d & 255, d >> 8
        kept[-1][0] = bytes(g)

    out = bytearray(header)
    for e in extensions:
        out += e
    for g, image, _ in kept:
        out += g + image
    out.append(0x3B)
    open(a.dst, 'wb').write(out)

    note = ''
    if a.start:
        got = np.asarray(Image.open(a.dst).convert('RGB')).astype(int)
        diff = np.abs(got - np.asarray(expected).astype(int)).max()
        note = ', first frame identical' if diff == 0 else f', first frame differs by up to {diff} (more than 256 colours)'
    dur = sum(f[2] for f in kept) if hold is None else sum(f[2] for f in kept[:-1]) + hold
    print(f'{a.dst}: {len(kept)} of {total} frames, {dur:.1f} s{note}')


if __name__ == '__main__':
    main()
