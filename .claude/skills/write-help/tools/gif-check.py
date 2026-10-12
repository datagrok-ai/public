"""Checks a help GIF recorded with rec-lib.mjs: click marks on the pointer tip, no jump between the first
frames, and the hold of the last frame. Exit code 1 when a check fails.

Usage: python gif-check.py <file.gif> [...]    Requires Pillow and numpy.

Click marks are checked in 2x recordings, rec-lib's default: in a GIF scaled down to half, the marks are too small
to be told from icons, and none are found. It does not judge sharpness: an upscaled recording is caught when it is
made (rec-lib warns when the frames are not twice the viewport).
"""
import sys
import numpy as np
from PIL import Image, ImageSequence

TIP_TOLERANCE = 6     # px between the mark's centre and the pointer tip
JUMP = 100            # px: a change between the first two frames larger than this is a jump
HOLD = (1.0, 6.0)     # s: expected hold of the last frame
STILL = 400           # changed px below which two frames count as the same picture (a blinking caret, a fading balloon)


def mark_mask(a):
    """Click marks: rec-lib's yellow (left button) or green (right button) disc at 75% opacity."""
    r, g, b = a[..., 0].astype(np.int32), a[..., 1].astype(np.int32), a[..., 2].astype(np.int32)
    yellow = (r > 170) & (g * 100 > r * 78) & (g * 100 < r * 93) & (b * 100 < r * 36)
    green = (g > 175) & (r * 100 < g * 55) & (b * 100 < g * 72) & (b * 100 > g * 45)
    return yellow | green


def circle(m, cx, cy):
    """The mark's centre and radius, from its pixels near (cx, cy). The pointer and dark text under the mark hide
    parts of it, but never its top and left edges: the arrow hangs down and to the right of its tip."""
    ys, xs = np.nonzero(m)
    near = (np.abs(xs - cx) < 45) & (np.abs(ys - cy) < 45)
    xs, ys = xs[near], ys[near]
    if len(xs) < 50:
        return cx, cy, 0, False
    left, right = np.percentile(xs, [0.5, 99.5])
    top, bottom = np.percentile(ys, [0.5, 99.5])
    r = max(right - left, bottom - top) / 2
    x0, y0 = left + r, top + r
    d2 = (xs - x0) ** 2 + (ys - y0) ** 2
    round_ = (d2 <= r * r).sum() >= 0.4 * np.pi * r * r and (d2 > (r + 2) ** 2).sum() < 0.1 * len(xs)
    return x0, y0, r, round_


def tip_offset(a, m, cx, cy):
    """Offset of the pointer tip from the mark's centre. Inside the mark only the arrow's fill stays white."""
    x0, y0, r, round_ = circle(m, cx, cy)
    # the mark of a 2x recording (rec-lib's default): a radius of about 24 px
    if not round_ or not 16 <= r <= 34:
        return 'not a mark'
    top, left = int(max(0, y0 - r)), int(max(0, x0 - r))
    win = a[top:int(y0 + r + 1), left:int(x0 + r + 1)]
    wy, wx = np.mgrid[top:top + win.shape[0], left:left + win.shape[1]]
    white = (win.min(axis=2) > 225) & ((wx - x0) ** 2 + (wy - y0) ** 2 <= (r - 3) ** 2)
    if white.sum() < 6:
        return None
    k = np.argmin(np.where(white, wx + wy, 1e9))
    # the fill starts about 2 px right of and 3 px below the tip
    return float(wx.flat[k] - 2 - x0), float(wy.flat[k] - 3 - y0)


def check(path):
    problems, notes = [], []
    im = Image.open(path)
    notes.append(f'{im.size[0]}x{im.size[1]}, {im.n_frames} frames')
    first, prev_m, prev_full, clicks, t, hold, last = [], None, None, 0, 0.0, 0.0, None
    for i, frame in enumerate(ImageSequence.Iterator(im)):
        full = np.asarray(frame.convert('RGB'))
        # scans run at half resolution; a click frame is measured at full resolution
        a = full[::2, ::2].astype(np.int16)
        d = frame.info.get('duration', 0) / 1000
        if i < 2:
            first.append(a)
        m = mark_mask(a)
        new = m if prev_m is None else m & ~prev_m
        if new.sum() > 75:
            ys, xs = np.nonzero(new)
            cx, cy = int(np.median(xs)), int(np.median(ys))
            was = prev_m is not None and prev_m[max(0, cy - 20):cy + 20, max(0, cx - 20):cx + 20].sum() > 20
            if not was:
                # only what appeared in this frame: yellow icons and highlights that were already there drop out
                appeared = mark_mask(full) if prev_full is None else mark_mask(full) & ~mark_mask(prev_full)
                off = tip_offset(full, appeared, 2 * cx, 2 * cy)
                # a click where the previous mark is still fading: only part of the mark is new
                if (off is None or off == 'not a mark') and prev_full is not None:
                    off = tip_offset(full, mark_mask(full), 2 * cx, 2 * cy)
                where = f'click at {t:.2f}s (frame {i})'
                if off == 'not a mark':
                    pass
                elif off is None:
                    problems.append(f'{where}: the pointer tip is not on the mark')
                elif max(abs(off[0]), abs(off[1])) > TIP_TOLERANCE:
                    problems.append(f'{where}: pointer tip is {off[0]:+.0f},{off[1]:+.0f} px from the mark centre')
                clicks += off != 'not a mark'
        prev_m, prev_full = m, full
        changed = last is None or (np.abs(a - last).sum(axis=2) > 30).sum() > STILL / 4
        hold = d if changed else hold + d
        t += d
        last = a
    notes.append(f'{t:.1f} s, {clicks} clicks')
    if len(first) == 2:
        diff = np.abs(first[0] - first[1]).sum(axis=2) > 30
        if diff.any():
            ys, xs = np.nonzero(diff)
            w, h = 2 * (xs.max() - xs.min()), 2 * (ys.max() - ys.min())
            if w > JUMP or h > JUMP:
                problems.append(f'the first two frames differ over {w}x{h} px: a jump at the start')
    notes.append(f'last frame held {hold:.1f} s')
    if not HOLD[0] <= hold <= HOLD[1]:
        problems.append(f'last frame held {hold:.1f} s, expected {HOLD[0]}-{HOLD[1]} s')
    return problems, notes


def main():
    failed = False
    for path in sys.argv[1:]:
        problems, notes = check(path)
        print(f'{path}: {"; ".join(notes)}', flush=True)
        for p in problems:
            print(f'  FAIL {p}', flush=True)
        failed |= bool(problems)
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
