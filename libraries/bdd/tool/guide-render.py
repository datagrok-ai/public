"""Renders a guide manifest (`steps.json` + the PNGs a `BDD_GUIDE` run wrote) into guide.mp4,
step-NN.png (the annotated picture of every step), steps.md, audit.png + audit.json, and with --gif
also guide.gif and guide-thumb.png. Pillow does the drawing, ffmpeg the encoding (FFMPEG env,
PATH, the imageio-ffmpeg package, or Playwright's own build, which can only write WebM). The video
and the stills are 1600×1000 with the caption strip: the page is filmed at 1080p and downsampled.

The pointer goes where the page says it went, and never skips: every frame's pointer continues
from the frame before, so an action starts where the one before it ended — a skip is an error,
not a frame. A press is marked where it landed, under the pointer's tip: a yellow ring for the
left button, a green one for the right. audit.png shows every press as the video shows it, and
audit.json how far the mark's centre and the pointer's tip are from the press.

    py tool/guide-render.py <scenario dir> [--gif] [--fps 30] [--hold 1.4] [--travel 0.7]
        [--zoom 1.8] [--zoom-time 0.8] [--ffmpeg <path>] [--quiet]
    py tool/guide-render.py --selftest <dir>      renders a synthetic guide there, and checks its audit
"""
import glob
import json
import math
import os
import shutil
import subprocess
import sys

try:
    from PIL import Image, ImageChops, ImageDraw, ImageFont
except ImportError:
    print('guide-render: Pillow is missing — py -m pip install pillow', file=sys.stderr)
    sys.exit(3)

ACCENT = (255, 90, 31)
CHECK = (46, 160, 67)
BAR = (20, 24, 32, 205)
LEFT_PRESS = (255, 200, 0)    # a left-button press: yellow
RIGHT_PRESS = (34, 197, 94)   # a right-button press: green
SMALL_TARGET_PX = 28   # a target no wider and no taller than this (an icon, a checkbox) is zoomed into; larger ones are lit and clicked in place
VIDEO_W, VIDEO_H = 1600, 1000   # the video and the stills, caption strip included; the page is filmed larger and downsampled
PULSE_S = 0.5   # how long a press ripples
MAX_STEP_PX = 70   # the farthest the pointer moves between two frames, in pixels of a 1080p page: a longer way takes more frames
AUDIT_TOLERANCE_PX = 2   # how far, in video pixels, a press mark's centre or the pointer's tip may sit from the press


def find_ffmpeg(explicit):
    if explicit:
        return explicit
    if os.environ.get('FFMPEG'):
        return os.environ['FFMPEG']
    on_path = shutil.which('ffmpeg')
    if on_path:
        return on_path
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        pass
    home = os.environ.get('PLAYWRIGHT_BROWSERS_PATH') or os.path.join(
        os.environ.get('LOCALAPPDATA', os.path.expanduser('~/.cache')), 'ms-playwright')
    found = sorted(glob.glob(os.path.join(home, 'ffmpeg-*', 'ffmpeg*')))
    return found[-1] if found else None


def has_encoder(ffmpeg, name):
    out = subprocess.run([ffmpeg, '-hide_banner', '-encoders'], capture_output=True, text=True).stdout
    return f' {name} ' in out


def load_font(size):
    for path in ['C:/Windows/Fonts/segoeui.ttf', 'C:/Windows/Fonts/arial.ttf',
                 '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', '/System/Library/Fonts/Helvetica.ttc']:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    try:
        return ImageFont.load_default(size)
    except TypeError:
        return ImageFont.load_default()


def ease(t):
    t = max(0.0, min(1.0, t))
    return t * t * (3 - 2 * t)


def cursor_sprite(scale):
    """The arrow, its tip at the sprite's (1, 1)."""
    pts = [(0, 0), (0, 16), (4, 12.5), (7, 19.5), (10, 18), (7, 11.5), (12.5, 11.5)]
    pts = [(x * scale, y * scale) for x, y in pts]
    w = int(14 * scale) + 4
    h = int(21 * scale) + 4
    sprite = Image.new('RGBA', (w, h), (0, 0, 0, 0))
    d = ImageDraw.Draw(sprite)
    d.polygon([(x + 2.5, y + 2.5) for x, y in pts], fill=(0, 0, 0, 90))
    d.polygon([(x + 1, y + 1) for x, y in pts], fill=(255, 255, 255, 255), outline=(0, 0, 0, 255))
    return sprite


def overlay(frame, sprite, x, y):
    """Composites the sprite with its top left at (x, y), cut to the frame."""
    x0, y0 = max(0, x), max(0, y)
    x1, y1 = min(frame.width, x + sprite.width), min(frame.height, y + sprite.height)
    if x1 > x0 and y1 > y0:
        frame.alpha_composite(sprite.crop((x0 - x, y0 - y, x1 - x, y1 - y)), (x0, y0))


def inside(box, x, y, pad=2):
    return box is not None and box['x'] - pad <= x <= box['x'] + box['width'] + pad and \
        box['y'] - pad <= y <= box['y'] + box['height'] + pad


def is_small(box):
    return box is not None and box['width'] <= SMALL_TARGET_PX and box['height'] <= SMALL_TARGET_PX


class Renderer:
    def __init__(self, folder, fps, hold, travel, zoom, zoom_time=0.8):
        self.folder = folder
        self.fps = fps
        self.hold = hold
        self.travel = travel
        self.zoom = zoom
        self.zoom_time = zoom_time
        with open(os.path.join(folder, 'steps.json'), encoding='utf-8') as f:
            self.manifest = json.load(f)
        self.steps = [s for s in self.manifest['steps'] if s['kind'] != 'setup']
        first = self._image(self.steps[0]['before'] if self.steps else self.manifest['steps'][0]['before'])
        self.w, self.h = first.width & ~1, first.height & ~1
        self.s = self.h / 1080
        self.big = load_font(max(20, self.h // 18))
        self.cursor = cursor_sprite(max(1.0, self.h / 700))
        # the caption sits above the page: nothing of the page is covered, and a player's timeline
        # lands on the bottom edge. The page keeps the pixels it was filmed at (every box, pointer
        # and zoom is in them); the strip is as tall as makes page plus strip the video's aspect
        # (120 px over a 1080p page), and the composed frame is downsampled to the video size —
        # never upsampled, so a small capture stays as it is
        self.bar_h = max(int(self.h * 0.085), round(self.w * VIDEO_H / VIDEO_W) - self.h) & ~1
        self.out_h = self.h + self.bar_h
        scale = min(1.0, VIDEO_W / self.w, VIDEO_H / self.out_h)
        self.video = (int(self.w * scale) & ~1, int(self.out_h * scale) & ~1)
        self.font = load_font(max(14, self.bar_h // 3))
        self.small = load_font(max(12, self.bar_h // 4))
        self.cur = (self.w / 2, self.h / 2)
        self.drawn = None
        self.frame_no = 0
        self.presses = []

    def _image(self, name):
        img = Image.open(os.path.join(self.folder, name)).convert('RGB')
        if hasattr(self, 'w') and (img.width != self.w or img.height != self.h):
            img = img.crop((0, 0, self.w, self.h))
        return img

    def spotlight(self, img, box, color=ACCENT):
        """The frame with everything but the box dimmed and the box outlined."""
        pad = 6
        x0, y0 = max(0, box['x'] - pad), max(0, box['y'] - pad)
        x1, y1 = min(self.w, box['x'] + box['width'] + pad), min(self.h, box['y'] + box['height'] + pad)
        layer = Image.new('RGBA', (self.w, self.h), (0, 0, 0, 105))
        d = ImageDraw.Draw(layer)
        d.rounded_rectangle((x0, y0, x1, y1), radius=6, fill=(0, 0, 0, 0))
        out = img.convert('RGBA')
        out.alpha_composite(layer)
        d = ImageDraw.Draw(out)
        d.rounded_rectangle((x0, y0, x1, y1), radius=6, outline=color + (255,), width=3)
        return out.convert('RGB')

    def caption_layer(self, step, index, total):
        layer = Image.new('RGBA', (self.w, self.bar_h), BAR)
        d = ImageDraw.Draw(layer)
        pad = self.bar_h // 3
        counter = f'{index} / {total}'
        d.text((pad, self.bar_h / 2), counter, font=self.small, fill=(170, 178, 190, 255), anchor='lm')
        x = pad + d.textlength(counter, font=self.small) + pad
        text = step['caption']
        if step['kind'] == 'check':
            s = self.font.size
            mid = self.bar_h / 2
            d.line([(x, mid), (x + s * 0.32, mid + s * 0.32), (x + s * 0.85, mid - s * 0.4)], fill=CHECK + (255,), width=3)
            x += s * 1.2
        labels = [key.replace('ControlOrMeta', 'Ctrl') for key in step.get('keys', [])[:4]]
        keys_w = sum(d.textlength(label, font=self.small) + 22 for label in labels)
        font, lines = self.fit_caption(d, text, self.w - x - pad - keys_w)
        gap = font.size * 0.62 if len(lines) > 1 else 0
        for i, line in enumerate(lines):
            y = self.bar_h / 2 + (2 * i - len(lines) + 1) * gap
            d.text((x, y), line, font=font, fill=(255, 255, 255, 255), anchor='lm')
        x += max(d.textlength(line, font=font) for line in lines) + pad
        for label in labels:
            tw = d.textlength(label, font=self.small)
            d.rounded_rectangle((x, self.bar_h * 0.22, x + tw + 14, self.bar_h * 0.78), radius=5,
                                fill=(60, 66, 78, 255), outline=(120, 128, 140, 255))
            d.text((x + 7, self.bar_h / 2), label, font=self.small, fill=(255, 255, 255, 255), anchor='lm')
            x += tw + 22
        return layer

    def fit_caption(self, d, text, width):
        """The caption at the bar's font size; if it is wider than the bar leaves it, the largest
        size down to seven tenths of that which fits on one line, else two lines at that size, the
        second cut short if it still overflows."""
        floor = max(12, round(self.font.size * 0.7))
        for size in range(self.font.size, floor - 1, -1):
            font = self.font if size == self.font.size else load_font(size)
            if d.textlength(text, font=font) <= width:
                return font, [text]
        words, first = text.split(' '), ''
        while words and d.textlength((first + ' ' + words[0]).lstrip(), font=font) <= width:
            first = (first + ' ' + words.pop(0)).lstrip()
        rest = ' '.join(words)
        if d.textlength(rest, font=font) > width:
            while rest and d.textlength(rest + '…', font=font) > width:
                rest = rest[:-1]
            rest = rest.rstrip() + '…'
        return font, [first, rest] if first else [rest]

    def press_mark(self, frame, at, t, button, ring=True):
        """A press `t` of the way through its ripple: a dot where the button went down and a ring
        spreading from it — yellow for the left button, green for the right. Both stay solid for
        most of it, so the page under them cannot tint them (a faded yellow over a blue button
        reads green), and fade at the end. Without the ring it is the button held through a drag.
        Drawn three times larger and shrunk, so the circles are smooth; a thin dark rim keeps both
        visible on a white page."""
        color = RIGHT_PRESS if button == 'right' else LEFT_PRESS
        s, k = self.s, 3
        half = int(56 * s) + 4
        patch = Image.new('RGBA', (2 * half * k, 2 * half * k), (0, 0, 0, 0))
        d = ImageDraw.Draw(patch)
        c = half * k
        e = ease(t)
        dot = 11 * s * k
        solid = int(235 * (1 - ease((t - 0.5) / 0.5))) if ring else 235
        d.ellipse((c - dot - k, c - dot - k, c + dot + k, c + dot + k), fill=(0, 0, 0, solid // 3))
        d.ellipse((c - dot, c - dot, c + dot, c + dot), fill=color + (solid,))
        if ring:
            r = (12 + 36 * e) * s * k
            w = max(k, round((5 - 2.5 * e) * s * k))
            a = int(255 * (1 - ease((t - 0.6) / 0.4)))
            d.ellipse((c - r - k, c - r - k, c + r + k, c + r + k), outline=(0, 0, 0, a // 3), width=w + 2 * k)
            d.ellipse((c - r, c - r, c + r, c + r), outline=color + (a,), width=w)
        patch = patch.resize((2 * half, 2 * half), Image.LANCZOS)
        overlay(frame, patch, round(at[0]) - half, round(at[1]) - half)

    def compose(self, base, cursor=None, press=None, held=None, zoom=1.0, anchor=None, caption=None):
        """A picture of the video: the page, a press and the pointer's tip on the same point, the
        zoom, the caption strip above."""
        frame = base.convert('RGBA')
        if press is not None:
            self.press_mark(frame, *press)
        if held is not None and cursor is not None:
            self.press_mark(frame, cursor, 0.0, held, ring=False)
        if cursor is not None:
            overlay(frame, self.cursor, round(cursor[0]) - 1, round(cursor[1]) - 1)
        frame = frame.convert('RGB')
        if zoom > 1.001 and anchor:
            x0, y0, x1, y1 = self.zoom_box(zoom, anchor)
            frame = frame.crop((x0, y0, x1, y1)).resize((self.w, self.h), Image.BILINEAR)
        canvas = Image.new('RGB', (self.w, self.out_h), BAR[:3])
        canvas.paste(frame, (0, self.bar_h))
        if caption is not None:
            canvas.paste(caption, (0, 0), caption)
        return canvas if canvas.size == self.video else canvas.resize(self.video, Image.LANCZOS)

    def zoom_box(self, zoom, anchor):
        cw, ch = self.w / zoom, self.h / zoom
        x0 = min(max(0, anchor[0] - cw / 2), self.w - cw)
        y0 = min(max(0, anchor[1] - ch / 2), self.h - ch)
        return int(x0), int(y0), int(x0 + cw), int(y0 + ch)

    def to_video(self, at, zoom=1.0, anchor=None):
        """Where a point of the page lands in a frame of the video."""
        x, y = at
        if zoom > 1.001 and anchor:
            x0, y0, x1, y1 = self.zoom_box(zoom, anchor)
            x, y = (x - x0) * self.w / (x1 - x0), (y - y0) * self.h / (y1 - y0)
        return x * self.video[0] / self.w, (y + self.bar_h) * self.video[1] / self.out_h

    def frame(self, base, cursor=None, **kw):
        """Every frame of the video comes through here: the pointer it draws continues from the
        pointer of the frame before, or the video is not made."""
        if cursor is not None:
            if self.drawn is not None:
                skip = math.hypot(cursor[0] - self.drawn[0], cursor[1] - self.drawn[1])
                if skip > MAX_STEP_PX * self.s + 0.5:
                    raise RuntimeError(f'the pointer skipped {skip:.0f} px at frame {self.frame_no}: '
                                       f'({self.drawn[0]:.0f}, {self.drawn[1]:.0f}) to ({cursor[0]:.0f}, {cursor[1]:.0f})')
            self.drawn = self.cur = cursor
        self.frame_no += 1
        return self.compose(base, cursor, **kw)

    def hold_frames(self, seconds, base, cursor=None, **kw):
        n = int(seconds * self.fps)
        if n <= 0:
            return
        still = self.frame(base, cursor, **kw)
        yield still
        for _ in range(n - 1):
            self.frame_no += 1
            yield still

    def title_frames(self, base):
        layer = Image.new('RGBA', (self.w, self.h), (10, 14, 20, 150))
        frame = Image.alpha_composite(base.convert('RGBA'), layer).convert('RGB')
        d = ImageDraw.Draw(frame)
        d.text((self.w / 2, self.h * 0.44), self.manifest['scenario'], font=self.big, fill=(255, 255, 255), anchor='mm')
        d.text((self.w / 2, self.h * 0.44 + self.h * 0.09), self.manifest['feature'], font=self.font,
               fill=(200, 206, 214), anchor='mm')
        yield from self.hold_frames(1.8, frame)

    def frames_for(self, distance, least):
        """Enough frames that the eased way never moves the pointer more than MAX_STEP_PX in one."""
        return max(1, least, math.ceil(1.5 * distance / (MAX_STEP_PX * self.s)))

    def travel_frames(self, base, lit, dest, caption, fade=False):
        """The pointer eased from where the last frame left it to dest — a longer way takes longer —
        with the stop lighting up over the second half when `fade`."""
        start = self.cur
        distance = math.hypot(dest[0] - start[0], dest[1] - start[1])
        if distance < 0.5:
            if fade:
                yield from self.light_frames(base, lit, caption)
            return
        reach = 1600 * self.s
        n = self.frames_for(distance, round(self.fps * self.travel * (0.5 + min(distance, reach) / reach)))
        for i in range(n):
            t = ease((i + 1) / n)
            pos = (start[0] + (dest[0] - start[0]) * t, start[1] + (dest[1] - start[1]) * t)
            mix = ease((t - 0.5) / 0.5) if fade else 1.0
            yield self.frame(lit if mix >= 1 else Image.blend(base, lit, mix) if mix > 0 else base, pos, caption=caption)

    def light_frames(self, base, lit, caption):
        """A stop lighting up where the pointer already is — one the step reached with the keyboard."""
        n = int(0.3 * self.fps)
        for i in range(n):
            yield self.frame(Image.blend(base, lit, ease((i + 1) / n)), self.cur, caption=caption)

    def zoom_frames(self, lit, at, z0, z1, caption):
        n = int(self.zoom_time * self.fps)
        for i in range(n):
            z = z0 + (z1 - z0) * ease((i + 1) / n)
            yield self.frame(lit, at, zoom=z, anchor=at, caption=caption)

    def press_frames(self, lit, at, press, target, step, index, caption, zoom, anchor):
        """The ripple of a press, the pointer's tip on its centre; the frame a fifth of the way in is
        kept for the audit, with the same frame bare and with the mark only, to measure both."""
        button = press.get('button', 'left')
        n = int((PULSE_S if (press.get('count') or 1) == 1 else PULSE_S * 0.6) * self.fps)
        peak = max(0, n // 5)
        for i in range(n):
            t = (i + 1) / n
            image = self.frame(lit, at, press=(at, t, button), zoom=zoom, anchor=anchor, caption=caption)
            if i == peak:
                self.presses.append({
                    'frame': self.frame_no - 1, 'step': index, 'caption': step['caption'], 'x': at[0], 'y': at[1],
                    'button': button, 'count': press.get('count') or 1, 'on': press.get('on'), 'target': target,
                    'inside': None if target is None else inside(target, at[0], at[1]),
                    'video': self.to_video(at, zoom, anchor),
                    'images': (image,
                               self.compose(lit, None, press=(at, t, button), zoom=zoom, anchor=anchor, caption=caption),
                               self.compose(lit, None, zoom=zoom, anchor=anchor, caption=caption))})
            yield image

    def drag_frames(self, lit, target, path, button, caption):
        """The pointer along a drag, the button held: each picture taken on the way fades in as the
        pointer reaches it, so what the drag draws grows under the pointer. Returns the last one."""
        keys = [p for p in path if p.get('shot')]
        if not path[-1].get('shot'):
            keys.append(path[-1])
        for p in keys:
            dest = (p['x'], p['y'])
            drawn = None
            if p.get('shot'):
                drawn = self._image(p['shot'])
                drawn = self.spotlight(drawn, target) if target else drawn
            start = self.cur
            n = self.frames_for(math.hypot(dest[0] - start[0], dest[1] - start[1]), int(0.3 * self.fps))
            for i in range(n):
                t = ease((i + 1) / n)
                pos = (start[0] + (dest[0] - start[0]) * t, start[1] + (dest[1] - start[1]) * t)
                mix = max(0.0, (t - 0.6) / 0.4) if drawn is not None else 0.0
                yield self.frame(Image.blend(lit, drawn, mix) if mix > 0 else lit, pos, held=button, caption=caption)
            if drawn is not None:
                lit = drawn
        yield from self.hold_frames(0.5, lit, self.cur, caption=caption)
        return lit

    def step_frames(self, step, index, total):
        fps = self.fps
        after = self._image(step['after'])
        caption = self.caption_layer(step, index, total)
        legs = step.get('legs', [])
        still_path = os.path.join(self.folder, f'step-{index:02d}.png')
        if step['kind'] == 'check':
            target = step.get('target')
            base = self.spotlight(after, target, CHECK) if target else after
            # a check that moved the pointer (a hover for a tooltip) shows it there
            rests = [(p['x'], p['y']) for leg in legs for p in leg.get('pointer', []) if p['type'] != 'up']
            if rests:
                yield from self.travel_frames(base, base, rests[-1], caption)
            self.compose(base, self.cur, caption=caption).save(still_path)
            yield from self.hold_frames(self.hold * 0.9, base, self.cur, caption=caption)
            return
        # a step with nothing to point at (a table opened through the API) shows the page it led to
        if not any(leg.get('target') or leg.get('pointer') for leg in legs):
            self.compose(after, self.cur, caption=caption).save(still_path)
            yield from self.hold_frames(self.hold, after, self.cur, caption=caption)
            return
        # the stops of the step in order — one, or a menu's group, each item, the leaf: the pointer
        # goes where the page says it went on each, over the page as it was then, and every press
        # ripples under its tip
        still = None
        zoom, zoom_at = 1.0, None
        lit, target = None, None
        for n_leg, leg in enumerate(legs):
            shot = self._image(leg['shot'])
            target = leg.get('target')
            lit = self.spotlight(shot, target) if target else shot
            events = leg.get('pointer', [])
            lit_up = False
            still = None
            i = 0
            while i < len(events):
                e = events[i]
                at = (e['x'], e['y'])
                if e['type'] == 'up' or (e['type'] == 'move' and e.get('held')):
                    i += 1
                    continue
                if zoom > 1 and math.hypot(at[0] - zoom_at[0], at[1] - zoom_at[1]) > 2:
                    yield from self.zoom_frames(lit, zoom_at, zoom, 1.0, caption)
                    zoom = 1.0
                yield from self.travel_frames(shot, lit, at, caption, fade=not lit_up)
                lit_up = True
                if e['type'] == 'move':
                    i += 1
                    continue
                path = []
                j = i + 1
                while j < len(events) and events[j]['type'] == 'move' and events[j].get('held'):
                    path.append(events[j])
                    j += 1
                # an icon-sized target is zoomed into, slowly enough to read as a zoom (a drag is not:
                # it leaves the zoomed picture); a press rests under the pointer before it lands
                if is_small(target) and zoom == 1.0 and self.zoom > 1 and not path:
                    yield from self.zoom_frames(lit, at, 1.0, self.zoom, caption)
                    zoom, zoom_at = self.zoom, at
                if (e.get('count') or 1) == 1:
                    yield from self.hold_frames(0.4, lit, at, zoom=zoom, anchor=zoom_at, caption=caption)
                yield from self.press_frames(lit, at, e, target, step, index, caption, zoom, zoom_at)
                still = self.compose(lit, at, press=(at, 0.3, e.get('button', 'left')), caption=caption)
                if path:
                    lit = yield from self.drag_frames(lit, target, path, e.get('button', 'left'), caption)
                    still = self.compose(lit, self.cur, held=e.get('button', 'left'), caption=caption)
                i = j
            if not lit_up and target:
                yield from self.light_frames(shot, lit, caption)
            if n_leg < len(legs) - 1:
                yield from self.hold_frames(0.4, lit, self.cur, zoom=zoom, anchor=zoom_at, caption=caption)
                if zoom > 1:
                    yield from self.zoom_frames(lit, zoom_at, zoom, 1.0, caption)
                    zoom = 1.0
        # the picture of the step for steps.md: the last stop lit, and where the pointer acted on it
        (still if still is not None else self.compose(lit, self.cur, caption=caption)).save(still_path)
        yield from self.hold_frames(0.3, lit, self.cur, zoom=zoom, anchor=zoom_at, caption=caption)
        # the page after the step: the light moves over to it, zooming out
        lit_after = self.spotlight(after, target) if target else after
        n = int((self.zoom_time if zoom > 1 else 0.35) * fps)
        for i in range(n):
            t = ease((i + 1) / n)
            z = zoom + (1 - zoom) * t
            base = Image.blend(lit, lit_after, min(1.0, t * 2)) if t < 0.5 else Image.blend(lit_after, after, (t - 0.5) * 2)
            yield self.frame(base, self.cur, zoom=z, anchor=zoom_at, caption=caption)
        yield from self.hold_frames(self.hold, after, self.cur, caption=caption)

    def frames(self):
        if not self.steps:
            return
        yield from self.title_frames(self._image(self.steps[0]['before']))
        total = len(self.steps)
        for i, step in enumerate(self.steps, 1):
            yield from self.step_frames(step, i, total)
        yield from self.hold_frames(0.8, self._image(self.steps[-1]['after']), self.cur)

    def write_markdown(self):
        lines = [f"# {self.manifest['scenario']}", '']
        if self.manifest.get('description'):
            lines += [self.manifest['description'], '']
        for i, step in enumerate(self.steps, 1):
            mark = '\u2713 ' if step['kind'] == 'check' else ''
            lines += [f"{i}. {mark}{step['caption']}", '', f"   ![Step {i}](step-{i:02d}.png)", '']
        with open(os.path.join(self.folder, 'steps.md'), 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        # a guide filmed again with fewer steps would keep the pictures of the steps it lost
        n = len(self.steps) + 1
        while os.path.exists(os.path.join(self.folder, f'step-{n:02d}.png')):
            os.remove(os.path.join(self.folder, f'step-{n:02d}.png'))
            n += 1

    def write_audit(self):
        """Measures every press on the frame the video shows: the centre of what the mark changed,
        and the tip of the pointer (the top left of what the arrow changed), against the press; and
        whether the press is on the element the stop lit. audit.png shows each, cut out of its frame
        as the video has it, then again with ticks aimed at the press."""
        report = []
        tiles = []
        for p in self.presses:
            image, marked, bare = p.pop('images')
            vx, vy = p['video']
            mark = ImageChops.difference(marked, bare).convert('L').point(lambda v: 255 if v > 12 else 0)
            tip = ImageChops.difference(image, marked).convert('L').point(lambda v: 255 if v > 24 else 0)
            p['mark_off'] = self.centre_off(mark, vx, vy)
            box = tip.getbbox()
            p['tip_off'] = None if box is None else round(math.hypot(box[0] - vx, box[1] - vy), 1)
            p['video'] = [round(vx, 1), round(vy, 1)]
            report.append(p)
            tiles.append(self.audit_tile(image, vx, vy, p))
        with open(os.path.join(self.folder, 'audit.json'), 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2)
        if tiles:
            cols = 3
            tw, th = tiles[0].size
            sheet = Image.new('RGB', (cols * tw, math.ceil(len(tiles) / cols) * th), (255, 255, 255))
            for i, tile in enumerate(tiles):
                sheet.paste(tile, ((i % cols) * tw, (i // cols) * th))
            sheet.save(os.path.join(self.folder, 'audit.png'))
        return report

    def centre_off(self, mask, vx, vy):
        """How far the centre of the mark is from the press; a ring cut by the page's edge has a
        centre elsewhere, so only a mark wholly on the page is measured."""
        box = mask.getbbox()
        if box is None:
            return None
        top = self.bar_h * self.video[1] / self.out_h
        if box[0] <= 0 or box[1] <= top + 1 or box[2] >= mask.width or box[3] >= mask.height:
            return 'clipped'
        # the ring is round: its centre is the middle of the pixels it changed
        return round(math.hypot((box[0] + box[2] - 1) / 2 - vx, (box[1] + box[3] - 1) / 2 - vy), 1)

    def audit_tile(self, image, vx, vy, p):
        """A press as the video shows it, and beside it the same with ticks aimed at the press."""
        cw, ch = 220, 150
        x0, y0 = int(vx - cw / 2), int(vy - ch / 2)
        crop = Image.new('RGB', (cw, ch), (255, 255, 255))
        crop.paste(image.crop((max(0, x0), max(0, y0), min(image.width, x0 + cw), min(image.height, y0 + ch))),
                   (max(0, -x0), max(0, -y0)))
        aimed = crop.copy()
        d = ImageDraw.Draw(aimed)
        cx, cy = vx - x0, vy - y0
        for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            d.line([(cx + dx * 34, cy + dy * 34), (cx + dx * 60, cy + dy * 60)], fill=(255, 0, 255), width=2)
        label_h = 34
        tile = Image.new('RGB', (2 * cw + 30, ch + label_h + 12), (255, 255, 255))
        tile.paste(crop, (6, label_h))
        tile.paste(aimed, (cw + 18, label_h))
        d = ImageDraw.Draw(tile)
        font = load_font(12)
        count = ' x2' if p['count'] > 1 else ''
        verdict = '' if p['inside'] is not False else '  NOT ON THE LIT ELEMENT'
        d.text((6, 2), f"{p['step']}. {p['caption']}"[:62], font=font, fill=(20, 20, 20))
        d.text((6, 17), f"{p['button']}{count} on {p['on'] or '?'}  mark {p['mark_off']} px, tip {p['tip_off']} px{verdict}"[:70],
               font=font, fill=(160, 0, 0) if verdict else (90, 90, 90))
        return tile

    def encode(self, ffmpeg, gif):
        h264 = has_encoder(ffmpeg, 'libx264')
        video = os.path.join(self.folder, 'guide.mp4' if h264 else 'guide.webm')
        codec = ['-c:v', 'libx264', '-preset', 'veryfast', '-crf', '22', '-pix_fmt', 'yuv420p', '-movflags', '+faststart'] \
            if h264 else ['-c:v', 'libvpx', '-b:v', '2M', '-pix_fmt', 'yuv420p']
        cmd = [ffmpeg, '-y', '-loglevel', 'error', '-f', 'rawvideo', '-pix_fmt', 'rgb24', '-s', f'{self.video[0]}x{self.video[1]}',
               '-r', str(self.fps), '-i', '-', *codec, video]
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
        count = 0
        try:
            for frame in self.frames():
                proc.stdin.write(frame.tobytes())
                count += 1
        finally:
            proc.stdin.close()
            code = proc.wait()
        if code != 0:
            raise RuntimeError(f'ffmpeg failed encoding {video}')
        outputs = [video]
        if gif and self.steps:
            out = os.path.join(self.folder, 'guide.gif')
            # a docs GIF: ten frames a second at 880 px, no dithering (which defeats GIF's run-length
            # compression on flat UI colors) — a seven-step guide lands around a megabyte
            width = min(880, self.video[0]) & ~1
            scale = f'fps=10,scale={width}:-2:flags=lanczos'
            palette = os.path.join(self.folder, 'palette.png')
            subprocess.run([ffmpeg, '-y', '-loglevel', 'error', '-i', video, '-vf',
                            f'{scale},palettegen=max_colors=96:stats_mode=diff', palette], check=True)
            with_press_colours(palette)
            subprocess.run([ffmpeg, '-y', '-loglevel', 'error', '-i', video, '-i', palette, '-lavfi',
                            f'[0:v]{scale}[v];[v][1:v]paletteuse=dither=none:diff_mode=rectangle', out], check=True)
            os.remove(palette)
            thumb = os.path.join(self.folder, 'guide-thumb.png')
            shutil.copyfile(os.path.join(self.folder, 'step-01.png'), thumb)
            outputs += [out, thumb]
        return outputs, count


def with_press_colours(path):
    """ffmpeg's GIF palette with the press colours put in: a ripple is half a second of the video,
    too little for the palette to keep its yellow or green, and the GIF would show it grey. Each
    colour goes in as it looks over a white, a grey and a dark page; the palette keeps its
    transparent cell last."""
    cells = list(Image.open(path).convert('RGBA').getdata())
    clear = [c for c in cells if c[3] == 0][:1]
    colours = [c for c in cells if c[3] == 255]
    for colour in (LEFT_PRESS, RIGHT_PRESS):
        for under, weights in (((255, 255, 255), (1.0, 0.8, 0.6, 0.4)), ((150, 150, 150), (0.8, 0.5)), ((0, 0, 0), (0.6,))):
            colours += [tuple(round(c * a + u * (1 - a)) for c, u in zip(colour, under)) + (255,) for a in weights]
    colours = list(dict.fromkeys(colours))[:256 - len(clear)]
    colours += [colours[-1]] * (256 - len(clear) - len(colours)) + clear
    palette = Image.new('RGBA', (16, 16))
    palette.putdata(colours)
    palette.save(path)


def selftest(folder):
    os.makedirs(folder, exist_ok=True)
    w, h = 640, 400
    names = ['01-before', '01-after', '02-before', '02-after', '03-before', '03-drag1', '03-drag2', '03-after',
             '04-before', '04-after', '05-before', '05-hop2', '05-after', '06-before', '06-after']
    for name in names:
        img = Image.new('RGB', (w, h), (240, 243, 247) if name == '01-before' else (224, 236, 250))
        d = ImageDraw.Draw(img)
        d.rectangle((40, 40, 120, 64), fill=(32, 131, 213))
        d.rectangle((600, 8, 620, 28), fill=(90, 90, 90))
        if name.startswith('03-drag'):
            d.rectangle((200, 100, 200 + 60 * int(name[7]), 100 + 40 * int(name[7])), fill=(180, 180, 180))
        if name == '05-hop2':
            d.rectangle((300, 150, 460, 250), fill=(255, 255, 255), outline=(150, 150, 150))
            d.text((312, 190), 'Item', fill=(40, 40, 40))
        d.text((300, 320), name, fill=(40, 40, 40))
        img.save(os.path.join(folder, name + '.png'))

    def press(x, y, button='left', count=1):
        return [{'type': 'down', 'x': x, 'y': y, 'button': button, 'count': count, 'on': 'selftest'},
                {'type': 'up', 'x': x, 'y': y, 'button': button}]

    def step(i, text, caption, kind, legs, keys=()):
        return {'index': i, 'line': i, 'keyword': 'Then' if kind == 'check' else 'When', 'text': text, 'caption': caption,
                'kind': kind, 'before': f'0{i}-before.png', 'after': f'0{i}-after.png', 'legs': legs,
                'keys': list(keys), 'typed': '', 'ms': 10}
    button = {'x': 40, 'y': 40, 'width': 80, 'height': 24}
    icon = {'x': 600, 'y': 8, 'width': 20, 'height': 20}
    manifest = {'feature': 'Self test', 'scenario': 'Two synthetic steps', 'description': '', 'tags': [],
                'viewport': {'width': w, 'height': h}, 'steps': [
                    # the press is off the button's centre, where the mark and the tip must be
                    step(1, 'user clicks on the blue button', 'Click on the blue button', 'action',
                         [{'shot': '01-before.png', 'target': button, 'pointer': [{'type': 'move', 'x': 52, 'y': 47}] + press(52, 47)}]),
                    step(2, 'the page should be blue', 'The page is blue', 'check', [], ['Enter']),
                    step(3, 'user drags across the page', 'Drag across the page', 'action',
                         [{'shot': '03-before.png', 'target': {'x': 180, 'y': 80, 'width': 300, 'height': 200}, 'pointer': [
                             {'type': 'move', 'x': 200, 'y': 100}, {'type': 'down', 'x': 200, 'y': 100, 'button': 'left', 'count': 1},
                             {'type': 'move', 'x': 230, 'y': 120, 'held': True},
                             {'type': 'move', 'x': 260, 'y': 140, 'held': True, 'shot': '03-drag1.png'},
                             {'type': 'move', 'x': 320, 'y': 180, 'held': True, 'shot': '03-drag2.png'},
                             {'type': 'up', 'x': 320, 'y': 180, 'button': 'left'}]}]),
                    # a right-click on an icon at the page's edge: zoomed into, green
                    step(4, 'user right-clicks on the grey icon', 'Right-click on the grey icon', 'action',
                         [{'shot': '04-before.png', 'target': icon, 'pointer': press(611, 19, 'right')}]),
                    # a menu: the group hovered, then the item clicked on the page the hover opened
                    step(5, 'user picks "Item" from the menu', 'Pick "Item" from the menu', 'action',
                         [{'shot': '05-before.png', 'target': button, 'pointer': [{'type': 'move', 'x': 80, 'y': 52}]},
                          {'shot': '05-hop2.png', 'target': {'x': 300, 'y': 180, 'width': 160, 'height': 30},
                           'pointer': [{'type': 'move', 'x': 330, 'y': 195}] + press(330, 195)}]),
                    step(6, 'user double-clicks on the blue button', 'Double-click on the blue button', 'action',
                         [{'shot': '06-before.png', 'target': button, 'pointer': press(100, 58) + press(100, 58, count=2)}])]}
    with open(os.path.join(folder, 'steps.json'), 'w', encoding='utf-8') as f:
        json.dump(manifest, f)


def main(argv):
    args = list(argv)

    def opt(name, default=None):
        if name in args:
            i = args.index(name)
            value = args[i + 1]
            del args[i:i + 2]
            return value
        return default

    def flag(name):
        if name in args:
            args.remove(name)
            return True
        return False

    gif = flag('--gif')
    quiet = flag('--quiet')
    fps = int(opt('--fps', 30))
    hold = float(opt('--hold', 1.4))
    travel = float(opt('--travel', 0.7))
    zoom = float(opt('--zoom', 1.8))
    zoom_time = float(opt('--zoom-time', 0.8))
    ffmpeg = find_ffmpeg(opt('--ffmpeg'))
    test_dir = opt('--selftest')
    if test_dir:
        selftest(test_dir)
        args.append(test_dir)
        gif = True
    if len(args) != 1:
        print(__doc__, file=sys.stderr)
        return 2
    if not ffmpeg:
        print('guide-render: no ffmpeg — py -m pip install imageio-ffmpeg, or set FFMPEG', file=sys.stderr)
        return 3
    renderer = Renderer(args[0], fps, hold, travel, zoom, zoom_time)
    renderer.write_markdown()
    outputs, count = renderer.encode(ffmpeg, gif)
    report = renderer.write_audit()
    # a press marked or pointed at away from where it landed is a rendering error; one outside the
    # element its stop lit is the recording's (the step lit one element and clicked another)
    off = [p for p in report if any(isinstance(v, float) and v > AUDIT_TOLERANCE_PX for v in (p['mark_off'], p['tip_off']))
           or p['mark_off'] is None or p['tip_off'] is None]
    astray = [p for p in report if p['inside'] is False]
    if not quiet:
        print(f'{count} frames, {count / fps:.1f} s: ' + ', '.join(outputs))
        print(f"audit: {len(report)} presses, {len(off)} marked or pointed away from the press, "
              f"{len(astray)} off the element lit for them — audit.png")
    for p in off + astray:
        print(f"guide-render: step {p['step']} \"{p['caption']}\": {p['button']} press at ({p['x']:.0f}, {p['y']:.0f}) "
              f"on {p['on']}: mark {p['mark_off']} px, tip {p['tip_off']} px, on the lit element: {p['inside']}", file=sys.stderr)
    return 4 if off and test_dir else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
