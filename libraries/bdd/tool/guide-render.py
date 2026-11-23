"""Renders a guide manifest (`steps.json` + the PNGs a `BDD_GUIDE` run wrote) into guide.mp4,
step-NN.png (the annotated picture of every step), steps.md, audit.png + audit.json, and with --gif
also guide.gif (made from the lossless frames), guide-thumb.png (the result) and guide-sheet.png
(twelve frames of the GIF, for review). Pillow does the drawing, ffmpeg the encoding (FFMPEG env,
PATH, the imageio-ffmpeg package, or Playwright's own build, which can only write WebM). The video
and the stills are at most 1600×1000 with the caption strip: a 1080p page is downsampled, a small
window is filmed at twice its density (the manifest's `scale`).

--help makes a help-page GIF: no title card, caption strip or dimming, no zoom, the set-up (the
Given steps) left out, a pause before the first action and a longer one on the result.

--demo makes a demo clip, one of several cut into a longer video: always 1920×1080, a title card
(the scenario, the feature, the description, wrapped to fit), then the page with no caption strip,
dimming or zoom, the presses still marked. The set-up (Given steps) and the checks are left out but
for those the feature captioned; a check's result shows as the page settles on it. A note (`the video
says …`, or a caption) is a line over the page for as long as it takes to read. The pace is a little
quicker, and a step rests as long as its change takes to take in. review-sheet.png (the last frame of
every step) and timeline.json (when each step starts and ends) are written for a review.

--tutorial makes a demo clip that teaches (a short how-to of a series): the element each stop is about to press is
outlined over a lightly dimmed page, the pointer rests on it a moment longer, and the notes rise into place.

--reel <out.mp4> cuts demo clips (the folders a --demo rendering wrote) together, in the order given, into
one video; each clip is copied beside it as `NN - <scenario>.mp4`, and `<out>.md` lists them with where
each starts in the reel.

The pointer goes where the page says it went, and never skips: every frame's pointer continues
from the frame before, so an action starts where the one before it ended — a skip is an error,
not a frame. A press is marked where it landed, under the pointer's tip: a yellow ring for the
left button, a green one for the right. audit.png shows every press as the video shows it, and
audit.json how far the mark's centre and the pointer's tip are from the press.

    py tool/guide-render.py <scenario dir> [--gif] [--help | --demo | --tutorial] [--fps 30] [--hold 1.4] [--travel 0.7]
        [--zoom 1.8] [--zoom-time 0.8] [--gif-fps 12] [--colors 256] [--gif-width <px>] [--end-hold <s>]
        [--ffmpeg <path>] [--quiet]
    py tool/guide-render.py --reel <out.mp4> <demo clip dir>…   cuts the clips together in that order
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
    from PIL import Image, ImageChops, ImageDraw, ImageFilter, ImageFont
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
BEAT_S = 0.15   # how long a picture taken in the middle of a step (typing, a wheel turn) is shown
DEMO_W, DEMO_H = 1920, 1080   # a demo clip, whatever window it was filmed in
DEMO_KICKER = (110, 180, 235)
NOTE_BOX = (16, 20, 28)


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


def load_font(size, bold=False):
    faces = [('C:/Windows/Fonts/segoeuib.ttf', 0), ('C:/Windows/Fonts/arialbd.ttf', 0),
             ('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 0), ('/System/Library/Fonts/Helvetica.ttc', 1)] if bold else \
        [('C:/Windows/Fonts/segoeui.ttf', 0), ('C:/Windows/Fonts/arial.ttf', 0),
         ('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 0), ('/System/Library/Fonts/Helvetica.ttc', 0)]
    for path, index in faces:
        if os.path.exists(path):
            return ImageFont.truetype(path, size, index=index)
    try:
        return ImageFont.load_default(size)
    except TypeError:
        return ImageFont.load_default()


def wrap(d, text, font, width):
    """The text in lines no wider than `width`, broken between words (a word wider than the line is a line)."""
    lines, line = [], ''
    for word in text.split():
        candidate = (line + ' ' + word).strip()
        if line and d.textlength(candidate, font=font) > width:
            lines.append(line)
            line = word
        else:
            line = candidate
    return lines + [line] if line else lines


def balanced(d, text, font, width):
    """The lines of `text` at the narrowest width that still needs no more lines than `width` does: lines of
    even length, no word left alone on the last."""
    lines = wrap(d, text, font, width)
    if len(lines) < 2:
        return lines
    lo, hi = max(d.textlength(w, font=font) for w in text.split()), width
    while hi - lo > 2:
        mid = (lo + hi) / 2
        if len(wrap(d, text, font, mid)) <= len(lines):
            hi = mid
        else:
            lo = mid
    return wrap(d, text, font, hi)


def fit_text(d, text, size, least, width, max_lines, bold=False):
    """The largest font from `size` down to `least` that sets the text in at most `max_lines` balanced lines
    of `width`; at the smallest size every line is kept, however many — a title card never cuts its text."""
    for s in range(size, least - 1, -2):
        font = load_font(s, bold)
        lines = balanced(d, text, font, width)
        if len(lines) <= max_lines:
            return font, lines
    return font, lines


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


def along(points, lengths, d):
    """The point `d` along the way through `points` (`lengths`: of each of its segments)."""
    for a, b, length in zip(points, points[1:], lengths):
        if d <= length:
            f = d / length if length > 0 else 1.0
            return a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f
        d -= length
    return points[-1]


# a step names Control and the gesture goes through ControlOrMeta, which a Mac's browser receives as Meta: a cap names the
# key the step meant, since a film made on a Mac is watched on every platform
KEY_NAMES = {'ControlOrMeta': 'Ctrl', 'Control': 'Ctrl', 'Meta': 'Ctrl', 'Escape': 'Esc', 'Delete': 'Del', ' ': 'Space'}


def key_label(keys):
    """Keys as a key cap names them: `ControlOrMeta+Z` is `Ctrl + Z`, `ArrowRight` is `Right`, and a key pressed
    again and again is named once, with the count."""
    runs = []
    for key in keys:
        if runs and runs[-1][0] == key:
            runs[-1][1] += 1
        else:
            runs.append([key, 1])
    named = (' + '.join(KEY_NAMES.get(p, p[5:] if p.startswith('Arrow') else p) for p in key.split('+')) for key, _ in runs)
    return '   '.join(name + (f' ×{n}' if n > 1 else '') for name, (_, n) in zip(named, runs))


class Renderer:
    def __init__(self, folder, fps, hold, travel, zoom, zoom_time=0.8, help=False, demo=False, tutorial=False):
        self.folder = folder
        self.fps = fps
        self.help = help
        # a tutorial is a demo clip that teaches: the element about to be pressed is outlined, and its notes rise
        # into place
        self.tutorial = tutorial and demo
        self.demo = demo
        self.hold = hold
        self.travel = travel
        self.zoom = 1.0 if help or demo else zoom
        self.zoom_time = zoom_time
        # a demo keeps a little quicker pace wherever a guide waits for its caption to be read: the rest
        # before a press, the ripple, a stop on a path's way, the step's picture before its result, a typing beat
        self.rest_s, self.pulse_s, self.leg_s, self.cut_s, self.beat_s = \
            (0.2, 0.42, 0.22, 0.12, 0.1) if demo else (0.4, PULSE_S, 0.4, 0.3, BEAT_S)
        # a tutorial is followed, not watched: the outlined element is up long enough to find before the press and
        # stays a moment after it, a menu's stops are read, and the pointer moves at a guide's pace
        if self.tutorial:
            self.rest_s, self.cut_s, self.leg_s = 0.75, 0.5, 0.4
        # but a menu's items (and a list's) are read as the pointer goes: a stop of a menu holds its opened group, the
        # rest before its leaf is pressed is longer, and an opened list stays up longer still
        self.menu_s, self.menu_rest_s, self.list_s = (0.38, 0.35, 0.45) if demo else (0.4, 0.4, 0.5)
        if self.tutorial:
            self.menu_s, self.menu_rest_s, self.list_s = 0.6, 0.6, 0.7
        with open(os.path.join(folder, 'steps.json'), encoding='utf-8') as f:
            self.manifest = json.load(f)
        # the page records CSS pixels; the pictures have `scale` times as many. A guide framed to a part
        # of the page (a sketcher's dialog: `frameTo` in guide.ts) shows that part alone: every picture is
        # cut to it, and every point and box is moved into it
        scale = self.manifest.get('scale') or 1
        self.frame_box = self.manifest.get('frame')
        if self.frame_box:
            self.frame_box = {k: round(v * scale) for k, v in self.frame_box.items()}
        fx, fy = (self.frame_box['x'], self.frame_box['y']) if self.frame_box else (0, 0)
        # the keyword each step stands under: an And or But takes the one before it
        kind = 'Given'
        for s in self.manifest['steps']:
            if s['keyword'] not in ('And', 'But', '*'):
                kind = s['keyword']
            s['type'] = kind
        self.steps = [s for s in self.manifest['steps'] if self.shown(s)]
        first_name = self.steps[0]['before'] if self.steps else self.manifest['steps'][0]['before']
        # a demo is drawn at the video's size: each picture is shrunk once, as it is read, and every point
        # and box with it
        self.k = 1.0
        self.cache = {}
        if demo:
            raw = self._raw(first_name)
            self.k = min(1.0, DEMO_W / raw.width, DEMO_H / raw.height)
        k = self.k

        def place(box):
            box.update(x=(box['x'] * scale - fx) * k, y=(box['y'] * scale - fy) * k, width=box['width'] * scale * k,
                       height=box['height'] * scale * k)
        for step in self.manifest['steps']:
            for box in [step.get('target')] + [leg.get('target') for leg in step.get('legs', [])]:
                if box:
                    place(box)
            for leg in step.get('legs', []):
                for p in leg.get('pointer', []):
                    p['x'] = (p['x'] * scale - fx) * k
                    p['y'] = (p['y'] * scale - fy) * k
                    if p.get('box'):
                        place(p['box'])
        self.scale = scale
        first = self._image(first_name)
        self.w, self.h = first.width & ~1, first.height & ~1
        # marks and the pointer keep their size on the page: twice the pixels where the page has twice the density
        self.s = max(self.h / 1080, scale * k)
        self.big = load_font(max(20, self.h // 18))
        self.cursor = cursor_sprite(max(1.0, self.h / 700, 1.2 * scale * k))
        # the caption sits above the page: nothing of the page is covered, and a player's timeline
        # lands on the bottom edge. The page keeps the pixels it was filmed at (every box, pointer
        # and zoom is in them); the strip is as tall as makes page plus strip the video's aspect
        # (120 px over a 1080p page), and the composed frame is downsampled to the video size —
        # never upsampled, so a small capture stays as it is. A demo has no strip: its page is
        # centred in the 1080p frame
        self.bar_h = 0 if help or demo else max(int(self.h * 0.085), round(self.w * VIDEO_H / VIDEO_W) - self.h) & ~1
        self.out_h = self.h + self.bar_h
        if demo:
            self.video = (DEMO_W, DEMO_H)
            self.pad = ((DEMO_W - self.w) // 2, (DEMO_H - self.h) // 2)
        else:
            scale = min(1.0, VIDEO_W / self.w, VIDEO_H / self.out_h)
            self.video = (int(self.w * scale) & ~1, int(self.out_h * scale) & ~1)
            self.pad = (0, 0)
        self.font = load_font(max(14, self.bar_h // 3))
        self.small = load_font(max(12, self.bar_h // 4))
        self.cur = (self.w / 2, self.h / 2)
        self.drawn = None
        self.frame_no = 0
        self.presses = []
        # the line a demo shows over the page and its key cap, while they last; the last frame composed without them
        self.note = None
        self.cap = None
        self.plain = None
        # the page the last frame was composed from, before its pointer, marks and note
        self.base_now = None
        # where each part of the video starts (the title, every step, the end): timeline.json, review-sheet.png
        self.marks = []
        if demo:
            self.settle_results()

    def shown(self, step):
        """What the video shows: never the set-up; in a help GIF or a demo no Given step either; a demo
        shows a check only when the feature captioned it, and a help GIF no note."""
        if step['kind'] == 'setup':
            return False
        if step['kind'] == 'note':
            return not self.help
        if (self.help or self.demo) and step['type'] == 'Given':
            return False
        return not (self.demo and step['kind'] == 'check' and not step.get('captioned'))

    def settle_results(self):
        """A demo leaves its checks and its set-up out, not what they waited for or arranged: the page an
        action led to once the steps after it passed (a column computed, an analysis drawn, the new column
        sized to show its sequences) is the action's result."""
        steps = self.manifest['steps']
        shown = {id(s) for s in self.steps}
        for i, step in enumerate(steps):
            if id(step) not in shown or step['kind'] not in ('action', 'note'):
                continue
            for later in steps[i + 1:]:
                # a table opened off camera is the next scene, not this step's result
                if id(later) in shown or later.get('view', step.get('view')) != step.get('view'):
                    break
                step['settled'] = later['after']

    def is_given(self, step):
        """A set-up step: a Given, or an And or But after one."""
        return step['type'] == 'Given'

    def _raw(self, name):
        img = Image.open(os.path.join(self.folder, name)).convert('RGB')
        if self.frame_box:
            b = self.frame_box
            img = img.crop((b['x'], b['y'], b['x'] + b['width'], b['y'] + b['height']))
        return img

    def _image(self, name):
        if name in self.cache:
            return self.cache[name]
        img = self._raw(name)
        if self.k < 1:
            img = img.resize((max(2, round(img.width * self.k)), max(2, round(img.height * self.k))), Image.LANCZOS)
        if not hasattr(self, 'w'):
            return img
        if img.width != self.w or img.height != self.h:
            img = img.crop((0, 0, self.w, self.h))
        # a picture is read again for every stop and hold of its step: the last few are kept, never drawn on
        if len(self.cache) >= 24:
            self.cache.pop(next(iter(self.cache)))
        self.cache[name] = img
        return img

    def spotlight(self, img, box, color=ACCENT):
        """The frame with everything but the box dimmed and the box outlined — none in a help GIF or a demo; a
        tutorial outlines the box over a lighter dimming."""
        if self.tutorial:
            return self.highlight(img, box)
        if self.help or self.demo:
            return img
        pad = 6
        x0, y0 = max(0, box['x'] - pad), max(0, box['y'] - pad)
        x1, y1 = min(self.w, box['x'] + box['width'] + pad), min(self.h, box['y'] + box['height'] + pad)
        if x1 - x0 < 2 or y1 - y0 < 2:
            return img
        layer = Image.new('RGBA', (self.w, self.h), (0, 0, 0, 105))
        d = ImageDraw.Draw(layer)
        d.rounded_rectangle((x0, y0, x1, y1), radius=6, fill=(0, 0, 0, 0))
        out = img.convert('RGBA')
        out.alpha_composite(layer)
        d = ImageDraw.Draw(out)
        d.rounded_rectangle((x0, y0, x1, y1), radius=6, outline=color + (255,), width=3)
        return out.convert('RGB')

    def highlight(self, img, box):
        """A tutorial's mark of where to press: the page dimmed a little round the box, the box outlined in the accent
        colour over a soft glow, so it reads on a white page and on a dark one."""
        # an icon is marked by a box larger than itself, so a status-bar toggle reads at a glance
        least = 44 * self.s
        cx, cy = box['x'] + box['width'] / 2, box['y'] + box['height'] / 2
        half_w = max(box['width'] / 2 + 5 * self.s, least / 2)
        half_h = max(box['height'] / 2 + 5 * self.s, least / 2)
        x0, y0 = max(0, cx - half_w), max(0, cy - half_h)
        x1, y1 = min(self.w - 1, cx + half_w), min(self.h - 1, cy + half_h)
        # a box the page had not scrolled into view yet when it was read: nothing of it to outline
        if x1 - x0 < 2 or y1 - y0 < 2:
            return img
        radius = round(6 * self.s)
        layer = Image.new('RGBA', (self.w, self.h), (0, 0, 0, 70))
        d = ImageDraw.Draw(layer)
        d.rounded_rectangle((x0, y0, x1, y1), radius=radius, fill=(0, 0, 0, 0))
        glow = Image.new('RGBA', (self.w, self.h), (0, 0, 0, 0))
        g = round(4 * self.s)
        ImageDraw.Draw(glow).rounded_rectangle((x0 - g, y0 - g, x1 + g, y1 + g), radius=radius + g,
                                               outline=ACCENT + (120,), width=2 * g)
        glow = glow.filter(ImageFilter.GaussianBlur(radius=g))
        out = img.convert('RGBA')
        out.alpha_composite(layer)
        out.alpha_composite(glow)
        ImageDraw.Draw(out).rounded_rectangle((x0, y0, x1, y1), radius=radius, outline=ACCENT + (255,),
                                              width=max(3, round(3 * self.s)))
        return out.convert('RGB')

    def caption_layer(self, step, index, total):
        if self.bar_h == 0:
            return None
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
        # the keys the page's keyboard pressed, else the ones the step pressed some other way (Escape on a dialog)
        labels = [key.replace('ControlOrMeta', 'Ctrl') for key in (step.get('keys') or step.get('pressed') or [])[:4]]
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
        if self.demo:
            if frame.size == self.video:
                return frame
            canvas = Image.new('RGB', self.video, NOTE_BOX)
            canvas.paste(frame, self.pad)
            return canvas
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
        if self.demo:
            return x + self.pad[0], y + self.pad[1]
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
        self.base_now = base
        self.plain = self.compose(base, cursor, **kw)
        return self.with_note(self.plain)

    def hold_frames(self, seconds, base, cursor=None, **kw):
        n = int(seconds * self.fps)
        if n <= 0:
            return
        yield self.frame(base, cursor, **kw)
        still = self.plain
        for _ in range(n - 1):
            self.frame_no += 1
            yield self.with_note(still)

    def title_frames(self, base):
        """The scenario's name over the first page, the feature's below it, each wrapped to fit."""
        layer = Image.new('RGBA', (self.w, self.h), (10, 14, 20, 150))
        frame = Image.alpha_composite(base.convert('RGBA'), layer).convert('RGB')
        d = ImageDraw.Draw(frame)
        width = self.w * 0.9
        big, title = fit_text(d, self.manifest['scenario'], self.big.size, max(14, round(self.big.size * 0.6)), width, 2)
        small, sub = fit_text(d, self.manifest['feature'], self.font.size, max(12, round(self.font.size * 0.7)), width, 2)
        big_lead, small_lead, gap = big.size * 1.2, small.size * 1.3, big.size * 0.5
        y = self.h * 0.47 - (len(title) * big_lead + gap + len(sub) * small_lead) / 2
        for line in title:
            d.text((self.w / 2, y + big_lead / 2), line, font=big, fill=(255, 255, 255), anchor='mm')
            y += big_lead
        y += gap
        for line in sub:
            d.text((self.w / 2, y + small_lead / 2), line, font=small, fill=(200, 206, 214), anchor='mm')
            y += small_lead
        yield from self.hold_frames(1.8, frame)

    def demo_card(self, base):
        """A demo's title card: the page it starts on, blurred and darkened, under the feature's name (small,
        above), the scenario's (large) and the description — each wrapped to fit, never cut."""
        vw, vh = self.video
        page = self.compose(base).filter(ImageFilter.GaussianBlur(radius=vh / 120))
        card = Image.alpha_composite(page.convert('RGBA'), Image.new('RGBA', self.video, (10, 14, 20, 190))).convert('RGB')
        d = ImageDraw.Draw(card)
        width = int(vw * 0.7)
        title = self.manifest['scenario']
        kicker = self.manifest.get('feature') or ''
        body = self.manifest.get('description') or ''
        # (font, lines, colour, leading, space after) from the top; a short rule between the title and the description
        blocks = []
        if kicker and kicker.lower() != title.lower():
            blocks.append(fit_text(d, kicker.upper(), round(vh * 0.022), round(vh * 0.017), width, 2, bold=True) +
                          (DEMO_KICKER, 1.4, round(vh * 0.02)))
        blocks.append(fit_text(d, title, round(vh * 0.056), round(vh * 0.034), width, 3, bold=True) +
                      ((255, 255, 255), 1.2, round(vh * 0.05)))
        if body:
            blocks.append(fit_text(d, body, round(vh * 0.029), round(vh * 0.019), width, 6) + ((212, 220, 230), 1.45, 0))
        height = sum(len(lines) * font.size * lead + space for font, lines, _, lead, space in blocks[:-1]) + \
            len(blocks[-1][1]) * blocks[-1][0].size * blocks[-1][3]
        y = (vh - height) / 2
        for i, (font, lines, colour, lead, space) in enumerate(blocks):
            for line in lines:
                d.text((vw / 2, y + font.size * lead / 2), line, font=font, fill=colour, anchor='mm')
                y += font.size * lead
            if i == len(blocks) - 2 and body:
                rule = round(vh * 0.035)
                d.rectangle((vw / 2 - rule, y + space / 2 - 2, vw / 2 + rule, y + space / 2 + 1), fill=DEMO_KICKER)
            y += space
        return card

    def demo_title_frames(self, base):
        """The card, faded in from dark (clips cut together meet in the dark), for as long as its words take to
        read at a glance (about three seconds: the card is the clip's headline, the clip shows the rest), then a
        dissolve into the page."""
        card = self.demo_card(base)
        dark = Image.new('RGB', self.video, NOTE_BOX)
        n = int(0.3 * self.fps)
        for i in range(n):
            self.frame_no += 1
            yield Image.blend(dark, card, ease((i + 1) / n))
        words = len((self.manifest['scenario'] + ' ' + (self.manifest.get('description') or '')).split())
        for _ in range(int(min(3.6, max(2.4, 1.0 + 0.1 * words)) * self.fps)):
            self.frame_no += 1
            yield card
        page = self.compose(base)
        n = int(0.45 * self.fps)
        for i in range(n):
            self.frame_no += 1
            yield Image.blend(card, page, ease((i + 1) / n))

    def reading(self, text, least=2.2, most=7.0):
        """Seconds to read a line of text."""
        return min(most, max(least, 1.2 + 0.3 * len(text.split())))

    def start_note(self, text, ahead=(), before=None, at=None, view=None):
        """A line over a demo's page from the next frame on, for as long as it takes to read: at the bottom or the
        top — where the feature pins it (`at`), else where it covers least of what the steps around it are about
        (`covers`): the step before it, which the note explains, and the steps it is read over. It belongs to the
        view it was put up in (`end_stale_note`)."""
        vw, vh = self.video
        d = ImageDraw.Draw(Image.new('RGB', (1, 1)))
        size = round(vh * 0.027)
        font, lines = fit_text(d, text, size, round(size * 0.75), int(vw * 0.6), 3)
        lead = font.size * 1.32
        px, py = round(font.size * 0.9), round(font.size * 0.55)
        bw = int(max(d.textlength(line, font=font) for line in lines) + 2 * px)
        bh = int(len(lines) * lead + 2 * py)
        layer = Image.new('RGBA', (bw, bh), (0, 0, 0, 0))
        ld = ImageDraw.Draw(layer)
        ld.rounded_rectangle((0, 0, bw - 1, bh - 1), radius=round(font.size * 0.45), fill=NOTE_BOX + (228,))
        for i, line in enumerate(lines):
            ld.text((bw / 2, py + (i + 0.5) * lead), line, font=font, fill=(255, 255, 255, 255), anchor='mm')
        x = (vw - bw) // 2
        ys = {'bottom': vh - bh - round(vh * 0.07), 'top': round(vh * 0.075)}
        steps = [s for s in (before, *ahead) if s is not None]
        place = at if at in ys else min(ys, key=lambda p: self.covers((x, ys[p], bw, bh), steps))
        self.note = {'rgb': layer.convert('RGB'), 'alpha': layer.getchannel('A'), 'at': (x, ys[place]),
                     'total': int(self.reading(text) * self.fps), 'shown': 0, 'view': view}

    def covers(self, rect, steps):
        """How much of what the steps are about a note at `rect` (video pixels) would cover: what they act on and
        where the pointer goes (their stops' targets, the pointer's way, as boxes) above all, then what they change
        on the page (a tooltip, a highlight) — counted only where it lies under the note. A step that brings a page
        of its own (a new view) changes everything and says nothing about where."""
        x0, y0 = rect[0] - self.pad[0], rect[1] - self.pad[1]
        x1, y1 = x0 + rect[2], y0 + rect[3]
        r = 28 * self.s
        score = 0.0
        for step in steps:
            boxes = [leg['target'] for leg in step.get('legs', []) if leg.get('target')]
            boxes += [{'x': p['x'] - r, 'y': p['y'] - r, 'width': 2 * r, 'height': 2 * r}
                      for leg in step.get('legs', []) for p in leg.get('pointer', []) if p['type'] in ('move', 'down')]
            for b in boxes:
                w = min(x1, b['x'] + b['width']) - max(x0, b['x'])
                h = min(y1, b['y'] + b['height']) - max(y0, b['y'])
                if w > 0 and h > 0:
                    score += 1000 * w * h / max(1.0, min(b['width'] * b['height'], rect[2] * rect[3]))
            legs = step.get('legs', [])
            first = self._image(legs[0]['shot'] if legs else step['before'])
            mask = self.change_mask(first, self._image(step.get('settled') or step['after']))
            if mask is not None:
                k = mask.width / self.w
                part = mask.crop((max(0, int(x0 * k)), max(0, int(y0 * k)), max(1, int(x1 * k)), max(1, int(y1 * k))))
                score += 100 * part.histogram()[255] / max(1, part.width * part.height)
        return score

    def change_mask(self, a, b):
        """Where two pictures of the page differ, on small grey copies; None when they differ almost everywhere."""
        size = (max(1, self.w // 8), max(1, self.h // 8))
        mask = ImageChops.difference(a.convert('L').resize(size, Image.BILINEAR), b.convert('L').resize(size, Image.BILINEAR))
        mask = mask.point(lambda v: 255 if v > 24 else 0)
        return None if mask.histogram()[255] > 0.4 * size[0] * size[1] else mask

    def before(self, step):
        """The shown step before this one: what a note put up now explains."""
        i = self.steps.index(step)
        return self.steps[i - 1] if i > 0 else None

    def ahead(self, step, n=2):
        """The step and the next shown ones a note will be read over."""
        i = self.steps.index(step)
        return self.steps[i:i + n + 1]

    def with_note(self, img):
        """The frame with the note and the key cap over it, each faded in and out; each counts the frame."""
        out = img
        for name, fade_s in (('note', 0.35 if self.tutorial else 0.25), ('cap', 0.1)):
            layer = getattr(self, name)
            if layer is None:
                continue
            fade = max(1, int(fade_s * self.fps))
            a = min(1.0, (layer['shown'] + 1) / fade, (layer['total'] - layer['shown']) / fade)
            x, y = layer['at']
            # a tutorial's note rises into place as it fades in
            if self.tutorial and name == 'note':
                y += round((1 - ease(min(1.0, (layer['shown'] + 1) / fade))) * 18 * self.video[1] / 1080)
            layer['shown'] += 1
            if layer['shown'] >= layer['total']:
                setattr(self, name, None)
            if out is img:
                out = img.copy()
            out.paste(layer['rgb'], (x, y), layer['alpha'] if a >= 1 else layer['alpha'].point(lambda v: int(v * a)))
        return out

    def change(self, a, b):
        """The share of the page two pictures differ in, read on small grey copies: a tooltip is a few
        hundredths, a dialog a tenth, a new view most of it."""
        size = (192, 108)
        diff = ImageChops.difference(a.convert('L').resize(size, Image.BILINEAR), b.convert('L').resize(size, Image.BILINEAR))
        return diff.point(lambda v: 255 if v > 20 else 0).histogram()[255] / (size[0] * size[1])

    def rest_for(self, change, least=0.3, most=1.3):
        """How long a demo rests on a change: a checkbox ticked hardly at all, a new view the longest."""
        return least + (most - least) * min(1.0, max(0.0, (change - 0.01) / 0.24))

    def demo_result(self, step, shown, after, lit=None):
        """The page a demo's step led to — once the steps after it settled it (a column computed, a viewer
        drawn, the new column sized), the page they left, cut to straight from the gesture: the page as the
        gesture left it (a file opened raw, a viewer docked before its set-up moved it) is never shown — held
        as long as its change takes to take in, a tooltip a hover brought up long enough to read, a dialog it
        opened long enough to read its inputs."""
        result = self._image(step['settled']) if step.get('settled') else after
        change = self.change(shown, result)
        # a tutorial rests a little longer on each result: the viewer follows along, and may pause to do the same
        hold = self.rest_for(change, 0.7, 2.2) if self.tutorial else self.rest_for(change, 0.3, 1.5)
        # and dissolves into it from the outlined page, so a changed page does not jump in
        if self.tutorial and lit is not None:
            n = int(0.35 * self.fps)
            for i in range(n):
                yield self.frame(Image.blend(lit, result, ease((i + 1) / n)), self.cur)
        events = [p['type'] for leg in step.get('legs', []) for p in leg.get('pointer', [])]
        # a hover brought up a tooltip: the bigger it is, the longer it takes to read
        if 'down' not in events and 'move' in events and change > 0.002:
            hold = max(hold, min(2.6, 1.2 + 7 * change))
        if step.get('dialog'):
            hold = max(hold, min(2.5, max(1.2, 0.9 + 0.12 * step['dialog']['inputs'])))
        yield from self.hold_frames(hold, result, self.cur)

    def finish_note(self):
        """A tutorial's note is read to its end before the next one comes up: the page holds where it is meanwhile."""
        if not (self.tutorial and self.note and self.base_now is not None):
            return
        left = (self.note['total'] - self.note['shown']) / self.fps
        if left > 0:
            yield from self.hold_frames(left, self.base_now, self.cur)

    def end_stale_note(self, step):
        """A note is about the view it was put up in: a step that brings another view (a file opened) fades it out as
        it goes, so the note is gone when the new view comes."""
        if self.note is not None and step.get('view') and self.note.get('view') and step['view'] != self.note['view']:
            self.note['total'] = min(self.note['total'], self.note['shown'] + max(1, int(0.2 * self.fps)))

    def frames_for(self, distance, least):
        """Enough frames that the eased way never moves the pointer more than MAX_STEP_PX in one."""
        return max(1, least, math.ceil(1.5 * distance / (MAX_STEP_PX * self.s)))

    def travel_frames(self, base, lit, dest, caption, fade=False, via=None):
        """The pointer eased from where the last frame left it to dest — a longer way takes longer —
        through the corner `via` when there is one, with the stop lighting up over the second half
        when `fade` (from `base`, the picture on screen as the pointer sets off)."""
        start = self.cur
        points = [start] + ([via] if via else []) + [dest]
        lengths = [math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(points, points[1:])]
        distance = sum(lengths)
        if distance < 0.5:
            if fade:
                yield from self.light_frames(base, lit, caption)
            return
        reach = 1600 * self.s
        n = self.frames_for(distance, round(self.fps * self.travel * (0.5 + min(distance, reach) / reach)))
        for i in range(n):
            t = ease((i + 1) / n)
            pos = along(points, lengths, t * distance)
            # a tutorial's outline comes up early in the way, so the eye finds the target before the pointer gets there
            mix = (ease(t / 0.35) if self.tutorial else ease((t - 0.5) / 0.5)) if fade and lit is not base else 1.0
            yield self.frame(lit if mix >= 1 else Image.blend(base, lit, mix) if mix > 0 else base, pos, caption=caption)

    def elbow(self, came_from, to, dest):
        """The corner of the way from an item of a menu into the flyout it opened beside it: along the item's own row
        first, then to the item of the flyout — a straight line would cross the items under the group, each of which
        opens a flyout of its own. None for a dropdown under a bar's group, or a way that starts elsewhere."""
        if came_from is None or to is None or not inside(came_from, self.cur[0], self.cur[1], pad=4):
            return None
        y = self.cur[1]
        if abs(dest[1] - y) < 2:
            return None
        inset = 10 * self.s
        right = came_from['x'] + came_from['width']
        if to['x'] >= right - 2:
            return min(right + inset, dest[0]), y
        if to['x'] + to['width'] <= came_from['x'] + 2:
            return max(came_from['x'] - inset, dest[0]), y
        return None

    def light_frames(self, base, lit, caption):
        """A stop lighting up where the pointer already is — one the step reached with the keyboard."""
        n = int(0.3 * self.fps)
        for i in range(n):
            yield self.frame(Image.blend(base, lit, ease((i + 1) / n)), self.cur, caption=caption)

    def arrive_frames(self, lit, arrived, zoom, anchor, caption):
        """What the pointer's arrival at a stop led to (the group opened under it, the item lit) cross-faded in where
        it came to rest, in a tenth of a second: the next stop's picture comes as the pointer gets there, not as it
        leaves — cause before effect."""
        n = max(1, round(0.1 * self.fps))
        for i in range(n):
            yield self.frame(Image.blend(lit, arrived, ease((i + 1) / n)), self.cur, zoom=zoom, anchor=anchor, caption=caption)

    def start_cap(self, label, seconds=0.6):
        """A key cap over a demo's page, at its bottom right, from the next frame on for `seconds`: the key a step
        pressed (Escape) or the one held through a press (Shift), which a demo has no caption strip to name."""
        vw, vh = self.video
        d = ImageDraw.Draw(Image.new('RGB', (1, 1)))
        font = load_font(round(vh * 0.024), bold=True)
        px, py = round(font.size * 0.8), round(font.size * 0.5)
        bw, bh = int(d.textlength(label, font=font) + 2 * px), int(font.size * 1.2 + 2 * py)
        depth, radius = max(3, round(vh * 0.004)), round(font.size * 0.35)
        layer = Image.new('RGBA', (bw, bh + depth), (0, 0, 0, 0))
        ld = ImageDraw.Draw(layer)
        # a key: its face over a darker edge
        ld.rounded_rectangle((0, depth, bw - 1, bh + depth - 1), radius=radius, fill=(6, 8, 12, 235))
        ld.rounded_rectangle((0, 0, bw - 1, bh - 1), radius=radius, fill=(40, 46, 58, 240), outline=(140, 148, 160, 255), width=2)
        ld.text((bw / 2, bh / 2), label, font=font, fill=(255, 255, 255, 255), anchor='mm')
        margin = round(vh * 0.05)
        self.cap = {'rgb': layer.convert('RGB'), 'alpha': layer.getchannel('A'), 'at': (vw - bw - margin, vh - bh - depth - margin),
                    'total': max(1, int(seconds * self.fps)), 'shown': 0}

    def key_frames(self, step, base):
        """The keys a demo's step pressed as its gesture, as a key cap shown a moment before the page they led to
        (the cap stays a moment over it)."""
        if not (self.demo and step.get('pressed')):
            return
        self.start_cap(key_label(step['pressed']))
        yield from self.hold_frames(0.3, base, self.cur)

    def zoom_frames(self, lit, at, z0, z1, caption):
        n = int(self.zoom_time * self.fps)
        for i in range(n):
            z = z0 + (z1 - z0) * ease((i + 1) / n)
            yield self.frame(lit, at, zoom=z, anchor=at, caption=caption)

    def press_frames(self, lit, at, press, target, step, index, caption, zoom, anchor):
        """The ripple of a press, the pointer's tip on its centre; the frame a fifth of the way in is
        kept for the audit, with the same frame bare and with the mark only, to measure both."""
        button = press.get('button', 'left')
        n = int((self.pulse_s if (press.get('count') or 1) == 1 else self.pulse_s * 0.6) * self.fps)
        peak = max(0, n // 5)
        # a press only the video shows (a native select's drawn list) is on the element it names, whatever is lit
        if press.get('shown') and press.get('box'):
            target = press['box']
        for i in range(n):
            t = (i + 1) / n
            image = self.frame(lit, at, press=(at, t, button), zoom=zoom, anchor=anchor, caption=caption)
            if i == peak:
                self.presses.append({
                    'frame': self.frame_no - 1, 'step': index, 'caption': step['caption'], 'x': at[0], 'y': at[1],
                    'button': button, 'count': press.get('count') or 1, 'on': press.get('on'), 'target': target,
                    'inside': None if target is None else inside(target, at[0], at[1]),
                    'video': self.to_video(at, zoom, anchor),
                    # measured without a demo's note over it
                    'images': (self.plain,
                               self.compose(lit, None, press=(at, t, button), zoom=zoom, anchor=anchor, caption=caption),
                               self.compose(lit, None, zoom=zoom, anchor=anchor, caption=caption))})
            yield image

    def drag_frames(self, lit, target, path, button, caption):
        """The pointer along a drag, the button held: each picture taken on the way is cut in as the
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
            n = self.frames_for(math.hypot(dest[0] - start[0], dest[1] - start[1]), int((0.2 if self.demo else 0.3) * self.fps))
            for i in range(n):
                t = ease((i + 1) / n)
                pos = (start[0] + (dest[0] - start[0]) * t, start[1] + (dest[1] - start[1]) * t)
                yield self.frame(drawn if drawn is not None and i == n - 1 else lit, pos, held=button, caption=caption)
            if drawn is not None:
                lit = drawn
        yield from self.hold_frames(0.3 if self.demo else 0.5, lit, self.cur, caption=caption)
        return lit

    def step_frames(self, step, index, total):
        fps = self.fps
        after = self._image(step['after'])
        caption = self.caption_layer(step, index, total)
        legs = step.get('legs', [])
        still_path = os.path.join(self.folder, f'step-{index:02d}.png')
        if self.demo:
            self.end_stale_note(step)
        # a line of text: a demo shows it over the page while the steps after it go on, a guide as its caption
        if step['kind'] == 'note':
            self.compose(after, self.cur, caption=caption).save(still_path)
            if self.demo:
                yield from self.finish_note()
                self.start_note(step['caption'], self.ahead(step), self.before(step), step.get('noteAt'),
                                step.get('view'))
                page = self._image(step['settled']) if step.get('settled') else after
                # a tutorial dissolves into a page the set-up after the note changed, rather than cutting to it
                if self.tutorial and self.base_now is not None and self.change(self.base_now, page) > 0.01:
                    n = int(0.35 * self.fps)
                    start = self.base_now
                    for i in range(n):
                        yield self.frame(Image.blend(start, page, ease((i + 1) / n)), self.cur)
                yield from self.hold_frames(0.8, page, self.cur)
            else:
                yield from self.hold_frames(self.reading(step['caption']), after, self.cur, caption=caption)
            return
        if self.demo and step.get('captioned'):
            yield from self.finish_note()
            self.start_note(step['caption'], self.ahead(step), self.before(step), step.get('noteAt'),
                            step.get('view'))
        if step['kind'] == 'check' and self.demo:
            self.compose(after, self.cur).save(still_path)
            yield from self.hold_frames(1.0, after, self.cur)
            return
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
        # a step with nothing to point at (a table opened through the API, a key pressed) shows the page it led to
        if not any(leg.get('target') or leg.get('pointer') for leg in legs):
            self.compose(after, self.cur, caption=caption).save(still_path)
            if self.demo:
                yield from self.key_frames(step, self._image(step['before']))
                yield from self.demo_result(step, self._image(step['before']), after, self._image(step['before']))
            else:
                yield from self.hold_frames(self.hold, after, self.cur, caption=caption)
            return
        # the stops of the step in order — one, or a menu's group, each item, the leaf: the pointer
        # goes where the page says it went on each, over the page as it was then, and every press
        # ripples under its tip. A stop's picture is the page as the pointer set off for it; what the
        # pointer's arrival led to (a group opened, an item lit) is the next stop's picture, cut to as
        # the pointer gets there
        still = None
        zoom, zoom_at = 1.0, None
        lit, target = None, None
        # the picture on screen as a stop begins (the one its arrival cut to), and the stop before it
        shown, came_from = None, None
        for n_leg, leg in enumerate(legs):
            shot = self._image(leg['shot'])
            target = leg.get('target')
            lit = self.spotlight(shot, target) if target else shot
            kind = leg.get('kind')
            events = list(leg.get('pointer', []))
            # a stop of a menu the pointer made no move for (a group already open) is gone to all the same: the
            # pointer travels to the item, which does not light up wherever the pointer happens to be
            if kind and target and not any(e['type'] in ('move', 'down') for e in events):
                events.insert(0, {'type': 'move', 'x': target['x'] + target['width'] / 2, 'y': target['y'] + target['height'] / 2})
            lit_up = False
            still = None
            i = 0
            while i < len(events):
                e = events[i]
                at = (e['x'], e['y'])
                # a picture taken while the step went on (text typed so far, a wheel turn): cut in, the pointer still
                if e['type'] == 'still':
                    shot = self._image(e['shot'])
                    lit = self.spotlight(shot, target) if target else shot
                    lit_up = True
                    yield from self.hold_frames(self.beat_s, lit, self.cur, zoom=zoom, anchor=zoom_at, caption=caption)
                    i += 1
                    continue
                if e['type'] == 'up' or (e['type'] == 'move' and e.get('held')):
                    i += 1
                    continue
                if zoom > 1 and math.hypot(at[0] - zoom_at[0], at[1] - zoom_at[1]) > 2:
                    yield from self.zoom_frames(lit, zoom_at, zoom, 1.0, caption)
                    zoom = 1.0
                via = self.elbow(came_from, target, at) if kind and not lit_up else None
                yield from self.travel_frames(shown if shown is not None and not lit_up else shot, lit, at, caption,
                                              fade=not lit_up, via=via)
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
                # it leaves the zoomed picture); a press rests under the pointer before it lands — a
                # demo names the keys held through it (Shift) on a key cap
                if is_small(target) and zoom == 1.0 and self.zoom > 1 and not path:
                    yield from self.zoom_frames(lit, at, 1.0, self.zoom, caption)
                    zoom, zoom_at = self.zoom, at
                if self.demo and e.get('mods'):
                    self.start_cap(key_label(['+'.join(e['mods'])]), self.rest_s + self.pulse_s + 0.2)
                if (e.get('count') or 1) == 1:
                    yield from self.hold_frames(self.menu_rest_s if kind == 'menu' else self.rest_s, lit, at, zoom=zoom,
                                                anchor=zoom_at, caption=caption)
                yield from self.press_frames(lit, at, e, target, step, index, caption, zoom, zoom_at)
                still = self.compose(lit, at, press=(at, 0.3, e.get('button', 'left')), caption=caption)
                if path:
                    lit = yield from self.drag_frames(lit, target, path, e.get('button', 'left'), caption)
                    still = self.compose(lit, self.cur, held=e.get('button', 'left'), caption=caption)
                i = j
            if not lit_up and target:
                yield from self.light_frames(shown if shown is not None else shot, lit, caption)
            if n_leg < len(legs) - 1:
                # the arrival: the next stop's picture, the stop just reached still lit, held — a menu's stop and a
                # list's long enough to read what opened
                nxt = self._image(legs[n_leg + 1]['shot'])
                shown = self.spotlight(nxt, target) if target else nxt
                yield from self.arrive_frames(lit, shown, zoom, zoom_at, caption)
                hold = self.menu_s if kind == 'menu' else self.list_s if kind == 'list' else self.leg_s
                yield from self.hold_frames(hold, shown, self.cur, zoom=zoom, anchor=zoom_at, caption=caption)
                if zoom > 1:
                    yield from self.zoom_frames(shown, zoom_at, zoom, 1.0, caption)
                    zoom = 1.0
            came_from = target
        # the picture of the step for steps.md: the last stop lit, and where the pointer acted on it
        (still if still is not None else self.compose(lit, self.cur, caption=caption)).save(still_path)
        yield from self.hold_frames(self.cut_s, lit, self.cur, zoom=zoom, anchor=zoom_at, caption=caption)
        if self.demo:
            yield from self.key_frames(step, lit)
            # the rest is measured on the page as filmed: a tutorial's outline and dimming are no change of the page
            yield from self.demo_result(step, shot, after, lit)
            return
        # the page after the step: zooming out over the step's picture, then cut to the page after it,
        # whose light fades — a cross-fade of two pictures of a changed page shows both at once
        if zoom > 1:
            yield from self.zoom_frames(lit, zoom_at, zoom, 1.0, caption)
        lit_after = self.spotlight(after, target) if target else after
        n = int(0.35 * fps) if lit_after is not after else 0
        for i in range(n):
            yield self.frame(Image.blend(lit_after, after, ease((i + 1) / n)), self.cur, caption=caption)
        yield from self.hold_frames(self.hold, after, self.cur, caption=caption)

    def frames(self, end_hold):
        if not self.steps:
            return
        first = self._image(self.steps[0]['before'])
        self.marks.append((self.frame_no, 'title'))
        if self.help:
            yield from self.hold_frames(1.0, first, self.cur)
        elif self.demo:
            yield from self.demo_title_frames(first)
        else:
            yield from self.title_frames(first)
        total = len(self.steps)
        for i, step in enumerate(self.steps, 1):
            self.marks.append((self.frame_no, f"{i}. {step['caption']}"))
            yield from self.step_frames(step, i, total)
        self.marks.append((self.frame_no, 'end'))
        last = self.steps[-1]
        # a note still being read keeps the last picture up until it is
        rest = (self.note['total'] - self.note['shown']) / self.fps if self.note else 0
        yield from self.hold_frames(max(end_hold, rest), self._image(last.get('settled') or last['after']), self.cur)
        # a demo ends in the dark its next clip starts from
        if self.demo:
            self.marks.append((self.frame_no, 'fade'))
            dark = Image.new('RGB', self.video, NOTE_BOX)
            n = int(0.35 * self.fps)
            for i in range(n):
                self.frame_no += 1
                yield Image.blend(self.plain, dark, ease((i + 1) / n))

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

    def encode(self, ffmpeg, gif, gif_fps=12, colors=256, gif_width=None, end_hold=0.8):
        h264 = has_encoder(ffmpeg, 'libx264')
        video = os.path.join(self.folder, 'guide.mp4' if h264 else 'guide.webm')
        # a demo is watched on a large screen and cut into a longer video: its text is kept sharper
        quality = ['-preset', 'medium', '-crf', '18'] if self.demo else ['-preset', 'veryfast', '-crf', '22']
        codec = ['-c:v', 'libx264', *quality, '-pix_fmt', 'yuv420p', '-movflags', '+faststart'] \
            if h264 else ['-c:v', 'libvpx', '-b:v', '2M', '-pix_fmt', 'yuv420p']
        cmd = [ffmpeg, '-y', '-loglevel', 'error', '-f', 'rawvideo', '-pix_fmt', 'rgb24', '-s', f'{self.video[0]}x{self.video[1]}',
               '-r', str(self.fps), '-i', '-', *codec, video]
        # the GIF is made from the frames themselves, at its own rate, not from the compressed video
        frames_dir = os.path.join(self.folder, 'gif-frames')
        if gif:
            shutil.rmtree(frames_dir, ignore_errors=True)
            os.makedirs(frames_dir)
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
        count = kept = 0
        # the last frame of every part (the title, each step, the end), for a review of a demo at a glance
        ends, seen, last = [], 0, None
        try:
            for frame in self.frames(end_hold):
                if len(self.marks) > seen:
                    if seen > 0:
                        ends += [last] * (len(self.marks) - seen)
                    seen = len(self.marks)
                proc.stdin.write(frame.tobytes())
                if gif and count * gif_fps // self.fps >= kept:
                    frame.save(os.path.join(frames_dir, f'{kept:05d}.png'), compress_level=1)
                    kept += 1
                if self.demo:
                    last = frame.resize((480, 270), Image.BILINEAR)
                count += 1
        finally:
            proc.stdin.close()
            code = proc.wait()
        if code != 0:
            raise RuntimeError(f'ffmpeg failed encoding {video}')
        outputs = [video]
        if self.demo and last is not None:
            ends.append(last)
            outputs += self.write_review(ends, count)
        if gif and self.steps:
            out = os.path.join(self.folder, 'guide.gif')
            width = min(gif_width or self.video[0], self.video[0]) & ~1
            frames = ['-framerate', str(gif_fps), '-i', os.path.join(frames_dir, '%05d.png')]
            scale = f'scale={width}:-2:flags=lanczos'
            palette = os.path.join(self.folder, 'palette.png')
            # no dithering: it defeats GIF's run-length compression on flat UI colours
            subprocess.run([ffmpeg, '-y', '-loglevel', 'error', *frames, '-vf',
                            f'{scale},palettegen=max_colors={colors}:stats_mode=diff', palette], check=True)
            with_press_colours(palette)
            subprocess.run([ffmpeg, '-y', '-loglevel', 'error', *frames, '-i', palette, '-lavfi',
                            f'[0:v]{scale}[v];[v][1:v]paletteuse=dither=none:diff_mode=rectangle', out], check=True)
            os.remove(palette)
            # the thumbnail is the result, without the strip, at the GIF's size
            result = self._image(self.steps[-1]['after'])
            thumb = os.path.join(self.folder, 'guide-thumb.png')
            result.resize((width, round(result.height * width / result.width)), Image.LANCZOS).save(thumb)
            sheet = os.path.join(self.folder, 'guide-sheet.png')
            self.contact_sheet(frames_dir, kept, sheet)
            shutil.rmtree(frames_dir, ignore_errors=True)
            self.gif = {'path': out, 'width': width, 'frames': kept, 'seconds': round(kept / gif_fps, 1),
                        'kb': round(os.path.getsize(out) / 1024)}
            outputs += [out, thumb, sheet]
        return outputs, count

    def write_review(self, ends, count):
        """timeline.json: when each part of the demo starts and how long it lasts; review-sheet.png: the
        last frame of each, with its start and its step."""
        parts = []
        for i, (start, label) in enumerate(self.marks):
            end = self.marks[i + 1][0] if i + 1 < len(self.marks) else count
            parts.append({'part': label, 'start': round(start / self.fps, 2), 'seconds': round((end - start) / self.fps, 2)})
        with open(os.path.join(self.folder, 'timeline.json'), 'w', encoding='utf-8') as f:
            json.dump({'seconds': round(count / self.fps, 2), 'parts': parts}, f, indent=2)
        cols, tw, th, label_h = 4, 480, 270, 40
        shown = [(thumb, part) for thumb, part in zip(ends, parts) if part['part'] != 'fade']
        sheet = Image.new('RGB', (cols * tw, math.ceil(len(shown) / cols) * (th + label_h)), (255, 255, 255))
        d = ImageDraw.Draw(sheet)
        font = load_font(13)
        for i, (thumb, part) in enumerate(shown):
            x, y = (i % cols) * tw, (i // cols) * (th + label_h)
            sheet.paste(thumb, (x, y + label_h))
            d.text((x + 6, y + 4), f"{part['start']:.1f} s +{part['seconds']:.1f} s", font=font, fill=(160, 0, 0))
            d.text((x + 6, y + 21), part['part'][:70], font=font, fill=(30, 30, 30))
        path = os.path.join(self.folder, 'review-sheet.png')
        sheet.save(path)
        return [os.path.join(self.folder, 'timeline.json'), path]

    def contact_sheet(self, frames_dir, count, path, n=12, cols=4, tile_w=400):
        """Twelve frames of the GIF, evenly spaced, for a review at a glance."""
        picks = sorted({round(i * (count - 1) / (n - 1)) for i in range(n)}) if count > 1 else [0]
        tiles = [Image.open(os.path.join(frames_dir, f'{i:05d}.png')).convert('RGB') for i in picks]
        tile_h = round(tiles[0].height * tile_w / tiles[0].width)
        sheet = Image.new('RGB', (cols * tile_w, math.ceil(len(tiles) / cols) * tile_h), (255, 255, 255))
        for k, tile in enumerate(tiles):
            sheet.paste(tile.resize((tile_w, tile_h), Image.LANCZOS), ((k % cols) * tile_w, (k // cols) * tile_h))
        sheet.save(path)


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


def reel(ffmpeg, out, folders):
    """Demo clips cut together in the order given: every demo is 1920×1080 at the same rate, so the cut is a
    stream copy. Each clip is copied beside the reel under its scenario's name, and `<out>.md` lists them."""
    clips = []
    for folder in folders:
        video, timeline = os.path.join(folder, 'guide.mp4'), os.path.join(folder, 'timeline.json')
        if not (os.path.exists(video) and os.path.exists(timeline)):
            print(f'guide-render: {folder} holds no demo clip (guide.mp4 and timeline.json): film it with @demo or --demo',
                  file=sys.stderr)
            return 2
        with open(os.path.join(folder, 'steps.json'), encoding='utf-8') as f:
            manifest = json.load(f)
        with open(timeline, encoding='utf-8') as f:
            clips.append((video, manifest, json.load(f)['seconds']))
    target = os.path.dirname(os.path.abspath(out))
    os.makedirs(target, exist_ok=True)
    listing, rows, start = [], [], 0.0
    for i, (video, manifest, seconds) in enumerate(clips, 1):
        title = manifest['scenario']
        # Windows takes no \ / : * ? " < > | in a file name, nor a trailing dot or space
        name = ''.join('-' if c in '\\/:*?"<>|' else c for c in title).rstrip('. ')
        shutil.copyfile(video, os.path.join(target, f'{i:02d} - {name}.mp4'))
        # the concat list quotes with ', so a ' in a path ("A sequence's structure") is closed, escaped and reopened
        listing.append("file '" + os.path.abspath(video).replace('\\', '/').replace("'", "'\\''") + "'")
        rows.append(f"| {i:02d} | {int(start // 60)}:{start % 60:04.1f} | {seconds:.1f} s | {title} | "
                    f"{manifest.get('description', '')} | {manifest['feature']} |")
        start += seconds
    stem = os.path.splitext(out)[0]
    with open(stem + '.concat.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(listing) + '\n')
    try:
        cut = subprocess.run([ffmpeg, '-y', '-loglevel', 'error', '-f', 'concat', '-safe', '0', '-i', stem + '.concat.txt',
                              '-c', 'copy', '-movflags', '+faststart', out])
    finally:
        os.remove(stem + '.concat.txt')
    if cut.returncode != 0:
        return cut.returncode
    with open(stem + '.md', 'w', encoding='utf-8') as f:
        f.write(f"# {os.path.basename(stem)}\n\n{len(clips)} clips, {start / 60:.1f} min; `{os.path.basename(out)}` is all of "
                f"them in this order, and each is beside it as `NN - <title>.mp4`.\n\n"
                "| # | Starts | Length | Title | Card text | Feature |\n|---|---|---|---|---|---|\n" + '\n'.join(rows) + '\n')
    print(f'{len(clips)} clips, {start:.1f} s: {out}, {stem}.md')
    return 0


def selftest(folder):
    os.makedirs(folder, exist_ok=True)
    w, h = 640, 400
    names = ['01-before', '01-after', '02-before', '02-after', '03-before', '03-drag1', '03-drag2', '03-after',
             '04-before', '04-after', '05-before', '05-hop2', '05-after', '06-before', '06-after',
             '08-before', '08-hop2', '08-hop3', '08-after', '09-before', '09-after', '10-before', '10-hop2', '10-hop3', '10-after']
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
        # a menu: a dropdown of three groups, the middle one opened into a flyout beside it, then its item lit
        if name.startswith('08-'):
            d.rectangle((300, 76, 420, 148), fill=(255, 255, 255), outline=(150, 150, 150))
            for k, label in enumerate(['Above', 'Group', 'Below']):
                if name in ('08-hop2', '08-hop3') and k == 1:
                    d.rectangle((301, 100, 419, 124), fill=(214, 228, 245))
                d.text((310, 82 + 24 * k), label, fill=(40, 40, 40))
            if name in ('08-hop2', '08-hop3'):
                d.rectangle((422, 100, 566, 196), fill=(255, 255, 255), outline=(150, 150, 150))
                if name == '08-hop3':
                    d.rectangle((423, 148, 565, 172), fill=(214, 228, 245))
                for k, label in enumerate(['One', 'Two', 'Item', 'Four']):
                    d.text((432, 106 + 24 * k), label, fill=(40, 40, 40))
        # a native select, its list drawn under it with the current option lit, then the option lit
        if name.startswith('10-'):
            d.rectangle((60, 200, 180, 222), fill=(255, 255, 255), outline=(150, 150, 150))
            d.text((66, 205), 'first' if name != '10-after' else 'third', fill=(40, 40, 40))
            if name in ('10-hop2', '10-hop3'):
                d.rectangle((60, 224, 180, 296), fill=(255, 255, 255), outline=(150, 150, 150))
                d.rectangle((61, 225 if name == '10-hop2' else 273, 179, 247 if name == '10-hop2' else 295), fill=(47, 111, 222))
                for k, label in enumerate(['first', 'second', 'third']):
                    d.text((70, 229 + 24 * k), label, fill=(40, 40, 40))
        d.text((300, 320), name, fill=(40, 40, 40))
        img.save(os.path.join(folder, name + '.png'))

    def press(x, y, button='left', count=1):
        return [{'type': 'down', 'x': x, 'y': y, 'button': button, 'count': count, 'on': 'selftest'},
                {'type': 'up', 'x': x, 'y': y, 'button': button}]

    def step(i, text, caption, kind, legs, keys=()):
        return {'index': i, 'line': i, 'keyword': 'Then' if kind == 'check' else 'When', 'text': text, 'caption': caption,
                'kind': kind, 'before': f'{i:02d}-before.png', 'after': f'{i:02d}-after.png', 'legs': legs,
                'keys': list(keys), 'typed': '', 'ms': 10}
    button = {'x': 40, 'y': 40, 'width': 80, 'height': 24}
    icon = {'x': 600, 'y': 8, 'width': 20, 'height': 20}
    group = {'x': 300, 'y': 100, 'width': 120, 'height': 24}
    item = {'x': 422, 'y': 148, 'width': 144, 'height': 24}
    select = {'x': 60, 'y': 200, 'width': 120, 'height': 22}
    option = {'x': 60, 'y': 272, 'width': 120, 'height': 24}

    def shown(box):
        x, y = box['x'] + box['width'] / 2, box['y'] + box['height'] / 2
        return [{'type': 'move', 'x': x, 'y': y}, {'type': 'down', 'x': x, 'y': y, 'button': 'left', 'count': 1, 'box': box,
                                                   'on': 'option', 'shown': True}, {'type': 'up', 'x': x, 'y': y, 'button': 'left'}]
    manifest = {'feature': 'Self test', 'scenario': 'Two synthetic steps', 'description': '', 'tags': [],
                'viewport': {'width': w, 'height': h}, 'steps': [
                    # the press is off the button's centre, where the mark and the tip must be; it opens a dialog
                    dict(step(1, 'user clicks on the blue button', 'Click on the blue button', 'action',
                              [{'shot': '01-before.png', 'target': button, 'pointer': [{'type': 'move', 'x': 52, 'y': 47}] + press(52, 47)}]),
                         dialog={'inputs': 3}),
                    step(2, 'the page should be blue', 'The page is blue', 'check', [], ['Enter']),
                    step(3, 'user drags across the page', 'Drag across the page', 'action',
                         [{'shot': '03-before.png', 'target': {'x': 180, 'y': 80, 'width': 300, 'height': 200}, 'pointer': [
                             {'type': 'move', 'x': 200, 'y': 100}, {'type': 'down', 'x': 200, 'y': 100, 'button': 'left', 'count': 1},
                             {'type': 'move', 'x': 230, 'y': 120, 'held': True},
                             {'type': 'move', 'x': 260, 'y': 140, 'held': True, 'shot': '03-drag1.png'},
                             {'type': 'move', 'x': 320, 'y': 180, 'held': True, 'shot': '03-drag2.png'},
                             {'type': 'up', 'x': 320, 'y': 180, 'button': 'left'}]}]),
                    # a right-click on an icon at the page's edge, Shift held: zoomed into, green
                    step(4, 'user right-clicks on the grey icon', 'Right-click on the grey icon', 'action',
                         [{'shot': '04-before.png', 'target': icon,
                           'pointer': [dict(p, mods=['Shift']) if p['type'] == 'down' else p for p in press(611, 19, 'right')]}]),
                    # a menu: the group hovered, then the item clicked on the page the hover opened
                    step(5, 'user picks "Item" from the menu', 'Pick "Item" from the menu', 'action',
                         [{'shot': '05-before.png', 'target': button, 'pointer': [{'type': 'move', 'x': 80, 'y': 52}]},
                          {'shot': '05-hop2.png', 'target': {'x': 300, 'y': 180, 'width': 160, 'height': 30},
                           'pointer': [{'type': 'move', 'x': 330, 'y': 195}] + press(330, 195)}]),
                    step(6, 'user double-clicks on the blue button', 'Double-click on the blue button', 'action',
                         [{'shot': '06-before.png', 'target': button, 'pointer': press(100, 58) + press(100, 58, count=2)}]),
                    # a line of text over the page, with nothing settled after it
                    dict(step(7, 'the video says "Blue is the button"', 'Blue is the button', 'note', []),
                         before='06-after.png', after='06-after.png', captioned=True),
                    # a menu walk: the group (its flyout opens as the pointer arrives), into the flyout along the
                    # group's row, the item lit before the press
                    step(8, 'user picks "Group > Item" from the menu', 'Pick "Group > Item" from the menu', 'action',
                         [{'shot': '08-before.png', 'target': group, 'kind': 'menu', 'pointer': [{'type': 'move', 'x': 360, 'y': 112}]},
                          {'shot': '08-hop2.png', 'target': item, 'kind': 'menu', 'pointer': [{'type': 'move', 'x': 494, 'y': 160}]},
                          {'shot': '08-hop3.png', 'target': item, 'kind': 'menu', 'pointer': press(494, 160)}]),
                    # a key the step pressed, nothing to point at
                    dict(step(9, 'user presses Escape', 'Press Escape', 'action', []), pressed=['Escape']),
                    # a native select: pressed, its drawn list up with the current option lit, the option lit as the
                    # pointer reaches it and pressed on film
                    step(10, 'user selects "third" in the select', 'Select "third" in the select', 'action',
                         [{'shot': '10-before.png', 'target': select, 'kind': 'list', 'pointer': [dict(p, on='select') for p in shown(select)]},
                          {'shot': '10-hop2.png', 'target': option, 'kind': 'menu',
                           'pointer': [{'type': 'move', 'x': 120, 'y': 284, 'kept': True}]},
                          {'shot': '10-hop3.png', 'target': option, 'kind': 'menu', 'pointer': shown(option)}])]}
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
    help = flag('--help')
    tutorial = flag('--tutorial') and not help
    demo = (flag('--demo') or tutorial) and not help
    fps = int(opt('--fps', 30))
    hold = float(opt('--hold', 0.9 if help else 1.4))
    gif_fps = int(opt('--gif-fps', 12))
    colors = int(opt('--colors', 256))
    gif_width = opt('--gif-width')
    gif_width = int(gif_width) if gif_width else None if help else 880
    end_hold = float(opt('--end-hold', 2.0 if help or demo else 0.8))
    travel = float(opt('--travel', 0.5 if demo and not tutorial else 0.7))
    zoom = float(opt('--zoom', 1.8))
    zoom_time = float(opt('--zoom-time', 0.8))
    ffmpeg = find_ffmpeg(opt('--ffmpeg'))
    reel_out = opt('--reel')
    test_dir = opt('--selftest')
    if test_dir:
        selftest(test_dir)
        args.append(test_dir)
        gif = True
    if len(args) != 1 and not (reel_out and args):
        print(__doc__, file=sys.stderr)
        return 2
    if not ffmpeg:
        print('guide-render: no ffmpeg — py -m pip install imageio-ffmpeg, or set FFMPEG', file=sys.stderr)
        return 3
    if reel_out:
        return reel(ffmpeg, reel_out, args)
    renderer = Renderer(args[0], fps, hold, travel, zoom, zoom_time, help, demo, tutorial)
    renderer.write_markdown()
    outputs, count = renderer.encode(ffmpeg, gif, gif_fps, colors, gif_width, end_hold)
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
    gif_info = getattr(renderer, 'gif', None)
    if gif_info and not quiet:
        print(f"gif: {gif_info['width']} px wide, {gif_info['frames']} frames, {gif_info['seconds']} s, {gif_info['kb']} KB")
    # a GIF for a help page shows a small window at twice its density; fewer pixels would be upscaled by the page
    soft = help and renderer.scale < 2
    if soft:
        print(f'guide-render: the pictures are {renderer.scale}x the window, not 2x: the help GIF would be soft', file=sys.stderr)
    return 4 if off or soft else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
