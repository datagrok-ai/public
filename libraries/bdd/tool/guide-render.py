"""Renders a guide manifest (`steps.json` + the PNGs a `BDD_GUIDE` run wrote) into guide.mp4,
step-NN.png (the annotated picture of every step), steps.md, and with --gif also guide.gif and
guide-thumb.png. Pillow does the drawing, ffmpeg the encoding (FFMPEG env, PATH, the
imageio-ffmpeg package, or Playwright's own build, which can only write WebM).

    py tool/guide-render.py <scenario dir> [--gif] [--fps 30] [--hold 1.4] [--travel 0.7]
        [--zoom 1.8] [--zoom-time 0.8] [--ffmpeg <path>] [--quiet]
    py tool/guide-render.py --selftest <dir>      renders a synthetic two-step guide there
"""
import glob
import json
import math
import os
import shutil
import subprocess
import sys

try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:
    print('guide-render: Pillow is missing — py -m pip install pillow', file=sys.stderr)
    sys.exit(3)

ACCENT = (255, 90, 31)
CHECK = (46, 160, 67)
BAR = (20, 24, 32, 205)
SMALL_TARGET_PX = 28   # a target no wider and no taller than this (an icon, a checkbox) is zoomed into; larger ones are lit and clicked in place


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
    pts = [(0, 0), (0, 16), (4, 12.5), (7, 19.5), (10, 18), (7, 11.5), (12.5, 11.5)]
    pts = [(x * scale, y * scale) for x, y in pts]
    w = int(14 * scale) + 4
    h = int(21 * scale) + 4
    sprite = Image.new('RGBA', (w, h), (0, 0, 0, 0))
    d = ImageDraw.Draw(sprite)
    d.polygon([(x + 2.5, y + 2.5) for x, y in pts], fill=(0, 0, 0, 90))
    d.polygon([(x + 1, y + 1) for x, y in pts], fill=(255, 255, 255, 255), outline=(0, 0, 0, 255))
    return sprite


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
        self.font = load_font(max(14, self.h // 34))
        self.small = load_font(max(12, self.h // 46))
        self.big = load_font(max(20, self.h // 18))
        self.cursor = cursor_sprite(max(1.0, self.h / 700))
        self.bar_h = int(self.h * 0.085) & ~1
        # the caption sits above the page: nothing of the page is covered, and a player's timeline
        # lands on the bottom edge
        self.out_h = self.h + self.bar_h
        self.cur = (self.w / 2, self.h / 2)

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
        d.text((x, self.bar_h / 2), text, font=self.font, fill=(255, 255, 255, 255), anchor='lm')
        x += d.textlength(text, font=self.font) + pad
        for key in step.get('keys', [])[:4]:
            label = key.replace('ControlOrMeta', 'Ctrl')
            tw = d.textlength(label, font=self.small)
            d.rounded_rectangle((x, self.bar_h * 0.22, x + tw + 14, self.bar_h * 0.78), radius=5,
                                fill=(60, 66, 78, 255), outline=(120, 128, 140, 255))
            d.text((x + 7, self.bar_h / 2), label, font=self.small, fill=(255, 255, 255, 255), anchor='lm')
            x += tw + 22
        return layer

    def compose(self, base, cursor=None, ripple=None, zoom=1.0, anchor=None, caption=None):
        frame = base.copy()
        if cursor:
            frame.paste(self.cursor, (int(cursor[0]) - 1, int(cursor[1]) - 1), self.cursor)
        if ripple is not None and cursor:
            t, r = ripple, 8 + 34 * ripple
            overlay = Image.new('RGBA', frame.size, (0, 0, 0, 0))
            ImageDraw.Draw(overlay).ellipse((cursor[0] - r, cursor[1] - r, cursor[0] + r, cursor[1] + r),
                                            outline=ACCENT + (int(255 * (1 - t)),), width=4)
            frame = Image.alpha_composite(frame.convert('RGBA'), overlay).convert('RGB')
        if zoom > 1.001 and anchor:
            cw, ch = self.w / zoom, self.h / zoom
            x0 = min(max(0, anchor[0] - cw / 2), self.w - cw)
            y0 = min(max(0, anchor[1] - ch / 2), self.h - ch)
            frame = frame.crop((int(x0), int(y0), int(x0 + cw), int(y0 + ch))).resize((self.w, self.h), Image.BILINEAR)
        canvas = Image.new('RGB', (self.w, self.out_h), BAR[:3])
        canvas.paste(frame, (0, self.bar_h))
        if caption is not None:
            canvas.paste(caption, (0, 0), caption)
        return canvas

    def title_frames(self, base):
        layer = Image.new('RGBA', (self.w, self.h), (10, 14, 20, 150))
        frame = Image.alpha_composite(base.convert('RGBA'), layer).convert('RGB')
        d = ImageDraw.Draw(frame)
        d.text((self.w / 2, self.h * 0.44), self.manifest['scenario'], font=self.big, fill=(255, 255, 255), anchor='mm')
        d.text((self.w / 2, self.h * 0.44 + self.h * 0.09), self.manifest['feature'], font=self.font,
               fill=(200, 206, 214), anchor='mm')
        frame = self.compose(frame)
        for _ in range(int(1.8 * self.fps)):
            yield frame

    def travel_frames(self, base, lit, dest, caption):
        """The pointer eased from where it rests to dest, the target lighting up over the second half."""
        start = self.cur
        n = max(1, int(self.travel * self.fps))
        for i in range(n):
            t = ease((i + 1) / n)
            pos = (start[0] + (dest[0] - start[0]) * t, start[1] + (dest[1] - start[1]) * t)
            mix = ease((t - 0.5) / 0.5)
            yield self.compose(Image.blend(base, lit, mix) if mix > 0 else base, pos, caption=caption)
        self.cur = dest

    def step_frames(self, step, index, total):
        fps = self.fps
        before = self._image(step['before'])
        after = self._image(step['after'])
        target = step.get('target')
        caption = self.caption_layer(step, index, total)
        clicks = step.get('clicks', [])
        pointer = step.get('pointer', [])
        if step['kind'] == 'check':
            base = self.spotlight(after, target, CHECK) if target else after
            still = self.compose(base, self.cur, caption=caption)
            still.save(os.path.join(self.folder, f'step-{index:02d}.png'))
            for _ in range(int(self.hold * 0.9 * fps)):
                yield still
            return
        # a path walked inside the step (a menu: the group, then each item): the pointer travels to
        # every stop over the page as it was then, the last stop being the step's own target
        legs = [(self._image(h['shot']), h['target']) for h in step.get('hops', [])] or [(before, target)]

        def hit(stop, c):
            return stop is not None and stop['x'] <= c['x'] <= stop['x'] + stop['width'] and stop['y'] <= c['y'] <= stop['y'] + stop['height']

        # the click ripples where the walk clicked: the selector that opened a picker, the leaf of a
        # menu, the row taken — a click no stop accounts for rides on the last one
        clicked = [any(hit(stop, c) for c in clicks) for _, stop in legs]
        clicks_last = bool(clicks) and (clicked[-1] or not any(clicked))
        for (shot, stop), here in zip(legs[:-1], clicked[:-1]):
            dest = (stop['x'] + stop['width'] / 2, stop['y'] + stop['height'] / 2)
            for c in clicks:
                if hit(stop, c):
                    dest = (c['x'], c['y'])
            lit = self.spotlight(shot, stop)
            yield from self.travel_frames(shot, lit, dest, caption)
            still = self.compose(lit, dest, caption=caption)
            for _ in range(int(0.5 * fps)):
                yield still
            if here:
                n = int(0.4 * fps)
                for i in range(n):
                    yield self.compose(lit, dest, ripple=(i + 1) / n, caption=caption)
        before, target = legs[-1]
        if target:
            dest = (target['x'] + target['width'] / 2, target['y'] + target['height'] / 2)
        elif clicks:
            dest = (clicks[-1]['x'], clicks[-1]['y'])
        elif pointer:
            dest = (pointer[-1]['x'], pointer[-1]['y'])
        else:
            dest = None
        if clicks and target:
            for c in clicks:
                if hit(target, c):
                    dest = (c['x'], c['y'])
        lit = self.spotlight(before, target) if target else before
        lit_after = self.spotlight(after, target) if target else after
        # a picture of the step for steps.md: the target lit, the pointer on it, no zoom; a step
        # with nothing to point at (a view switch through the API) shows the page it led to
        if dest is None:
            still = self.compose(after, self.cur, caption=caption)
            still.save(os.path.join(self.folder, f'step-{index:02d}.png'))
            for _ in range(int(self.hold * fps)):
                yield still
            return
        # the pages pictured while a button was held: a drag's still shows what it drew, at its end
        drags = [p for p in pointer if p.get('shot')]
        if drags:
            drawn = self._image(drags[-1]['shot'])
            still = self.compose(self.spotlight(drawn, target) if target else drawn, (drags[-1]['x'], drags[-1]['y']), caption=caption)
        else:
            still = self.compose(lit, dest, caption=caption)
        still.save(os.path.join(self.folder, f'step-{index:02d}.png'))
        small = target is not None and target['width'] <= SMALL_TARGET_PX and target['height'] <= SMALL_TARGET_PX
        zoom = self.zoom if small and self.zoom > 1 else 1.0
        yield from self.travel_frames(before, lit, dest, caption)
        # an icon-sized target is zoomed into, slowly enough to be read as a zoom; any other is lit
        # where it is — either way the target rests under the pointer before the click lands
        if zoom > 1:
            n = int(self.zoom_time * fps)
            for i in range(n):
                z = 1 + (zoom - 1) * ease((i + 1) / n)
                yield self.compose(lit, dest, zoom=z, anchor=dest, caption=caption)
        rest = self.compose(lit, dest, zoom=zoom, anchor=dest, caption=caption)
        for _ in range(int(0.4 * fps)):
            yield rest
        if clicks_last:
            n = int(0.4 * fps)
            for i in range(n):
                yield self.compose(lit, dest, ripple=(i + 1) / n, zoom=zoom, anchor=dest, caption=caption)
        # the drag: the pointer travels to each pictured point over the page as it was, the next
        # picture fading in as it arrives, so the box grows under the pointer; the step then ends
        # from the last picture, the pointer where the button was released
        for p in drags:
            drawn = self.spotlight(self._image(p['shot']), target) if target else self._image(p['shot'])
            to = (p['x'], p['y'])
            n = max(1, int(0.3 * fps))
            for i in range(n):
                t = ease((i + 1) / n)
                pos = (dest[0] + (to[0] - dest[0]) * t, dest[1] + (to[1] - dest[1]) * t)
                mix = max(0.0, (t - 0.6) / 0.4)
                yield self.compose(Image.blend(lit, drawn, mix) if mix > 0 else lit, pos, zoom=zoom, anchor=dest, caption=caption)
            lit, dest = drawn, to
        if drags:
            rest = self.compose(lit, dest, zoom=zoom, anchor=dest, caption=caption)
            for _ in range(int(0.5 * fps)):
                yield rest
        n = int((self.zoom_time if zoom > 1 else 0.35) * fps)
        for i in range(n):
            t = ease((i + 1) / n)
            z = zoom + (1 - zoom) * t
            base = Image.blend(lit, lit_after, min(1.0, t * 2)) if t < 0.5 else Image.blend(lit_after, after, (t - 0.5) * 2)
            yield self.compose(base, dest, zoom=z, anchor=dest, caption=caption)
        tail = (pointer[-1]['x'], pointer[-1]['y']) if pointer else dest
        if math.hypot(tail[0] - dest[0], tail[1] - dest[1]) > 30:
            n = int(0.45 * fps)
            for i in range(n):
                t = ease((i + 1) / n)
                pos = (dest[0] + (tail[0] - dest[0]) * t, dest[1] + (tail[1] - dest[1]) * t)
                yield self.compose(after, pos, caption=caption)
            self.cur = tail
        still = self.compose(after, self.cur, caption=caption)
        for _ in range(int(self.hold * fps)):
            yield still

    def frames(self):
        if not self.steps:
            return
        yield from self.title_frames(self._image(self.steps[0]['before']))
        total = len(self.steps)
        for i, step in enumerate(self.steps, 1):
            yield from self.step_frames(step, i, total)
        last = self.compose(self._image(self.steps[-1]['after']), self.cur)
        for _ in range(int(0.8 * self.fps)):
            yield last

    def write_markdown(self):
        lines = [f"# {self.manifest['scenario']}", '']
        if self.manifest.get('description'):
            lines += [self.manifest['description'], '']
        for i, step in enumerate(self.steps, 1):
            mark = '\u2713 ' if step['kind'] == 'check' else ''
            lines += [f"{i}. {mark}{step['caption']}", '', f"   ![Step {i}](step-{i:02d}.png)", '']
        with open(os.path.join(self.folder, 'steps.md'), 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))

    def encode(self, ffmpeg, gif):
        h264 = has_encoder(ffmpeg, 'libx264')
        video = os.path.join(self.folder, 'guide.mp4' if h264 else 'guide.webm')
        codec = ['-c:v', 'libx264', '-preset', 'veryfast', '-crf', '22', '-pix_fmt', 'yuv420p', '-movflags', '+faststart'] \
            if h264 else ['-c:v', 'libvpx', '-b:v', '2M', '-pix_fmt', 'yuv420p']
        cmd = [ffmpeg, '-y', '-loglevel', 'error', '-f', 'rawvideo', '-pix_fmt', 'rgb24', '-s', f'{self.w}x{self.out_h}',
               '-r', str(self.fps), '-i', '-', *codec, video]
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
        count = 0
        for frame in self.frames():
            proc.stdin.write(frame.tobytes())
            count += 1
        proc.stdin.close()
        if proc.wait() != 0:
            raise RuntimeError(f'ffmpeg failed encoding {video}')
        outputs = [video]
        if gif and self.steps:
            out = os.path.join(self.folder, 'guide.gif')
            # a docs GIF: ten frames a second at 880 px, no dithering (which defeats GIF's run-length
            # compression on flat UI colors) — a seven-step guide lands around a megabyte
            width = min(880, self.w) & ~1
            vf = (f'fps=10,scale={width}:-2:flags=lanczos,split[a][b];[a]palettegen=max_colors=96:stats_mode=diff[p];'
                  '[b][p]paletteuse=dither=none:diff_mode=rectangle')
            subprocess.run([ffmpeg, '-y', '-loglevel', 'error', '-i', video, '-vf', vf, out], check=True)
            thumb = os.path.join(self.folder, 'guide-thumb.png')
            shutil.copyfile(os.path.join(self.folder, 'step-01.png'), thumb)
            outputs += [out, thumb]
        return outputs, count


def selftest(folder):
    os.makedirs(folder, exist_ok=True)
    w, h = 640, 400
    for name, color in [('01-before.png', (240, 243, 247)), ('01-after.png', (224, 236, 250)),
                        ('02-before.png', (224, 236, 250)), ('02-after.png', (224, 236, 250)),
                        ('03-before.png', (224, 236, 250)), ('03-drag1.png', (224, 236, 250)),
                        ('03-drag2.png', (224, 236, 250)), ('03-after.png', (224, 236, 250))]:
        img = Image.new('RGB', (w, h), color)
        d = ImageDraw.Draw(img)
        d.rectangle((40, 40, 120, 64), fill=(32, 131, 213))
        if name.startswith('03-drag'):
            d.rectangle((200, 100, 200 + 60 * int(name[7]), 100 + 40 * int(name[7])), fill=(180, 180, 180))
        d.text((300, 200), name, fill=(40, 40, 40))
        img.save(os.path.join(folder, name))
    manifest = {'feature': 'Self test', 'scenario': 'Two synthetic steps', 'description': '', 'tags': [],
                'viewport': {'width': w, 'height': h}, 'steps': [
                    {'index': 1, 'line': 1, 'keyword': 'When', 'text': 'user clicks on the blue button', 'caption': 'Click on the blue button',
                     'kind': 'action', 'before': '01-before.png', 'after': '01-after.png',
                     'target': {'x': 40, 'y': 40, 'width': 80, 'height': 24}, 'pointer': [{'x': 80, 'y': 52}],
                     'clicks': [{'x': 80, 'y': 52, 'button': 'left'}], 'keys': [], 'typed': '', 'ms': 10},
                    {'index': 2, 'line': 2, 'keyword': 'Then', 'text': 'the page should be blue', 'caption': 'The page is blue',
                     'kind': 'check', 'before': '02-before.png', 'after': '02-after.png', 'pointer': [], 'clicks': [],
                     'keys': ['Enter'], 'typed': '', 'ms': 5},
                    {'index': 3, 'line': 3, 'keyword': 'When', 'text': 'user drags across the page', 'caption': 'Drag across the page',
                     'kind': 'action', 'before': '03-before.png', 'after': '03-after.png',
                     'target': {'x': 180, 'y': 80, 'width': 300, 'height': 200},
                     'pointer': [{'x': 200, 'y': 100}, {'x': 260, 'y': 140, 'shot': '03-drag1.png'}, {'x': 320, 'y': 180, 'shot': '03-drag2.png'}],
                     'clicks': [{'x': 200, 'y': 100, 'button': 'left'}], 'keys': [], 'typed': '', 'ms': 10}]}
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
    if not quiet:
        print(f'{count} frames, {count / fps:.1f} s: ' + ', '.join(outputs))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
