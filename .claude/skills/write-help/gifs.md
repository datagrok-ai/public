# Help GIFs

Part of the `write-help` skill. What a help GIF must look like, which tool records it, and how.

## What a help GIF shows

**Content**

1. **One action and its result**, planned before recording: open, set, see the outcome. Cut
   everything before the first meaningful action (start about 1 s before it) and after the
   result.
2. **Show the whole interaction**: set, change, and undo. A filter is applied and reset, a pin is
   added and removed, a slider goes there and back, legend categories are toggled on and off.
3. **Cause before effect.** A selection, filter, or highlight never appears before the gesture
   that causes it. Reorder the steps rather than explain in words.
4. **Meaningful data** where the effect is visible: demo datasets the reader can open (`demog`,
   Northwind, ChEMBL demo, package samples), or a subset of them. No placeholder columns, no test
   entity names in titles or dialogs, no personal pins or favorites.
5. **The last frame is held about 2 s** so the reader sees what changed.
6. **A step that looks broken on video is dropped** (a clipped popover, a half-visible selector,
   a dialog off the frame). Better to leave the action out than show it broken.

**Frame**

7. **Only what the GIF is about.** Keep the shell (menus, tabs). Close Browse, Toolbox, Context
   Panel, and Help unless a step needs one. Open a panel only for that step and close it again.
   No full-screen mode (`Alt+F`): it is recorded shrunk with black bars. Move a splitter or close
   the grid instead.
8. **A viewer GIF shows the whole viewer.** Show the grid only when the result is visible in the
   data (selection, filter, pinned rows), and keep the status bar with **Selected** and
   **Filtered** in the frame for such scenarios.
9. **Change settings from the viewer's own context menu** where one exists: the Context Panel
   squeezes the view and draws attention away. When a setting lives only in the Context Panel,
   open it with the viewer's gear, filter it with its search box, and widen it with the splitter
   if labels are cut.
10. **Readable at help width.** Help shows GIFs about 800 px wide. Record at the final size,
    800×500, as `help/CLAUDE.md` requires. The one exception: when a panel must stay open, record at
    up to 1000×625 and scale down to 800×500 (1.25×). For a small
    preview (a 400×300 gallery tile), make the viewer itself small, for example docked next to the
    grid, and crop to it, instead of shrinking a large one.
11. **No spinners, no black edges**, light theme, default zoom, no browser chrome.

**Pointer**

12. **Every click is marked**: yellow for the left button, green for the right. The mark stays
    while the button is held, moves with the pointer during a drag, and lingers about 250 ms after
    release. A double click shows two marks.
13. **Movements are continuous.** Each action starts where the previous one ended, the pointer
    travels to the target on screen, and the button is pressed only when it has arrived. The drawn
    pointer never lags behind the element it moves.
14. **Drags go in short steps** so the viewer redraws on the way. Mouse-wheel scrolling needs a
    frame per wheel step, or it reads as a jump.
15. **Clicks are real UI gestures**, not API calls that imitate them.

**Captions and files**

16. Step captions over the frame are allowed, in the reader's words: no internal names,
    coordinates, or IDs. Dimming around a highlighted element is allowed where it helps find it.
17. **Replacing a GIF is an edit.** Analyse the original frame by frame (`git show HEAD:<path>`
    and a contact sheet) and repeat its whole scenario, fixing only what was asked. Keep its size
    unless the size is the complaint.
18. File name, thumbnail, and alt text: SKILL.md, "Images and GIFs".
19. **Leave the stand as you found it.** Delete every entity the recording created and verify
    it. Never save shared settings: stop before **Save and apply**.

## Choosing the tool

| Tool | Use when |
|---|---|
| `grok-bdd guide` (the `bdd-answer` skill, `libraries/bdd`) | The scenario can be written as a BDD feature and its output meets the rules above. The feature also stays as a regression test. `--help-pages` films every `@help:<page dir>` feature and copies the GIF and thumbnail into the page's `img/` |
| `tools/rec-lib.mjs` (this skill) | The output of `grok-bdd guide` breaks a rule above: a crop to one viewer, frames without captions, panels per step, a held last frame, a thumbnail of the result |

`grok-bdd guide` films at 1920×1080 unless `BDD_GUIDE_VIEWPORT=<w>x<h>` is set (below 1920 px the
top menu folds its last groups, such as Chem and Bio, into "more"). The GIF is at most 880 px
wide, and the thumbnail is the first step (`step-01.png`), not the result. For a help GIF, set the
viewport to the final size (800×500) and check the thumbnail; switch to `rec-lib.mjs` when the
table above says so.

## rec-lib.mjs

**Dependencies** (not in this folder):

- Node.js 18 or later.
- Playwright (`@playwright/test`) and its Chromium. By default it is resolved from `libraries/bdd`
  of this repository: `npm install` there, then `npx playwright install chromium`. Set
  `REC_PLAYWRIGHT_FROM` to another `package.json` to use another copy.
- A full ffmpeg build (Playwright's writes WebM only): `FFMPEG`, else the binary of the Python
  package `imageio-ffmpeg`, else `ffmpeg` on `PATH`.
- The `grok` CLI (`datagrok-tools`) with the stand as a server alias in `~/.grok/config.yaml`,
  signed in with `grok login <alias>`. The recorder reads the alias's `url` and gets a session
  token from `grok s token`. Nothing secret is printed.

**Environment**: `REC_HOST` (config alias, `dev` by default), `REC_W` / `REC_H` (window, 800×500
by default), `REC_BROWSER_ARGS` (Chromium flags; the default enables the hardware GPU, which 3D
viewers need), `REC_GROK` (the grok command), `REC_OUT` (output folder of the examples).

**API**

- `film(name, setup, scene, opts)`:
  - `setup` is JavaScript run in the page, not recorded.
  - `opts.pre(page)` runs Playwright steps, not recorded.
  - `scene(page, mouse)` is the recorded part.
  - Options: `out` (path without extension; default `./out/<name>` in the current directory, keep
    it outside the repository), `crop` (`[x, y, w, h]` in CSS pixels, or `async (page) => box`),
    `size` (default 800×500; with a crop, at most 800 px wide in the crop's aspect ratio), `fps`
    (12), `colors` (palette size, 128; 64 keeps long GIFs small), `thumbAt` (fraction of the
    timeline for the thumbnail), `noTooltips`, `start` (pointer's starting point).
  - It writes `<out>.gif`, `<out>-thumb.png`, and `<out>-sheet.png` (12 frames for review). A
    final `waitForTimeout(2000)` in the scene holds the last frame. The browser closes even when the
    scene fails.
- `Mouse`: `moveTo`, `click` (`{button: 'right'}`), `dblclick`, `drag`, `clickEl`, `hoverEl`,
  `type`. Every action travels to its target first.
- `walkMenu(mouse, labels)` walks a context-menu path: hovers each group, clicks the leaf.
- `open()` and `cleanShell(page)` for read-only probes without recording.

Run a scene with `node <scene>.mjs` from the folder that holds `rec-lib.mjs`, or import it by
relative path as the examples do.

**Examples** (`tools/examples/`):

| Script | Shows |
|---|---|
| `rec-sc-reg.mjs` | Viewer only, context-menu path, column selector popup, legend toggles, the property panel opened by the gear for one step |
| `rec-sc-ma.mjs` | Slider there and back, switching an axis column, a native `<select>` with the keyboard |
| `rec-lc-small.mjs` | A 400×300 gallery preview: viewer docked next to the grid, cropped to the viewer |
| `rec-helm.mjs` | A full-screen editor addressed by `data-testid`, waiting for slow tabs |
| `rec-cw.mjs` | A settings page, stopping before **Save and apply** |
| `rec-dd.mjs` | Opening the Browse tree in `pre`, saving a project and reopening it from **Dashboards** |
| `dd-projects.mjs` | `node dd-projects.mjs "<project name>" [--delete]`: lists or deletes the current user's projects with that name. Run it after every take, failed ones too |

## Process

1. **Probe before filming.** A short read-only script opens the page, prints labels
   (`.d4-menu-item-label`, `[name=…]`, a dialog's `innerText`), and takes screenshots. Choose
   locators by `name` and `data-testid`, not by coordinates.
2. **Keep the recorded part short.** Loading data, opening the tree, and closing panels go to
   `setup` or `pre`.
3. **Review every take** by its contact sheet, then by single frames at GIF size
   (`ffmpeg -ss <t> -i <gif> -frames:v 1 frame.png`): nothing cut off, labels readable, cause
   before effect, last frame held.
4. **Clean up the stand** and verify.
5. Place the GIF and check the page as in SKILL.md, "Check and hand over".

## Pitfalls

- A dialog input may keep its last value: press Ctrl+A before typing.
- Double-clicking a query in Browse opens a preview. Use **Run** from its context menu to get a
  table view with the **Source** pane and **REFRESH**.
- Clicking the sidebar Browse icon while Browse is open closes it.
- Tree nodes load lazily: wait a few seconds after expanding, and scroll with the mouse wheel over
  the tree until the node is in the window.
- Native `<select>` lists are not captured by the screencast: click the field, then use the arrows
  and Enter.
- Balloons a step triggers that explain nothing: hide them with
  `.d4-balloon, .d4-balloon-container { display: none !important; }`.
- Column selector popups open above or below the selector depending on space: read the focused
  input's position and click relative to it.
- A hover-only element (a title bar icon) appears only while its parent is hovered: move to the
  parent first.
- To give a viewer the whole view, drag the grid's splitter to the edge, then close the grid.
  Closing the grid alone leaves the remaining panel at its old size.
- The line chart's range slider is drawn on the canvas along the top edge of the X axis and only
  while the pointer is over the axis. Hover the axis, wait for it to draw, then grab a handle a few
  pixels inside the axis ends.
- The sunburst reports no hit areas: compute points in polar coordinates from the viewer's box.
  Small windows put inner rings in the empty center, so aim at the outer ring.
- A button label in capitals may be lowercase in the DOM (`Save` shown as SAVE): match it
  case-insensitively.
