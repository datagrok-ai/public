/* The progress tracker of a Hit Design / PeptiHit campaign: the platform Tile Viewer over the campaign's table, one
   lane per stage (and one per other stage a row holds), in a modal dialog over the campaign view — pinned to it:
   hidden while another view is current, closed with it. A card dragged to another lane moves its row to that stage at
   once; closing the tracker, whichever way, keeps the moves and saves them. The card layout a user makes (Edit
   Form...) is kept with the campaign. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Dialog, RecordsInput, TextInput, computed, div, signal, statusLine, subscribe} from '@datagrok-libraries/u2';
import type {ViewHost} from '@datagrok-libraries/u2/dg';
import type {HitDesignApp} from '../hit-design-app';
import {TileCategoriesColName} from '../consts';

/** A lane as wide as its card, the lanes scrolling sideways when they do not fit: the platform's own card (no
 * layout kept with the campaign) needs 368 px, so five lanes fit at 1920; a kept layout reaches its rightmost element,
 * plus the card's frame. */
const LANE_WIDTH = 368;
const CARD_FRAME = 40;

/** The progress tracker over the campaign of `host`'s view, shown pinned to the view; answers its dialog. */
export function showProgressTracker(app: HitDesignApp, host: ViewHost): Dialog {
  const view = host.view as DG.TableView;
  const df = view.dataFrame;
  const dialog = new Dialog('Progress Tracker', {name: 'progress-tracker', closeOnly: true});
  const lanes = div([], 'hit-design-tracker-lanes');
  let viewer: DG.Viewer | null = null;
  // the stages when the tracker opened, or since Modify Stages saved them: closing saves when they differ
  let stagesBefore = stages();
  const moved = signal(0);
  const needed = signal(0);
  const room = signal(0);

  dialog.build(() => {
    dialog.add(lanes);
    dialog.addButton('Modify Stages', () => modifyStages(), {start: true});
    // a scrollbar the system may hide (macOS) is no sign that more lanes are there: the line says so
    const fit = new ResizeObserver(() => room.value = lanes.clientWidth);
    fit.observe(lanes);
    dialog.own(() => fit.disconnect());
    const more = computed(() => room.value > 0 && needed.value > room.value ? ' · scroll sideways for more stages' : '');
    dialog.status(computed(() => (moved.value === 0 ? 'Drag a card to another lane to change its stage' :
      `${moved.value} ${moved.value === 1 ? 'molecule' : 'molecules'} moved · saved when the tracker closes`) +
      more.value));
    dialog.onClose(() => closed());
    subscribe(df.onValuesChanged, () => moved.value = countMoved());
  });
  mount();
  return host.show(dialog,
    {modal: true, x: 12, y: 12, width: window.innerWidth - 24, height: window.innerHeight - 24});

  function stages(): (string | null)[] {
    return df.col(TileCategoriesColName)?.toList() ?? [];
  }

  function countMoved(): number {
    return stages().reduce((n, stage, i) => i < stagesBefore.length && stage !== stagesBefore[i] ? n + 1 : n, 0);
  }

  /** The campaign's card layout over this campaign's table, with what names a column the table no longer has dropped.
   * The layout names the table it was made over ("Table"), which another open campaign may hold by now: the cards'
   * fields would read that table's columns. */
  function sketchState(): any | null {
    const stored = app.campaign?.tilesViewerFormSketch;
    if (!stored)
      return null;
    let state: any;
    try {
      state = JSON.parse(stored);
    }
    catch (e) {
      console.error('Failed to parse sketch state string', e);
      return null;
    }
    const table = view.dataFrame;
    state.table = table.name;
    state.elementStates = (state.elementStates ?? [])
      .filter((s: any) => s?.viewerSettings?.column && table.col(s.viewerSettings.column));
    for (const s of state.elementStates) {
      if (s?.viewerSettings?.table)
        s.viewerSettings.table = table.name;
    }
    return state.elementStates.length > 0 ? state : null;
  }

  function laneWidth(sketch: any | null): number {
    const right = Math.max(0, ...(sketch?.elementStates ?? [])
      .map((s: any) => (Number(s?.left) || 0) + (Number(s?.width) || 0)));
    return right > 0 ? right + CARD_FRAME : LANE_WIDTH;
  }

  /** The template's stages, then any other stage a row holds (a typo typed into the grid): no row is left out. */
  function laneNames(): string[] {
    const known = new Set(app.stages);
    return [...app.stages, ...new Set(stages().filter((s): s is string => !!s && !known.has(s)))];
  }

  /** The Tile Viewer over the stages as they are now. */
  function mount(): void {
    unmount();
    const names = laneNames();
    const sketch = sketchState();
    let tiles: DG.Viewer;
    try {
      tiles = view.addViewer(DG.VIEWER.TILE_VIEWER,
        {lanesColumnName: TileCategoriesColName, lanes: names, ...sketch ? {sketchState: sketch} : {}});
      const options = tiles.getOptions();
      options.look.sketchState = sketch;
      tiles.setOptions(options);
      // the Tile Viewer applies a sketch state only through copying a look
      tiles.copyViewersLook(tiles);
    } catch (e) {
      grok.shell.error('Failed to create tiles viewer with preset layout.\n Falling back to default layout.');
      console.error('Failed to create tiles viewer', e);
      tiles = view.addViewer(DG.VIEWER.TILE_VIEWER, {lanesColumnName: TileCategoriesColName, lanes: names});
    }
    needed.value = names.length * laneWidth(sketch);
    Object.assign(tiles.root.style, {minWidth: `${needed.value}px`, height: '100%'});
    lanes.replaceChildren(tiles.root);
    viewer = tiles;
  }

  /** Off the view, so the campaign's saved layout never holds it. */
  function unmount(): void {
    const tiles = viewer;
    viewer = null;
    if (!tiles)
      return;
    try {
      tiles.detach();
      tiles.close();
    }
    catch {
      // the view went first, and the viewer with it
    }
  }

  /** Keeps the card layout with the campaign and saves the moved stages, however the tracker closed (its view too) —
   * unless the app opened another campaign meanwhile: a save would write that one. */
  function closed(): void {
    try {
      const sketch = viewer?.props.sketchState;
      if (sketch && typeof sketch === 'object' && app.campaign)
        app.campaign.tilesViewerFormSketch = JSON.stringify(sketch);
    } catch (e) {
      grok.shell.error('Failed to save sketch state. check the console for more details.');
      console.error('Failed to save sketch state', e);
    }
    unmount();
    if (app.dataFrame === df && countMoved() > 0)
      void app.saveCampaign();
  }

  /** The stages as rows, reordered, added and removed; OK saves them and the tracker shows their lanes. Rows in a
   * removed stage move to the first one, which the dialog says before OK. */
  function modifyStages(): void {
    const modify = new Dialog('Modify Stages', {name: 'modify-stages'});
    modify.build(() => {
      const records = new RecordsInput<{name: string}>({
        name: 'stages', orderable: true, addText: 'Add stage', newItem: () => ({name: ''}),
        value: app.stages.map((name) => ({name})), rowName: (_, i) => `Stage ${i + 1}`,
        columns: [{name: 'name', header: 'Name', input: (value) => new TextInput({value: value ?? '', inline: true})}],
      });
      const kept = computed(() => [...new Set((records.value ?? []).map((r) => r.name?.trim())
        .filter((s): s is string => !!s))]);
      const moves = computed(() => {
        const names = new Set(kept.value);
        const first = kept.value[0];
        const counts = new Map<string, number>();
        for (const stage of stages()) {
          if (stage && !names.has(stage))
            counts.set(stage, (counts.get(stage) ?? 0) + 1);
        }
        return first === undefined || counts.size === 0 ? null : [...counts]
          .map(([stage, n]) => `${n} ${n === 1 ? 'molecule' : 'molecules'} in ${stage} move to ${first}`).join('; ');
      });
      modify.add(records);
      modify.add(statusLine(moves));
      modify.okDisabled(computed(() => kept.value.length === 0 ? 'A campaign needs at least one stage' : null),
        {footer: true});
      modify.onOK(async () => {
        // off the view while the campaign saves its layout
        unmount();
        try {
          // unsaved, the moves stay counted: closing the tracker saves them
          if (await app.setStages(kept.peek())) {
            stagesBefore = stages();
            moved.value = 0;
          }
        } finally {
          mount();
        }
      });
    });
    host.show(modify, {width: 420});
  }
}
