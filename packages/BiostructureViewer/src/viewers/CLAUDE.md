# Viewers

`MolstarViewer` (`molstar-viewer/`, registered as **Biostructure**) and `NglViewer` (`ngl-viewer.ts`,
registered as **NGL**) are `DG.JsViewer`s that queue every load and update on one `PendingSyncer`
(`pending-syncer.ts`, a `PromiseSyncer` that counts what is queued and what its debounces hold).

## Representation (Mol*)

The `representation` property is applied as one extra representation (`sequence-visual-<ref>`) over
each loaded data structure — the viewer's own data, found in `plugin.managers.structure.hierarchy`
— queued on the syncer. A load re-applies a non-default representation to the structure it brings
(data change, current-row reload); the load preset already draws the default. Ligand overlays keep
their row-colored ball-and-stick: they are rebuilt on every current / mouse-over row change.

## Automation surface

`getWidgetStatus()` extends `super.getWidgetStatus()` (the Dart status and any `addStatusProvider`
contributions) with the WebGL `canvas` part once the engine exists and these `values`, all computed
when asked from what the viewer and the engine hold. Mol*'s are built by
`molstar-viewer/molstar-status.ts`.

| Reading | Biostructure (Mol*) | NGL |
|---|---|---|
| `structure loaded` | the plugin holds a structure cell made from the viewer's data | `compList[0]` is a structure component |
| `ligands shown` | distinct rows whose ligands are loaded (selected + current + mouse-over, per the Behaviour switches) | same, as stage components |
| `ligand rows` | those rows, 1-based, ascending, `, `-separated | same |
| `layout expanded` | `plugin.layout.state.isExpanded` | — |
| `controls shown` | `plugin.layout.state.showControls` | — |
| `binding site shown` | the binding-site component and its representation are in the state tree | — |
| `binding site atoms` | the atoms of the binding-site components Mol* built (the radius and Whole Residues decide them); absent while none is shown | — |
| `components` | — | `stage.compList.length` |
| `representation` | type of the representation the viewer applied over the data structure if it built, else of the first one the load preset made | type of the first representation of `compList[0]` |

**Render pending.** `isRenderPending` is true while anything is queued on the syncer (data load,
current-row reload, ligand rebuild, representation apply, layout update, `invalidate`), while one of
its debounces (`setData`, ligand rebuild, current-row reload) holds a request its handler has not
queued yet — a request lives only as long as its handler's subscription, so detaching clears it —
and while the engine is still working (Mol* `behaviors.state.isUpdating`, NGL `stage.tasks.count`).
`onRendered` fires on `invalidate` and every time the syncer's queue drains.

**Overlay toggles.** Mol* keeps the on/off state of its viewport buttons (`Screenshot / State
Snapshot`, `Toggle Controls Panel`, `Toggle Expanded Viewport`, `Settings / Controls Info`,
`Toggle Selection Mode`) only in `msp-btn-link-toggle-on|off`, and the package's own `Binding site`
button follows the same classes; the bdd library reads `msp-btn-link-toggle-on` as `selected` when
asked, so nothing watches the DOM.

**Named elements.** The NGL file preview and the NGL double-click view put
`data-u2-name="ngl-host"` on their `.d4-ngl-viewer` host. The preview view is named after the file;
the double-click view stays `NGL`, since a file importer receives only the content.
