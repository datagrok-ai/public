# Globe

A three-globe WebGL globe with one bar per row at its latitude and longitude, sized and colored by
a numeric column. Registered by Charts as `Globe`, so a feature reaches it as **`globe viewer`**.

## Automation surface (`getWidgetStatus`)

`parts`: `root`, and `canvas` — the WebGL canvas, once attached. A WebGL canvas has no pixels a test
can read, so the reading is the picture.

| Reading | What |
|---|---|
| `points` | point meshes the points layer shows; a row whose color value is missing gets a transparent color and a hidden mesh. `0` until the globe is ready and while the error message is shown |

## `isRenderPending` and `onRendered`

Both come from `RenderSignals` (`utils/utils.ts`), capped at 120 frames here. `pointsData` places
nothing: three-globe digests the data on a 1 ms debounce and binds each point to its mesh as
`__threeObj`, and it hides its whole scene until the globe image has loaded. The flag is held from
`render` until the globe is ready and every point it was handed has its mesh — not until the mesh
count matches, which a new data set of the same size would satisfy at once.
