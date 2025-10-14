// Use splitV and splitH to divide the area into resizable containers
// If you want a non-even split, specify it (in pixels - this is important) in the element width or height style.

ui.splitV([
  ui.splitH([
    ui.divText('Top Left', { style: { width: '100px'}}),
    ui.divText('Top Right'),
  ], { style: { height: '40px'}}, true),
  ui.splitH([
    ui.divText('Panel 1'),
    ui.divText('Panel 2')
  ], {}, true)
], {}, true)