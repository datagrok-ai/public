[//]: # (doc max 10 lines: public/help/visualize/viewers/grid#color: color cell renderer)

# Color Cell Renderer

Displays a filled color swatch in the cell based on the stored color value.

## Supported formats

All standard CSS color formats are supported:

- **Hex**: `#F00`, `#FF5733`, `#FF573380` (3, 6, or 8 digits)
- **RGB / RGBA**: `rgb(255, 87, 51)`, `rgba(255, 87, 51, 0.5)`
- **HSL / HSLA**: `hsl(14, 100%, 60%)`, `hsla(14, 100%, 60%, 0.5)`
- **Named CSS colors**: `red`, `cornflowerblue`, `teal`

## Auto-detection

The renderer is applied automatically when Datagrok detects a string column whose values are
color codes (hex, rgb, hsl). Named CSS colors are also detected when all sampled values are
valid CSS color names.

## Behavior

- Empty cells and invalid color strings are rendered as blank
- Hovering over a cell shows the raw stored value in a tooltip
- The swatch is padded inside the cell; padding scales down proportionally for small cell sizes so
  the swatch never overlaps neighboring columns
- A subtle border is drawn around the swatch so light colors (e.g., white, pale yellow) remain
  visible against the grid background
