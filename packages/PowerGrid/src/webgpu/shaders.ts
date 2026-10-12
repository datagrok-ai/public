/* eslint-disable max-len */
/** WGSL of the WebGPU scatter plot renderer: the WebGL shaders (../webgl/shaders.ts) ported line by line, so that both
 * backends place, shape and color markers alike. Rows come from a list of row indexes; per-row values from storage
 * buffers addressed by row index, or the one value every row shares. */

/** The bits of `Frame.flags`. */
export const FLAGS = {rowColors: 1, rowSizes: 2, rowTypes: 4, selection: 8, markerCache: 16};
/** The size of the `Frame` struct. */
export const FRAME_BYTES = 96;
export const PICK_WORKGROUP = 256;

const PRELUDE = /* wgsl */`
struct Frame {
  xform: vec4f,        // screen x = xform.x + xform.y * x, screen y = xform.z + xform.w * y
  viewBox: vec4f,      // left, top, right, bottom
  origin: vec2f,       // the GPU canvas' top-left in device pixels of the viewer canvas
  size: vec2f,         // the GPU canvas in device pixels
  mouse: vec2f,
  dpr: f32,
  alpha: f32,
  border: f32,
  markerSize: f32,
  color: u32,          // ARGB, every row's when colors are not per row
  markerType: u32,
  selectedColor: u32,
  flags: u32,
  count: u32,
  pad: u32,
};

const ROW_COLORS = ${FLAGS.rowColors}u;
const ROW_SIZES = ${FLAGS.rowSizes}u;
const ROW_TYPES = ${FLAGS.rowTypes}u;
const SELECTION = ${FLAGS.selection}u;
const MARKER_CACHE = ${FLAGS.markerCache}u;

const CIRCLE = 0u; const CIRCLE_BORDER = 1u; const SQUARE = 2u; const SQUARE_BORDER = 3u; const CROSS_BORDER = 4u;
const CROSS_X_BORDER = 5u; const OUTLIER = 6u; const DIAMOND = 7u; const DIAMOND_BORDER = 8u; const TRIANGLE_TOP = 9u;
const TRIANGLE_RIGHT = 10u; const TRIANGLE_BOTTOM = 11u; const TRIANGLE_LEFT = 12u; const TRIANGLE_TOP_BORDER = 13u;
const TRIANGLE_RIGHT_BORDER = 14u; const TRIANGLE_BOTTOM_BORDER = 15u; const TRIANGLE_LEFT_BORDER = 16u;
const STAR = 17u; const DOT = 18u; const GRADIENT = 19u;

@group(0) @binding(0) var<uniform> frame: Frame;
@group(0) @binding(1) var<storage, read> rows: array<u32>;
@group(0) @binding(2) var<storage, read> xs: array<f32>;
@group(0) @binding(3) var<storage, read> ys: array<f32>;
@group(0) @binding(5) var<storage, read> sizes: array<f32>;
@group(0) @binding(6) var<storage, read> types: array<u32>;

fn screen(row: u32) -> vec2f {
  return vec2f(frame.xform.x + frame.xform.y * xs[row], frame.xform.z + frame.xform.w * ys[row]);
}
fn rowSize(row: u32) -> f32 {
  if ((frame.flags & ROW_SIZES) != 0u) { return sizes[row]; }
  return frame.markerSize;
}
fn rowType(row: u32) -> u32 {
  if ((frame.flags & ROW_TYPES) != 0u) { return types[row]; }
  return frame.markerType;
}
`;

const DRAW_PRELUDE = /* wgsl */`
${PRELUDE}
@group(0) @binding(4) var<storage, read> colors: array<u32>;
@group(0) @binding(7) var<storage, read> selection: array<u32>;

fn inViewBox(p: vec2f) -> bool {
  return p.x >= frame.viewBox.x && p.x <= frame.viewBox.z && p.y >= frame.viewBox.y && p.y <= frame.viewBox.w;
}
fn clip(device: vec2f) -> vec4f {
  let p = (device - frame.origin) / frame.size * 2.0 - 1.0;
  return vec4f(p.x, -p.y, 0.0, 1.0);
}
const CULLED = vec4f(2.0, 2.0, 2.0, 1.0);

// selected rows take the selection color, as ScatterPlotCore.getColor paints them
fn rowColor(row: u32) -> vec3f {
  var c = frame.color;
  if ((frame.flags & ROW_COLORS) != 0u) { c = colors[row]; }
  if ((frame.flags & SELECTION) != 0u && ((selection[row >> 5u] >> (row & 31u)) & 1u) != 0u) { c = frame.selectedColor; }
  return vec3f(f32((c >> 16u) & 255u), f32((c >> 8u) & 255u), f32(c & 255u)) / 255.0;
}

// the markers whose body is a 3-pixel stroke rather than a fill
fn isStroked(kind: u32) -> bool {
  return kind == CIRCLE_BORDER || kind == SQUARE_BORDER || kind == CROSS_BORDER || kind == CROSS_X_BORDER ||
    kind == DIAMOND_BORDER || (kind >= TRIANGLE_TOP_BORDER && kind <= TRIANGLE_LEFT_BORDER);
}
`;

/** Draws markers: the dot renderer's one opaque device pixel per row, and shapes as quads around the markers. */
export const DRAW_SHADER = /* wgsl */`
${DRAW_PRELUDE}

struct DotVarying {
  @builtin(position) position: vec4f,
  @location(0) @interpolate(flat) color: vec3f,
};

@vertex fn dotVertex(@builtin(vertex_index) i: u32) -> DotVarying {
  let row = rows[i];
  let p = screen(row);
  var v: DotVarying;
  v.position = select(CULLED, clip(floor(p * frame.dpr) + 0.5), inViewBox(p));
  v.color = rowColor(row);
  return v;
}

@fragment fn dotFragment(v: DotVarying) -> @location(0) vec4f {
  return vec4f(v.color, 1.0);
}

struct ShapeVarying {
  @builtin(position) position: vec4f,
  @location(0) @interpolate(flat) color: vec3f,
  @location(1) @interpolate(flat) center: vec2f,
  @location(2) @interpolate(flat) size: f32,
  @location(3) @interpolate(flat) kind: u32,
  @location(4) @interpolate(flat) dotOffset: vec2f,
};

@vertex fn shapeVertex(@builtin(vertex_index) corner: u32, @builtin(instance_index) i: u32) -> ShapeVarying {
  let row = rows[i];
  let p = screen(row);
  let rawSize = rowSize(row);
  let kind = rowType(row);
  // Marker.draw truncates the size, and blits the markers it caches, under 50 device pixels, at whole CSS pixels
  let center = select(p, floor(p + 0.5), (frame.flags & MARKER_CACHE) != 0u && rawSize * frame.dpr < 50.0);
  // the widest stroke, mitred at the sharpest corner (60 degrees), reaches its full width past the outline
  let stroke = max(frame.border, select(0.0, 3.0, isStroked(kind)));
  let half = ceil((floor(rawSize) * 0.5 + stroke) * frame.dpr + 1.0);
  let offset = vec2f(f32(corner & 1u), f32(corner >> 1u)) * 2.0 - 1.0;
  var v: ShapeVarying;
  v.position = select(CULLED, clip(center * frame.dpr + offset * half), inViewBox(p));
  v.color = rowColor(row);
  v.center = center * frame.dpr - frame.origin;
  v.size = floor(rawSize);
  v.kind = kind;
  v.dotOffset = ceil(p) - center;
  return v;
}

// The shapes below are measured in CSS pixels, y down, from the marker center; negative is inside.
const NONE = 1e6;
fn box(p: vec2f, b: vec2f) -> f32 { let d = abs(p) - b; return length(max(d, vec2f(0.0))) + min(max(d.x, d.y), 0.0); }

// Convex outlines as the max of edge distances: their offsets have the mitred corners of a canvas stroke.
fn square(p: vec2f, h: f32) -> f32 { return max(abs(p.x), abs(p.y)) - h; }
fn diamond(p: vec2f, h: f32) -> f32 { return (abs(p.x) + abs(p.y) - h) * 0.70710678; }
fn triangleTop(p: vec2f, s: f32) -> f32 {
  let h = s * 0.8660254;
  return max(p.y - h * 0.5, 0.8660254 * abs(p.x) - 0.5 * (p.y + h * 0.5));
}
fn triangle(kind: u32, p: vec2f, s: f32) -> f32 {
  if (kind == TRIANGLE_RIGHT || kind == TRIANGLE_RIGHT_BORDER) { return triangleTop(vec2f(p.y, -p.x), s); }
  if (kind == TRIANGLE_BOTTOM || kind == TRIANGLE_BOTTOM_BORDER) { return triangleTop(-p, s); }
  if (kind == TRIANGLE_LEFT || kind == TRIANGLE_LEFT_BORDER) { return triangleTop(vec2f(-p.y, p.x), s); }
  return triangleTop(p, s);
}

// Five spikes from the top, the inner radius half the outer one.
fn star(p: vec2f, r: f32) -> f32 {
  var dirs = array<vec2f, 10>(vec2f(0.0, -1.0), vec2f(0.58778525, -0.80901699), vec2f(0.95105652, -0.30901699),
    vec2f(0.95105652, 0.30901699), vec2f(0.58778525, 0.80901699), vec2f(0.0, 1.0), vec2f(-0.58778525, 0.80901699),
    vec2f(-0.95105652, 0.30901699), vec2f(-0.95105652, -0.30901699), vec2f(-0.58778525, -0.80901699));
  var d = 1e9;
  var s = 1.0;
  var j = 9u;
  for (var i = 0u; i < 10u; i++) {
    let vi = dirs[i] * select(r * 0.5, r, i % 2u == 0u);
    let vj = dirs[j] * select(r * 0.5, r, j % 2u == 0u);
    let e = vj - vi;
    let w = p - vi;
    let b = w - e * clamp(dot(w, e) / dot(e, e), 0.0, 1.0);
    d = min(d, dot(b, b));
    let c0 = p.y >= vi.y;
    let c1 = p.y < vj.y;
    let c2 = e.x * w.y > e.y * w.x;
    if ((c0 && c1 && c2) || (!c0 && !c1 && !c2)) { s = -s; }
    j = i;
  }
  return s * sqrt(d);
}

// Two butt-capped strokes of the given half width: the arms of a cross, or of an x when rotated.
fn arms(p: vec2f, h: f32, halfWidth: f32) -> f32 { return min(box(p, vec2f(halfWidth, h)), box(p, vec2f(h, halfWidth))); }
fn armsX(p: vec2f, h: f32, halfWidth: f32) -> f32 {
  return arms(vec2f(p.x + p.y, p.y - p.x) * 0.70710678, h * 1.41421356, halfWidth);
}

fn outline(kind: u32, p: vec2f, s: f32) -> f32 {
  let h = s * 0.5;
  if (kind == SQUARE || kind == SQUARE_BORDER) { return square(p, h); }
  if (kind == DIAMOND || kind == DIAMOND_BORDER) { return diamond(p, h); }
  if (kind >= TRIANGLE_TOP && kind <= TRIANGLE_LEFT_BORDER) { return triangle(kind, p, s); }
  if (kind == STAR) { return star(p, h); }
  return length(p) - h;
}

// Signed distances from [p] to the marker body and to its border.
fn marker(v: ShapeVarying, p: vec2f) -> vec2f {
  let s = v.size;
  let h = s * 0.5;
  if (v.kind == DOT) { return vec2f(box(p - v.dotOffset - 0.5, vec2f(0.5)), NONE); }
  if (v.kind == OUTLIER) { return vec2f(armsX(p, s * 0.4, 0.5), NONE); }
  let w = frame.border * 0.5;
  if (v.kind == CROSS_BORDER || v.kind == CROSS_X_BORDER) {
    if (v.kind == CROSS_X_BORDER) { return vec2f(armsX(p, h, 1.5), select(armsX(p, h, w), NONE, w == 0.0)); }
    return vec2f(arms(p, h, 1.5), select(arms(p, h, w), NONE, w == 0.0));
  }
  let d = outline(v.kind, p, s);
  return vec2f(select(d, abs(d) - 1.5, isStroked(v.kind)), select(abs(d) - w, NONE, w == 0.0));
}

// Pixel coverage at signed distances [d], ramping over a pixel, or over a quarter pixel's share of it.
fn coverage(d: vec2f, ramps: f32) -> vec2f { return clamp(0.5 - d * frame.dpr * ramps, vec2f(0.0), vec2f(1.0)); }

@fragment fn shapeFragment(v: ShapeVarying) -> @location(0) vec4f {
  let p = (v.position.xy - v.center) / frame.dpr;
  var o: vec4f;
  if (v.kind == GRADIENT) {
    // a canvas gradient fades to transparent black interpolating unpremultiplied colors, which darkens it
    let g = max(1.0 - length(p) / (v.size * 0.5), 0.0);
    o = vec4f(v.color * g, 1.0) * (g * frame.alpha);
  }
  else {
    let d = marker(v, p);
    // beyond the reach of any sample: no sample moves more than 0.36 pixels and none ramps over more than 0.5
    if (min(d.x, d.y) * frame.dpr > 0.86) { discard; }
    // The canvas anti-aliases rectangles and straight strokes by the distance to the edge, one ramp per pixel, and
    // curves and paths closer to the covered area, which four samples a pixel approximate near the edges.
    var c: vec2f;
    let area = v.kind == CIRCLE || v.kind == CIRCLE_BORDER || v.kind == STAR || v.kind == CROSS_BORDER;
    if (area && min(abs(d.x), abs(d.y)) * frame.dpr < 0.86) {
      let q = 0.25 / frame.dpr;
      c = (coverage(marker(v, p + vec2f(-q, -q)), 2.0) + coverage(marker(v, p + vec2f(q, -q)), 2.0) +
        coverage(marker(v, p + vec2f(-q, q)), 2.0) + coverage(marker(v, p + vec2f(q, q)), 2.0)) * 0.25;
    }
    else {
      c = coverage(d, 1.0);
    }
    // The border is the body color darkened by 50 per channel (Color.darken), stroked over the body.
    let body = vec4f(v.color, 1.0) * (c.x * frame.alpha);
    let edge = vec4f(max(v.color - 50.0 / 255.0, vec3f(0.0)), 1.0) * (c.y * frame.alpha);
    o = edge + body * (1.0 - edge.a);
  }
  if (o.a <= 0.0) { discard; }
  return o;
}
`;

/** Hit testing: of the rows whose marker contains the mouse point, by `Marker.hitTest`, keeps the highest, as the
 * core's CPU hit test, scanning rows from the last, does. */
export const PICK_SHADER = /* wgsl */`
${PRELUDE}
@group(0) @binding(8) var<storage, read_write> result: atomic<i32>;

@compute @workgroup_size(${PICK_WORKGROUP})
fn pick(@builtin(global_invocation_id) id: vec3u, @builtin(num_workgroups) groups: vec3u) {
  let i = id.x + id.y * groups.x * ${PICK_WORKGROUP}u;
  if (i >= frame.count) { return; }
  let row = rows[i];
  let d = abs(frame.mouse - screen(row));
  let h = rowSize(row) * 0.5;
  let kind = rowType(row);
  var hit = d.x <= h && d.y <= h;
  if (kind == DOT) { hit = d.x <= 1.0 && d.y <= 1.0; }
  else if (kind == CIRCLE || kind == CIRCLE_BORDER || kind == GRADIENT || kind == STAR || kind == OUTLIER) { hit = dot(d, d) <= h * h; }
  if (hit) { atomicMax(&result, i32(row)); }
}
`;
