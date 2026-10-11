/* eslint-disable max-len */
/** Shaders of the WebGL scatter plot renderer. Positions arrive as world values relative to an origin, and
 * `u_xform` maps them to the viewer canvas in CSS pixels exactly as the core's `toScreen` does. */

const VERTEX_COMMON = `#version 300 es
precision highp float;
precision highp int;
layout(location = 0) in float a_x;
layout(location = 1) in float a_y;
layout(location = 2) in vec4 a_color;  // the core's ARGB color read as bytes: (b, g, r, a)
layout(location = 3) in float a_size;
layout(location = 4) in float a_type;
layout(location = 5) in float a_selected;
uniform vec4 u_xform;    // screen x = xform.x + xform.y * a_x, screen y = xform.z + xform.w * a_y
uniform vec4 u_viewBox;  // left, top, right, bottom
uniform vec2 u_origin;   // the GL canvas' top-left in device pixels of the viewer canvas
uniform vec2 u_size;     // the GL canvas in device pixels
uniform float u_dpr;
uniform bool u_selection;     // selected rows take the selection color, as ScatterPlotCore.getColor paints them
uniform vec3 u_selectedColor;

vec2 screen() { return vec2(u_xform.x + u_xform.y * a_x, u_xform.z + u_xform.w * a_y); }
vec3 rowColor() { return u_selection && a_selected > 0.5 ? u_selectedColor : a_color.zyx; }
bool inViewBox(vec2 p) { return p.x >= u_viewBox.x && p.x <= u_viewBox.z && p.y >= u_viewBox.y && p.y <= u_viewBox.w; }
vec4 clip(vec2 device) {
  vec2 p = (device - u_origin) / u_size * 2.0 - 1.0;
  return vec4(p.x, -p.y, 0.0, 1.0);
}
const vec4 CULLED = vec4(2.0, 2.0, 2.0, 1.0);
`;

/** The dot renderer: one opaque device pixel per row, where the core's pixel-buffer renderer puts it. */
export const DOT_VERTEX = `${VERTEX_COMMON}
out vec4 v_color;
void main() {
  vec2 p = screen();
  gl_PointSize = 1.0;
  gl_Position = inViewBox(p) ? clip(floor(p * u_dpr) + 0.5) : CULLED;
  v_color = vec4(rowColor(), 1.0);
}`;

export const DOT_FRAGMENT = `#version 300 es
precision highp float;
in vec4 v_color;
out vec4 o;
void main() { o = v_color; }`;

const MARKER_TYPES = `
const int CIRCLE = 0, CIRCLE_BORDER = 1, SQUARE = 2, SQUARE_BORDER = 3, CROSS_BORDER = 4, CROSS_X_BORDER = 5,
  OUTLIER = 6, DIAMOND = 7, DIAMOND_BORDER = 8, TRIANGLE_TOP = 9, TRIANGLE_RIGHT = 10, TRIANGLE_BOTTOM = 11,
  TRIANGLE_LEFT = 12, TRIANGLE_TOP_BORDER = 13, TRIANGLE_RIGHT_BORDER = 14, TRIANGLE_BOTTOM_BORDER = 15,
  TRIANGLE_LEFT_BORDER = 16, STAR = 17, DOT = 18, GRADIENT = 19;

// The markers whose body is a 3-pixel stroke rather than a fill.
bool isStroked(int type) {
  return type == CIRCLE_BORDER || type == SQUARE_BORDER || type == CROSS_BORDER || type == CROSS_X_BORDER ||
    type == DIAMOND_BORDER || (type >= TRIANGLE_TOP_BORDER && type <= TRIANGLE_LEFT_BORDER);
}
`;

/** Marker shapes, drawn as point sprites: the fragment shader reproduces `Marker.draw` with signed distances. */
export const SHAPE_VERTEX = `${VERTEX_COMMON}
${MARKER_TYPES}
uniform float u_border;
uniform bool u_markerCache;  // the viewer canvas is large enough for Marker.draw to cache markers
flat out vec3 v_color;
flat out vec2 v_center;
flat out float v_size;
flat out int v_type;
flat out vec2 v_dot;
void main() {
  vec2 p = screen();
  // Marker.draw truncates the size, and blits the markers it caches, under 50 device pixels, at whole CSS pixels
  vec2 center = u_markerCache && a_size * u_dpr < 50.0 ? floor(p + 0.5) : p;
  v_center = center * u_dpr - u_origin;
  v_size = floor(a_size);
  v_type = int(a_type);
  v_dot = ceil(p) - center;
  v_color = rowColor();
  // the widest stroke, mitred at the sharpest corner (60 degrees), reaches its full width past the outline
  float stroke = max(u_border, isStroked(v_type) ? 3.0 : 0.0);
  gl_PointSize = 2.0 * ceil((v_size * 0.5 + stroke) * u_dpr + 1.0);
  gl_Position = inViewBox(p) ? clip(center * u_dpr) : CULLED;
}`;

export const SHAPE_FRAGMENT = `#version 300 es
precision highp float;
precision highp int;
${MARKER_TYPES}
uniform vec2 u_size;
uniform float u_dpr;
uniform float u_border;
uniform float u_alpha;
flat in vec3 v_color;
flat in vec2 v_center;
flat in float v_size;
flat in int v_type;
flat in vec2 v_dot;
out vec4 o;

// The shapes below are measured in CSS pixels, y down, from the marker center; negative is inside.
const float NONE = 1e6;
float box(vec2 p, vec2 b) { vec2 d = abs(p) - b; return length(max(d, 0.0)) + min(max(d.x, d.y), 0.0); }

// Convex outlines as the max of edge distances: their offsets have the mitred corners of a canvas stroke.
float square(vec2 p, float h) { return max(abs(p.x), abs(p.y)) - h; }
float diamond(vec2 p, float h) { return (abs(p.x) + abs(p.y) - h) * 0.70710678; }
float triangleTop(vec2 p, float s) {
  float h = s * 0.8660254;
  return max(p.y - h * 0.5, 0.8660254 * abs(p.x) - 0.5 * (p.y + h * 0.5));
}
float triangle(int type, vec2 p, float s) {
  if (type == TRIANGLE_RIGHT || type == TRIANGLE_RIGHT_BORDER) return triangleTop(vec2(p.y, -p.x), s);
  if (type == TRIANGLE_BOTTOM || type == TRIANGLE_BOTTOM_BORDER) return triangleTop(-p, s);
  if (type == TRIANGLE_LEFT || type == TRIANGLE_LEFT_BORDER) return triangleTop(vec2(-p.y, p.x), s);
  return triangleTop(p, s);
}

// Five spikes from the top, the inner radius half the outer one.
const vec2 STAR_DIRS[10] = vec2[](vec2(0.0, -1.0), vec2(0.58778525, -0.80901699), vec2(0.95105652, -0.30901699),
  vec2(0.95105652, 0.30901699), vec2(0.58778525, 0.80901699), vec2(0.0, 1.0), vec2(-0.58778525, 0.80901699),
  vec2(-0.95105652, 0.30901699), vec2(-0.95105652, -0.30901699), vec2(-0.58778525, -0.80901699));
float star(vec2 p, float r) {
  float d = 1e9;
  float s = 1.0;
  for (int i = 0, j = 9; i < 10; j = i, i++) {
    vec2 vi = STAR_DIRS[i] * (i % 2 == 0 ? r : r * 0.5);
    vec2 vj = STAR_DIRS[j] * (j % 2 == 0 ? r : r * 0.5);
    vec2 e = vj - vi;
    vec2 w = p - vi;
    vec2 b = w - e * clamp(dot(w, e) / dot(e, e), 0.0, 1.0);
    d = min(d, dot(b, b));
    bvec3 c = bvec3(p.y >= vi.y, p.y < vj.y, e.x * w.y > e.y * w.x);
    if (all(c) || all(not(c))) s = -s;
  }
  return s * sqrt(d);
}

// Two butt-capped strokes of the given half width: the arms of a cross, or of an x when rotated.
float arms(vec2 p, float h, float halfWidth) { return min(box(p, vec2(halfWidth, h)), box(p, vec2(h, halfWidth))); }
float armsX(vec2 p, float h, float halfWidth) {
  return arms(vec2(p.x + p.y, p.y - p.x) * 0.70710678, h * 1.41421356, halfWidth);
}

float outline(int type, vec2 p, float s) {
  float h = s * 0.5;
  if (type == SQUARE || type == SQUARE_BORDER) return square(p, h);
  if (type == DIAMOND || type == DIAMOND_BORDER) return diamond(p, h);
  if (type >= TRIANGLE_TOP && type <= TRIANGLE_LEFT_BORDER) return triangle(type, p, s);
  if (type == STAR) return star(p, h);
  return length(p) - h;
}

/** Signed distances from [p] to the marker body and to its border. */
vec2 marker(vec2 p) {
  float s = v_size;
  float h = s * 0.5;
  if (v_type == DOT)
    return vec2(box(p - v_dot - 0.5, vec2(0.5)), NONE);
  if (v_type == OUTLIER)
    return vec2(armsX(p, s * 0.4, 0.5), NONE);
  float w = u_border * 0.5;
  if (v_type == CROSS_BORDER || v_type == CROSS_X_BORDER) {
    bool x = v_type == CROSS_X_BORDER;
    return vec2(x ? armsX(p, h, 1.5) : arms(p, h, 1.5), w == 0.0 ? NONE : x ? armsX(p, h, w) : arms(p, h, w));
  }
  float d = outline(v_type, p, s);
  return vec2(isStroked(v_type) ? abs(d) - 1.5 : d, w == 0.0 ? NONE : abs(d) - w);
}

/** Pixel coverage at signed distances [d], ramping over a pixel, or over a quarter pixel's share of it. */
vec2 coverage(vec2 d, float ramps) { return clamp(0.5 - d * u_dpr * ramps, 0.0, 1.0); }

void main() {
  vec2 p = (vec2(gl_FragCoord.x, u_size.y - gl_FragCoord.y) - v_center) / u_dpr;
  if (v_type == GRADIENT) {
    // a canvas gradient fades to transparent black interpolating unpremultiplied colors, which darkens it
    float g = max(1.0 - length(p) / (v_size * 0.5), 0.0);
    o = vec4(v_color * g, 1.0) * (g * u_alpha);
  }
  else {
    vec2 d = marker(p);
    // beyond the reach of any sample: no sample moves more than 0.36 pixels and none ramps over more than 0.5
    if (min(d.x, d.y) * u_dpr > 0.86)
      discard;
    // The canvas anti-aliases rectangles and straight strokes by the distance to the edge, one ramp per pixel, and
    // curves and paths closer to the covered area, which four samples a pixel approximate near the edges.
    vec2 c;
    bool area = v_type == CIRCLE || v_type == CIRCLE_BORDER || v_type == STAR || v_type == CROSS_BORDER;
    if (area && min(abs(d.x), abs(d.y)) * u_dpr < 0.86) {
      float q = 0.25 / u_dpr;
      c = (coverage(marker(p + vec2(-q, -q)), 2.0) + coverage(marker(p + vec2(q, -q)), 2.0) +
        coverage(marker(p + vec2(-q, q)), 2.0) + coverage(marker(p + vec2(q, q)), 2.0)) * 0.25;
    }
    else
      c = coverage(d, 1.0);
    // The border is the body color darkened by 50 per channel (Color.darken), stroked over the body.
    vec4 bodyColor = vec4(v_color, 1.0) * (c.x * u_alpha);
    vec4 edgeColor = vec4(max(v_color - 50.0 / 255.0, 0.0), 1.0) * (c.y * u_alpha);
    o = edgeColor + bodyColor * (1.0 - edgeColor.a);
  }
  if (o.a <= 0.0)
    discard;
}`;

/** Hit testing: draws the rows into an offscreen target scissored to the mouse pixel and keeps the last row whose
 * marker contains the mouse point, by `Marker.hitTest`. Rows go in ascending order, so the highest one wins. */
export const PICK_VERTEX = `${VERTEX_COMMON}
flat out vec2 v_p;
flat out float v_size;
flat out int v_type;
flat out uint v_row;
void main() {
  vec2 p = screen();
  v_p = p;
  v_size = a_size;
  v_type = int(a_type);
  v_row = uint(gl_VertexID);
  gl_PointSize = 2.0 * ceil((max(a_size * 0.5, 1.0) + 2.0) * u_dpr);
  gl_Position = clip(p * u_dpr);
}`;

export const PICK_FRAGMENT = `#version 300 es
precision highp float;
precision highp int;
${MARKER_TYPES}
uniform vec2 u_mouse;
flat in vec2 v_p;
flat in float v_size;
flat in int v_type;
flat in uint v_row;
out vec4 o;
void main() {
  vec2 d = abs(u_mouse - v_p);
  float h = v_size * 0.5;
  bool isRound = v_type == CIRCLE || v_type == CIRCLE_BORDER || v_type == GRADIENT || v_type == STAR || v_type == OUTLIER;
  bool hit = v_type == DOT ? d.x <= 1.0 && d.y <= 1.0 : isRound ? dot(d, d) <= h * h : d.x <= h && d.y <= h;
  if (!hit)
    discard;
  uint id = v_row + 1u;
  o = vec4(uvec4(id, id >> 8, id >> 16, id >> 24) & 255u) / 255.0;
}`;
