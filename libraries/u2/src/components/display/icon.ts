/* Font Awesome icons, in the platform's own class convention (`grok-icon fal fa-{name}`), so
   in-app they render with the FA Pro sheet the app already loads. Standalone hosts mark
   `<html class="u2-standalone">` and load css/icons.css, which maps the same classes onto the
   vendored FA Pro webfonts. */
import {Scope} from '../../core/scope.js';
import {Tooltip} from '../../core/tooltip.js';

const STANDALONE_FACES = [
  '300 12px "Font Awesome 5 Pro Light"',
  '400 12px "Font Awesome 5 Pro Regular"',
  '900 12px "Font Awesome 5 Pro Solid"',
  '400 12px "Font Awesome 5 Brands U2"',
];
let facesRequested = false;

// font-display: block keeps ::before glyphs invisible until the face loads, and some hosts
// never trigger the load implicitly — request the vendored faces once, on first icon use
function requestStandaloneFaces(): void {
  if (facesRequested || !document.fonts || !document.documentElement.classList.contains('u2-standalone'))
    return;
  facesRequested = true;
  for (const f of STANDALONE_FACES)
    document.fonts.load(f).catch(() => undefined);
}

export type IconVariant = 'light' | 'regular' | 'solid';

export interface IconOptions {
  variant?: IconVariant;
  size?: 'small' | 'medium' | 'large';
  /** Bound through the ambient Scope when there is one, so it dies with its owner; otherwise
   * it falls back to the native `title` attribute. */
  tooltip?: string;
  cls?: string;
}

const VARIANT_CLASS: Record<IconVariant, string> = {light: 'fal', regular: 'far', solid: 'fas'};

export function icon(name: string, options?: IconOptions): HTMLElement {
  requestStandaloneFaces();
  const el = document.createElement('i');
  el.className = ['grok-icon', 'u2-icon', VARIANT_CLASS[options?.variant ?? 'light'], `fa-${name}`,
    options?.size ? `u2-icon-${options.size}` : undefined, options?.cls]
    .filter((c) => c).join(' ');
  el.dataset.u2 = 'icon';
  const tooltip = options?.tooltip;
  if (tooltip) {
    el.setAttribute('aria-label', tooltip);
    const scope = Scope.ambient;
    if (scope)
      Tooltip.bind(el, tooltip, scope);
    else
      el.title = tooltip;
  }
  return el;
}
