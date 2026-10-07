/** The portal's three icons, one place for all of them: the app, a portfolio and a program. Single-colour
 * Font Awesome glyphs, as the platform draws its own, so Browse, the Workspace on Home and the app agree. */
import * as ui from 'datagrok-api/ui';

export const PORTAL_ICONS = {
  app: 'project-diagram',
  portfolio: 'briefcase',
  program: 'flask',
} as const;

export type PortalIconKind = keyof typeof PORTAL_ICONS;

/** The icon of [kind], muted and not clickable; [tooltip] says what it stands for. */
export function portalIcon(kind: PortalIconKind, tooltip: string | null = null): HTMLElement {
  const i = ui.iconFA(PORTAL_ICONS[kind], null, tooltip);
  i.classList.add('pp-icon', `pp-icon-${kind}`);
  return i;
}
