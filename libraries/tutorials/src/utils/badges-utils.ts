import * as grok from 'datagrok-api/grok';
import {Tutorial} from '../tutorial';

export interface AchievementBadge {
  icon: string;
  description: string;
  unlocked: boolean;
}

export interface BadgeRegistry {
  [id: string]: string;
}

export const STORAGE_NAME = 'achievement-badges';
export const AWARD_ASSIGNED_EVENT = 'd4-award-assigned-event';

export function createBadgeRegistry(tutorials: Tutorial[]): BadgeRegistry {
  return Object.fromEntries(
    tutorials.map((tutorial) => [
      tutorial.name,
      JSON.stringify({
        icon: tutorial.icon,
        description: `Completed ${tutorial.name} tutorial`,
        unlocked: false,
      }),
    ]),
  );
}

export function loadBadges(tutorials: Tutorial[]): void {
  const stored = grok.userSettings.get(STORAGE_NAME);
  if (stored) return;
  const badgeRegistry = createBadgeRegistry(tutorials);
  grok.userSettings.addAll(STORAGE_NAME, badgeRegistry);
}

export function awardBadge(badgeName: string): AchievementBadge | null {
  const all = grok.userSettings.get(STORAGE_NAME);
  if (!all) return null;

  const raw = all[badgeName];
  if (!raw) return null;

  const storedBadge: AchievementBadge = JSON.parse(raw);
  if (storedBadge.unlocked) return storedBadge;
  storedBadge.unlocked = true;
  all[badgeName] = JSON.stringify(storedBadge);

  grok.userSettings.put(STORAGE_NAME, all);
  fireAwardAssignedEvent(storedBadge.icon);
  return storedBadge;
}

export function fireAwardAssignedEvent(icon: string) {
  grok.events.fireCustomEvent(AWARD_ASSIGNED_EVENT, {icon: icon});
}
