export type Platform = 'mac' | 'other';

export const platformKeyMap: Record<string, { mac: string; other: string }> = {
  Alt: { mac: 'Option', other: 'Alt' },
  Ctrl: { mac: 'Command', other: 'Ctrl' },
  Shift: { mac: 'Shift', other: 'Shift' },
  Enter: { mac: 'Return', other: 'Enter' },

  // Composite actions
  'Ctrl + Home':     { mac: 'Command + ↑', other: 'Ctrl + Home' },  // Jump to first row
  'Ctrl + End':      { mac: 'Command + ↓', other: 'Ctrl + End' },   // Jump to last row
  'End':             { mac: 'Command + →', other: 'End' },          // Last column
  'Home':            { mac: 'Command + ←', other: 'Home' },         // First column
  'Ctrl + Shift + Home': { mac: 'Command + Shift + ↑', other: 'Ctrl + Shift + Home' }, // Select rows above
  'Ctrl + Shift + End':  { mac: 'Command + Shift + ↓', other: 'Ctrl + Shift + End' },  // Select rows below
};

export function getPlatform(): Platform {
  return navigator.platform.toLowerCase().includes('mac') ? 'mac' : 'other';
}