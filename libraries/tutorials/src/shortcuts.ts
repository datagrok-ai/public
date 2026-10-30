import * as ui from 'datagrok-api/ui';


type Keyboard = 'Windows' | 'Mac';

const IS_MAC = navigator.platform.startsWith('Mac');

const MAC_MODIFIERS: Record<string, string> = {Ctrl: '⌘', Alt: '⌥'};

const MAC_NAVIGATION: Record<string, string> = {
  'Home': 'Ctrl+←',
  'End': 'Ctrl+→',
  'Ctrl+Home': 'Ctrl+↑',
  'Ctrl+End': 'Ctrl+↓',
  'Ctrl+Shift+Home': 'Ctrl+Shift+↑',
  'Ctrl+Shift+End': 'Ctrl+Shift+↓',
};

export class Shortcuts {
  private static readonly STORAGE_KEY = 'tutorials-keyboard';
  private static _keyboard: Keyboard | null = null;

  private static get keyboard(): Keyboard {
    if (this._keyboard === null) {
      const saved = localStorage.getItem(this.STORAGE_KEY);
      this._keyboard = saved === 'Windows' || saved === 'Mac' ? saved : IS_MAC ? 'Mac' : 'Windows';
    }
    return this._keyboard;
  }

  private static get otherKeyboard(): Keyboard {
    return this.keyboard === 'Mac' ? 'Windows' : 'Mac';
  }

  private static text(keys: string, keyboard: Keyboard): string {
    const windowsKeys = keys.split(/\s*\+\s*/).join('+');
    if (keyboard === 'Windows')
      return windowsKeys;
    const macKeys = (IS_MAC ? MAC_NAVIGATION[windowsKeys] : undefined) ?? windowsKeys;
    return macKeys.split('+').map((k) => MAC_MODIFIERS[k] ?? k).join('+');
  }

  static render(root: ParentNode, onSwitch: () => void): void {
    for (const kbd of root.querySelectorAll<HTMLElement>('kbd')) {
      const bound = kbd.dataset.keys !== undefined;
      const keys = kbd.dataset.keys ??= kbd.textContent ?? '';
      kbd.textContent = this.text(keys, this.keyboard);
      if (bound)
        continue;
      ui.tooltip.bind(kbd, () => `On ${this.otherKeyboard}: ${this.text(keys, this.otherKeyboard)}. ` +
        `Click to show ${this.otherKeyboard} shortcuts.`);
      kbd.onclick = () => {
        this._keyboard = this.otherKeyboard;
        localStorage.setItem(this.STORAGE_KEY, this._keyboard);
        ui.tooltip.hide();
        onSwitch();
      };
    }
  }
}
