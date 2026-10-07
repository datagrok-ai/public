import {Viewer} from './src/viewer';

/** A div with the children. */
export function div(children: HTMLElement[] = [], className?: string): HTMLDivElement {
  return document.createElement('div');
}

export function button(text: string): HTMLButtonElement;
export function button(text: string, onClick: () => void): HTMLButtonElement;
export function button(text: string, onClick?: () => void): HTMLButtonElement {
  return document.createElement('button');
}

export function longSignature(aVeryLongParameterName: string, anotherVeryLongParameterName: number, yetAnotherVeryLongParameterName: boolean,
  oneMoreVeryLongParameterName: string[], theLastVeryLongParameterName: Record<string, string>, viewer: Viewer): void {
}

export namespace input {
  export function string(name: string, value?: string): HTMLInputElement {
    return document.createElement('input');
  }
}

function internal(): void {
}
