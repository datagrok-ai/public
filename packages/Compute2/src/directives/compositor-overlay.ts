import * as ui from 'datagrok-api/ui';
import * as Vue from 'vue';
import './compositor-overlay.css';

interface OverlayState {
  overlay: HTMLElement;
  isActive: boolean;
  allowInput: boolean;
  swallow: (e: Event) => void;
}

// The overlay stops the pointer but not the keyboard; swallowing also keeps the focus where it is,
// which a mousedown on the overlay would otherwise move.
const SWALLOWED_EVENTS = [
  'pointerdown', 'mousedown', 'click', 'dblclick', 'auxclick', 'contextmenu',
  'keydown', 'keypress', 'beforeinput',
];

const stateMap = new WeakMap<HTMLElement, OverlayState>();

function show(state: OverlayState, el: HTMLElement): void {
  if (state.isActive) return;
  state.isActive = true;

  const overlay = state.overlay;
  overlay.style.opacity = '';
  overlay.classList.add('c2-compositor-overlay--active');

  // Force synchronous layout commit so the compositor has the element's
  // geometry before heavy JS starts. Without this, the browser may batch
  // the DOM mutation with the heavy work and never paint the overlay.
  el.getBoundingClientRect();
}

function hide(state: OverlayState): void {
  if (!state.isActive) return;
  state.isActive = false;

  const overlay = state.overlay;
  overlay.classList.remove('c2-compositor-overlay--active');
  overlay.style.opacity = '0';
}

/** `v-compositor-overlay.allowInput` lets the covered inputs keep taking input while active. */
export const compositorOverlay: Vue.Directive<HTMLElement, boolean, 'allowInput'> = {
  mounted(el, binding) {
    const overlay = ui.div([
      ui.loader(),
    ], 'c2-compositor-overlay');

    const position = getComputedStyle(el).position;
    if (position === 'static' || position === '')
      el.style.position = 'relative';

    el.appendChild(overlay);

    const swallow = (e: Event) => {
      if (state.isActive && !state.allowInput) {
        e.preventDefault();
        e.stopImmediatePropagation();
      }
    };
    const state: OverlayState = {overlay, isActive: false, allowInput: !!binding.modifiers.allowInput, swallow};
    stateMap.set(el, state);
    for (const type of SWALLOWED_EVENTS)
      el.addEventListener(type, swallow, true);

    if (binding.value)
      show(state, el);
  },

  updated(el, binding) {
    const state = stateMap.get(el);
    if (!state) return;

    state.allowInput = !!binding.modifiers.allowInput;
    if (binding.value)
      show(state, el);
    else
      hide(state);
  },

  beforeUnmount(el) {
    const state = stateMap.get(el);
    if (state) {
      for (const type of SWALLOWED_EVENTS)
        el.removeEventListener(type, state.swallow, true);
      state.overlay.remove();
      stateMap.delete(el);
    }
  },
};
