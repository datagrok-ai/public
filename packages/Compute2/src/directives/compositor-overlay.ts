import * as ui from 'datagrok-api/ui';
import * as Vue from 'vue';
import './compositor-overlay.css';

interface OverlayState {
  overlay: HTMLElement;
  isActive: boolean;
}

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

export const compositorOverlay: Vue.Directive<HTMLElement, boolean> = {
  mounted(el, binding) {
    const overlay = ui.div([
      ui.loader(),
    ], 'c2-compositor-overlay');

    const position = getComputedStyle(el).position;
    if (position === 'static' || position === '')
      el.style.position = 'relative';

    el.appendChild(overlay);

    const state: OverlayState = {overlay, isActive: false};
    stateMap.set(el, state);

    if (binding.value)
      show(state, el);
  },

  updated(el, binding) {
    const state = stateMap.get(el);
    if (!state) return;

    if (binding.value)
      show(state, el);
    else
      hide(state);
  },

  beforeUnmount(el) {
    const state = stateMap.get(el);
    if (state) {
      state.overlay.remove();
      stateMap.delete(el);
    }
  },
};
