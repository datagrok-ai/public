interface GuardState {
  onWheel: (e: WheelEvent) => void;
  hint: HTMLElement;
  hideTimer: number;
}

const states = new WeakMap<HTMLElement, GuardState>();

/**
 * Keeps a zoomable viewer from hijacking wheel events inside a scrollable container:
 * plain wheel scrolls the page (events are stopped before the viewer sees them),
 * Ctrl/Cmd + wheel passes through to the viewer (zoom). A short hint is shown on plain wheel.
 */
export const wheelGuard = {
  mounted: (el: HTMLElement) => {
    if (getComputedStyle(el).position === 'static')
      el.style.position = 'relative';

    const hint = document.createElement('div');
    hint.textContent = 'Ctrl + scroll to zoom';
    hint.style.cssText = `
      position: absolute; left: 50%; top: 50%; transform: translate(-50%, -50%);
      background: rgba(60, 60, 60, 0.75); color: white; padding: 6px 14px; border-radius: 4px;
      font-size: 13px; pointer-events: none; z-index: 100; opacity: 0; transition: opacity 0.2s;`;
    el.appendChild(hint);

    const state: GuardState = {
      hint,
      hideTimer: 0,
      onWheel: (e: WheelEvent) => {
        if (e.ctrlKey || e.metaKey) {
          // keep the browser from zooming the page; the event still reaches the viewer
          e.preventDefault();
          return;
        }
        e.stopPropagation();
        hint.style.opacity = '1';
        clearTimeout(state.hideTimer);
        state.hideTimer = window.setTimeout(() => hint.style.opacity = '0', 800);
      },
    };
    el.addEventListener('wheel', state.onWheel, {capture: true, passive: false});
    states.set(el, state);
  },
  unmounted: (el: HTMLElement) => {
    const state = states.get(el);
    if (!state)
      return;
    el.removeEventListener('wheel', state.onWheel, {capture: true});
    clearTimeout(state.hideTimer);
    state.hint.remove();
    states.delete(el);
  },
};
