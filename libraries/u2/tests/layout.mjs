// the tests' DOM lays nothing out: a test gives an element its sizes, as a browser test double does
export function layout(el, {rect, ...metrics}) {
  for (const [name, value] of Object.entries(metrics))
    Object.defineProperty(el, name, {value, configurable: true, writable: true});
  if (rect !== undefined) {
    el.getBoundingClientRect = () => rect;
    el.getClientRects = () => [rect];
  }
  return el;
}

if (ResizeObserver.instances === undefined) {
  const Base = ResizeObserver;
  globalThis.ResizeObserver = class ResizeObserver extends Base {
    static instances = [];

    constructor(callback) {
      super(callback);
      this.callback = callback;
      this.targets = [];
      this.disconnected = false;
      ResizeObserver.instances.push(this);
    }

    observe(target, options) {
      if (!this.targets.includes(target))
        this.targets.push(target);
      super.observe(target, options);
    }

    unobserve(target) {
      this.targets = this.targets.filter((t) => t !== target);
      super.unobserve(target);
    }

    disconnect() {
      this.disconnected = true;
      this.targets = [];
      super.disconnect();
    }
  };
}
