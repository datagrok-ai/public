import * as Vue from 'vue';

export interface CompositorOverlayService {
  isActive: Vue.Ref<boolean>;
  show(): Promise<void>;
  hide(): void;
}

export function createCompositorOverlayService(): CompositorOverlayService {
  const isActive = Vue.ref(false);

  async function show() {
    isActive.value = true;
    await new Promise((r) => setTimeout(r, 0));
  }

  function hide() {
    isActive.value = false;
  }

  return {isActive, show, hide};
}
