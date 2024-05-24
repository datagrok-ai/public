
import { defineComponent } from 'vue';
import type { FormT } from '@datagrok-libraries/webcomponents';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-form': FormT
    }
  }
}

export const Form = defineComponent({
  name: 'Form',
  setup() {
    return () => (
      <keep-alive>
        <dg-form></dg-form>
      </keep-alive>
    );
  }
});
