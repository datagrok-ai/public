import {Control, DateInput, Form, TextInput, button, component, notify} from '@datagrok-libraries/u2';

const EMAIL = /^[^@\s]+@[^@\s]+\.[^@\s]+$/;

export function bookingForm(): Control {
  return component(() => {
    const name = new TextInput({label: 'Name', required: true});
    const email = new TextInput({label: 'Email', required: true, placeholder: 'name@company.com'});
    email.addValidator((v) => EMAIL.test(v) ? null : 'Not an email address');
    const from = new DateInput({label: 'From', required: true});
    const to = new DateInput({label: 'To', required: true});
    from.addValidator((d) => d && to.value && d > to.value ? 'Must be before To' : null);
    const form = new Form().addAll([name, email, from, to]);
    return form.addButtons((row) => row.append(
      button('Save', () => notify.info(`Saved ${name.value}`), {primary: true, disabled: form.problem}),
      button('Check', () => form.validate()),
    ));
  }, {name: 'booking-form'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = bookingForm();
  host.append(demo.root);
  return demo;
}
