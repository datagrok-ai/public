/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/display/messaging.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.message-input]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {pressKeyIn, shouldBe, shouldContainText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The message input", () => {
  const session = feature(test, "features/demo/display/messaging.feature", import.meta.url);
  test("Composing and sending", {tag: ["@demo", "@realizes:u2.message-input"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(5, "Given user opens the \"Message input\" demo page", () => openDemoPage(page, "Message input"));
    enter(page, "U2 Demo");
    await session.step(6, "Then Send button should be disabled", () => shouldBe(page, el("Send button"), "disabled"));
    await session.step(7, "When user types \"Hello from bdd\" into message input", () => typeInto(page, "Hello from bdd", el("message input")));
    await session.step(8, "Then Send button should be enabled", () => shouldBe(page, el("Send button"), "enabled"));
    await session.step(9, "When user presses Control+Enter in message input", () => pressKeyIn(page, "Control+Enter", el("message input")));
    await session.step(10, "Then value of sent readout should contain text \"Hello from bdd\"", () => shouldContainText(page, el("value of sent readout"), "Hello from bdd"));
  });
});
