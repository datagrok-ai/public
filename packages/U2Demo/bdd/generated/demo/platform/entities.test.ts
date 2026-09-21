/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/platform/entities.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.dg.entities]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {selectIn, shouldBe, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Entity pickers and chips", () => {
  const session = feature(test, "features/demo/platform/entities.feature", import.meta.url);
  test("Picking a user and a group renders them as chips", {tag: ["@demo", "@realizes:u2.dg.entities"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(6, "Given user opens the \"Entities\" demo page", () => openDemoPage(page, "Entities"));
    enter(page, "U2 Demo");
    await session.step(7, "Then Admin chip should be visible", () => shouldBe(page, el("Admin chip"), "visible"));
    await session.step(8, "When user types \"adm\" into user picker", () => typeInto(page, "adm", el("user picker")));
    await session.step(9, "And user selects \"Admin\" in user picker", () => selectIn(page, "Admin", el("user picker")));
    await session.step(10, "Then value of user readout should have text \"Admin\"", () => shouldHaveText(page, el("value of user readout"), "Admin"));
    await session.step(11, "When user types \"All users\" into group picker", () => typeInto(page, "All users", el("group picker")));
    await session.step(12, "And user selects \"All users\" in group picker", () => selectIn(page, "All users", el("group picker")));
    await session.step(13, "Then value of group readout should have text \"All users\"", () => shouldHaveText(page, el("value of group readout"), "All users"));
    await session.step(14, "And \"All users\" chip should be visible", () => shouldBe(page, el("\"All users\" chip"), "visible"));
  });
});
