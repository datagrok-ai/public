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
  const session = feature(test);
  test("Picking a user and a group renders them as chips", {tag: ["@demo", "@realizes:u2.dg.entities"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Entities\" demo page", () => openDemoPage(page, "Entities"));
    enter(page, "U2 Demo");
    await test.step("Then Admin chip should be visible", () => shouldBe(page, el("Admin chip"), "visible"));
    await test.step("When user types \"adm\" into user picker", () => typeInto(page, "adm", el("user picker")));
    await test.step("And user selects \"Admin\" in user picker", () => selectIn(page, "Admin", el("user picker")));
    await test.step("Then value of user readout should have text \"Admin\"", () => shouldHaveText(page, el("value of user readout"), "Admin"));
    await test.step("When user types \"All users\" into group picker", () => typeInto(page, "All users", el("group picker")));
    await test.step("And user selects \"All users\" in group picker", () => selectIn(page, "All users", el("group picker")));
    await test.step("Then value of group readout should have text \"All users\"", () => shouldHaveText(page, el("value of group readout"), "All users"));
    await test.step("And \"All users\" chip should be visible", () => shouldBe(page, el("\"All users\" chip"), "visible"));
  });
});
