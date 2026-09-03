/* ---
generated: features/demo/platform/molecules.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.dg.molecules]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {enterInto, selectIn, shouldBe, shouldContainText, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Bridged chemistry inputs", () => {
  const session = feature(test);
  test("The structure input, the property form and the structure typeahead", {tag: ["@demo", "@realizes:u2.dg.molecules"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Molecules\" demo page", () => openDemoPage(page, "Molecules"));
    enter(page, "U2 Demo");
    await test.step("Then value of smiles readout should have text \"CC(=O)OC1=CC=CC=C1C(=O)O\"", () => shouldHaveText(page, el("value of smiles readout"), "CC(=O)OC1=CC=CC=C1C(=O)O"));
    await test.step("And Structure input should be visible", () => shouldBe(page, el("Structure input"), "visible"));
    await test.step("When user types \"Aspirin acetate\" into Name input in object form", () => typeInto(page, "Aspirin acetate", el("Name input in object form")));
    await test.step("Then value of compound readout should contain text \"Aspirin acetate\"", () => shouldContainText(page, el("value of compound readout"), "Aspirin acetate"));
    await test.step("When user enters \"181\" into \"MW, Da\" input in object form", () => enterInto(page, "181", el("\"MW, Da\" input in object form")));
    await test.step("Then value of compound readout should contain text \"\\\"mw\\\":181\"", () => shouldContainText(page, el("value of compound readout"), "\"mw\":181"));
    await test.step("When user types \"Caf\" into compound picker", () => typeInto(page, "Caf", el("compound picker")));
    await test.step("And user selects \"Caffeine\" in compound picker", () => selectIn(page, "Caffeine", el("compound picker")));
    await test.step("Then value of picked readout should have text \"Caffeine\"", () => shouldHaveText(page, el("value of picked readout"), "Caffeine"));
  });
});
