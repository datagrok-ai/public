/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
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
  const session = feature(test, "features/demo/platform/molecules.feature", import.meta.url);
  test("The structure input, the property form and the structure typeahead", {tag: ["@demo", "@realizes:u2.dg.molecules"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user opens the \"Molecules\" demo page", () => openDemoPage(page, "Molecules"));
    enter(page, "U2 Demo");
    await session.step(8, "Then value of smiles readout should have text \"CC(=O)OC1=CC=CC=C1C(=O)O\"", () => shouldHaveText(page, el("value of smiles readout"), "CC(=O)OC1=CC=CC=C1C(=O)O"));
    await session.step(9, "And Structure input should be visible", () => shouldBe(page, el("Structure input"), "visible"));
    await session.step(10, "When user types \"Aspirin acetate\" into Name input in object form", () => typeInto(page, "Aspirin acetate", el("Name input in object form")));
    await session.step(11, "Then value of compound readout should contain text \"Aspirin acetate\"", () => shouldContainText(page, el("value of compound readout"), "Aspirin acetate"));
    await session.step(12, "When user enters \"181\" into \"MW, Da\" input in object form", () => enterInto(page, "181", el("\"MW, Da\" input in object form")));
    await session.step(13, "Then value of compound readout should contain text \"\\\"mw\\\":181\"", () => shouldContainText(page, el("value of compound readout"), "\"mw\":181"));
    await session.step(14, "When user types \"Caf\" into compound picker", () => typeInto(page, "Caf", el("compound picker")));
    await session.step(15, "And user selects \"Caffeine\" in compound picker", () => selectIn(page, "Caffeine", el("compound picker")));
    await session.step(16, "Then value of picked readout should have text \"Caffeine\"", () => shouldHaveText(page, el("value of picked readout"), "Caffeine"));
  });
});
