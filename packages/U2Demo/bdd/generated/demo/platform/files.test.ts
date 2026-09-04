/* ---
generated: features/demo/platform/files.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.dg.file-input]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {enterInto, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("File pickers", () => {
  const session = feature(test);
  test("A typed path resolves against the file shares", {tag: ["@demo", "@realizes:u2.dg.file-input"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Files\" demo page", () => openDemoPage(page, "Files"));
    enter(page, "U2 Demo");
    await test.step("Then value of file readout should have text \"(none)\"", () => shouldHaveText(page, el("value of file readout"), "(none)"));
    await test.step("When user enters \"System:DemoFiles/demog.csv\" into File input", () => enterInto(page, "System:DemoFiles/demog.csv", el("File input")));
    await test.step("Then value of file readout should have text \"demog.csv\"", () => shouldHaveText(page, el("value of file readout"), "demog.csv"));
  });
});
