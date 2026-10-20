/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/workspace.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openApp, rootFolder} from '../bindings/steps.js';
import {clickOn, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Parking orphan folders", () => {
  const session = feature(test, "features/workspace.feature", import.meta.url);
  test("An approver previews parking, parks a folder no program claims, and unparks it", {tag: ["@portal"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(8, "Given user opens the Project Portal app", () => openApp(page));
    enter(page, "Project Portal app");
    await session.step(11, "Given a folder \"BddParkScratch\" no program claims under the workspace root of the \"EX-SMD\" portfolio", () => rootFolder(page, "BddParkScratch", "EX-SMD"));
    await session.step(12, "When user clicks on small molecule discovery portfolio", () => clickOn(page, el("small molecule discovery portfolio")));
    await session.step(13, "And user clicks on workspace button", () => clickOn(page, el("workspace button")));
    await session.step(14, "Then workspace dialog should be visible", () => shouldBe(page, el("workspace dialog"), "visible"));
    await session.step(15, "And orphans should contain text \"BddParkScratch\"", () => shouldContainText(page, el("orphans"), "BddParkScratch"));
    await session.step(16, "And orphans should contain text \"EXS-ION4 is terminated\"", () => shouldContainText(page, el("orphans"), "EXS-ION4 is terminated"));
    await session.step(17, "When user clicks on park preview button", () => clickOn(page, el("park preview button")));
    await session.step(18, "Then park result should contain text \"would be parked in _Parked/\"", () => shouldContainText(page, el("park result"), "would be parked in _Parked/"));
    await session.step(19, "And park result should contain text \"nothing deleted\"", () => shouldContainText(page, el("park result"), "nothing deleted"));
    await session.step(20, "When user clicks on park scratch folder link", () => clickOn(page, el("park scratch folder link")));
    await session.step(21, "Then park result should contain text \"1 folder(s) parked\"", () => shouldContainText(page, el("park result"), "1 folder(s) parked"));
    await session.step(22, "And parked folders should contain text \"_Parked/\"", () => shouldContainText(page, el("parked folders"), "_Parked/"));
    await session.step(23, "And orphans should not contain text \"BddParkScratch\"", () => shouldNotContainText(page, el("orphans"), "BddParkScratch"));
    await session.step(24, "When user clicks on unpark scratch folder button", () => clickOn(page, el("unpark scratch folder button")));
    await session.step(25, "Then park result should contain text \"unparked\"", () => shouldContainText(page, el("park result"), "unparked"));
    await session.step(26, "And no parked folders should be visible", () => shouldBe(page, el("no parked folders"), "visible"));
    await session.step(27, "And orphans should contain text \"BddParkScratch\"", () => shouldContainText(page, el("orphans"), "BddParkScratch"));
  });
});
