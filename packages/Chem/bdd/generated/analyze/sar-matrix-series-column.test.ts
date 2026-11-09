/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyze/sar-matrix-series-column.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.sar-matrix]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {matricesBuilt} from '../../bindings/sar-matrix.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("SAR Matrix series from the user's own Series column", () => {
  const session = feature(test, "features/analyze/sar-matrix-series-column.feature", import.meta.url);
  test("A Series column names one matrix per value", {tag: ["@realizes:chem.cp.sar-matrix"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "When user picks \"Chem > Analyze > SAR Matrix...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > SAR Matrix..."));
    await session.step(17, "Then \"SAR Matrix\" dialog should be visible", () => shouldBe(page, el("\"SAR Matrix\" dialog"), "visible"));
    await session.step(18, "When user selects \"Average Mass\" in Activity input in \"SAR Matrix\" dialog", () => selectIn(page, "Average Mass", el("Activity input in \"SAR Matrix\" dialog")));
    await session.step(19, "And user selects \"none\" in Scaling input in \"SAR Matrix\" dialog", () => selectIn(page, "none", el("Scaling input in \"SAR Matrix\" dialog")));
    await session.step(20, "And user checks \"Use existing R-groups\" input in \"SAR Matrix\" dialog", () => check(page, el("\"Use existing R-groups\" input in \"SAR Matrix\" dialog")));
    await session.step(21, "And user selects \"Core\" in Core input in \"SAR Matrix\" dialog", () => selectIn(page, "Core", el("Core input in \"SAR Matrix\" dialog")));
    await session.step(22, "And user clicks on editor of R-groups input in \"SAR Matrix\" dialog", () => clickOn(page, el("editor of R-groups input in \"SAR Matrix\" dialog")));
    await session.step(23, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(24, "When user clicks on None label in \"Select columns...\" dialog", () => clickOn(page, el("None label in \"Select columns...\" dialog")));
    await session.step(25, "And user clicks on the \"cell 2 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 2 of x", el("grid viewer in \"Select columns...\" dialog")));
    await session.step(26, "And user clicks on the \"cell 3 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 3 of x", el("grid viewer in \"Select columns...\" dialog")));
    await session.step(27, "And user clicks on the \"cell 4 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 4 of x", el("grid viewer in \"Select columns...\" dialog")));
    await session.step(28, "And user clicks on the \"cell 5 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 5 of x", el("grid viewer in \"Select columns...\" dialog")));
    await session.step(29, "And user clicks on the \"cell 6 of x\" area of grid viewer in \"Select columns...\" dialog", () => clickArea(page, "cell 6 of x", el("grid viewer in \"Select columns...\" dialog")));
    await session.step(30, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(31, "And user selects \"Series\" in Series input in \"SAR Matrix\" dialog", () => selectIn(page, "Series", el("Series input in \"SAR Matrix\" dialog")));
    await session.step(32, "Then \"SAR Matrix\" dialog should contain text \"One matrix per Series\"", () => shouldContainText(page, el("\"SAR Matrix\" dialog"), "One matrix per Series"));
    await session.step(33, "When user clicks on OK button in \"SAR Matrix\" dialog", () => clickOn(page, el("OK button in \"SAR Matrix\" dialog")));
    await session.step(34, "Then SAR Matrix Viewer viewer should have built its matrices", () => matricesBuilt(page, el("SAR Matrix Viewer viewer")));
    await session.step(35, "And the \"matrices\" reading of SAR Matrix Viewer viewer should be 1", () => readingIs(page, "matrices", el("SAR Matrix Viewer viewer"), 1));
    await session.step(36, "And the \"series\" reading of SAR Matrix Viewer viewer should be \"Triazoles\"", () => readingReads(page, "series", el("SAR Matrix Viewer viewer"), "Triazoles"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
});
