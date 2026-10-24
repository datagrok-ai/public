/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/mpo/data-driven-training-data.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, close, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelOpen, openDataset, packageInstalled, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A data-driven profile retrains into its full layout after its Training Data panel was closed", () => {
  const session = feature(test, "features/mpo/data-driven-training-data.feature", import.meta.url);
  test("Picking another dataset after closing Training Data brings the panel back", async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(11, "And the \"EDA\" package is installed", () => packageInstalled(page, "EDA"));
    await session.step(12, "And user opens drugs-props-train dataset", () => openDataset(page, ds("drugs-props-train")));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(15, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(16, "And Apps---Chem tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Chem tree node inside browse tree")));
    await session.step(19, "When user clicks on Apps---Chem---MPO-profiles tree node inside browse tree", () => clickOn(page, el("Apps---Chem---MPO-profiles tree node inside browse tree")));
    await session.step(20, "And user clicks on \"Create profile\" button", () => clickOn(page, el("\"Create profile\" button")));
    await session.step(21, "Then the \"Untitled Profile\" view should be current", () => viewIsCurrent(page, "Untitled Profile"));
    await session.step(22, "When user selects \"Data-driven\" in Method input", () => selectIn(page, "Data-driven", el("Method input")));
    await session.step(23, "And user selects \"drugs-props-train\" in Dataset input", () => selectIn(page, "drugs-props-train", el("Dataset input")));
    await session.step(24, "Then \"Training Data\" tab should be visible", () => shouldBe(page, el("\"Training Data\" tab"), "visible"));
    await session.step(25, "When user closes \"Training Data\" tab", () => close(page, el("\"Training Data\" tab")));
    await session.step(26, "Then \"Training Data\" tab should be absent", () => shouldBe(page, el("\"Training Data\" tab"), "absent"));
    await session.step(27, "When user selects \"demog-1000\" in Dataset input", () => selectIn(page, "demog-1000", el("Dataset input")));
    await session.step(28, "Then \"Training Data\" tab should be visible", () => shouldBe(page, el("\"Training Data\" tab"), "visible"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
  });
});
