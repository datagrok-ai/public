/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/service/surface.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [helm.service-surface]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {helmInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {call, callWith, resultColumnLength, resultEveryStarts, resultHasMethods} from '@datagrok-libraries/bdd/bindings/platform/functions';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The service surface other packages call", () => {
  const session = feature(test, "features/service/surface.feature", import.meta.url);
  test("The helper exposes the methods other packages call", {tag: ["@realizes:helm.service-surface"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the Helm package is initialized", () => helmInitialized(page));
    await session.step(16, "When user calls \"Helm:getHelmHelper\" function", () => call(page, "Helm:getHelmHelper"));
    await session.step(17, "Then the result should have methods \"parse, removeGaps, getMolfiles, getHoveredAtom, createHelmInput, createHelmWebEditor, createWebEditorApp, overrideMonomersFuncs, revertOriginalMonomersFuncs, buildMonomersFuncsFromLib\"", () => resultHasMethods(page, "parse, removeGaps, getMolfiles, getHoveredAtom, createHelmInput, createHelmWebEditor, createWebEditorApp, overrideMonomersFuncs, revertOriginalMonomersFuncs, buildMonomersFuncsFromLib"));
    await session.step(18, "And no errors should have been logged", () => noErrors(page));
  });
  test("getMolfiles converts every row of a HELM column", {tag: ["@realizes:helm.service-surface"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the Helm package is initialized", () => helmInitialized(page));
    await session.step(21, "Given user opens helm-showcase dataset", () => openDataset(page, ds("helm-showcase")));
    await session.step(22, "When user calls \"Helm:getMolfiles\" function with:", () => callWith(page, "Helm:getMolfiles", [["col","column:HELM"]]), [["col","column:HELM"]]);
    await session.step(24, "Then the result should be a column of 53 values", () => resultColumnLength(page, 53));
    await session.step(25, "And every value of the result should start with \"HWE pseudo-molfile\"", () => resultEveryStarts(page, "HWE pseudo-molfile"));
    await session.step(26, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(27, "And no errors should have been logged", () => noErrors(page));
  });
});
