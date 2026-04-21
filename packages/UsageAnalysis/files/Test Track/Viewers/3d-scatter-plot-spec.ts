import {test, expect} from '@playwright/test';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e.message ?? String(e)}); }
}

test('3D Scatter Plot tests', async ({page, baseURL}) => {
  test.setTimeout(300_000);

  await page.goto(baseURL ?? '/');
  await page.locator('[name="Toolbox"], [name="Browse"], .d4-sidebar').first().waitFor({timeout: 60000});
  await page.waitForFunction(() => {
    try {
      if (typeof grok === 'undefined' || !grok.shell) return false;
      grok.shell.settings.showFiltersIconsConstantly;
      return true;
    } catch { return false; }
  }, {timeout: 30000});

  // Setup
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Add 3D Scatter Plot
  await page.evaluate(async () => {
    const tv = grok.shell.tv;
    tv.addViewer('3d scatter plot');
    await new Promise(r => setTimeout(r, 1000));
  });

  // Helper: get 3d scatter plot viewer and setOptions
  const set3dOptions = (opts: Record<string, any>) =>
    page.evaluate(async (o: Record<string, any>) => {
      const tv = grok.shell.tv;
      const v3d = Array.from(tv.viewers).find((v: any) => v.type === '3d scatter plot') as any;
      v3d.setOptions(o);
      await new Promise(r => setTimeout(r, 400));
    }, opts);

  // ── Axis column assignment ──────────────────────────────────────────────────
  await softStep('S1: Set X=AGE, Y=HEIGHT, Z=WEIGHT', async () => {
    await set3dOptions({xColumnName: 'AGE', yColumnName: 'HEIGHT', zColumnName: 'WEIGHT'});
  });

  await softStep('S1: Set X=WEIGHT, Y=AGE, Z=HEIGHT', async () => {
    await set3dOptions({xColumnName: 'WEIGHT', yColumnName: 'AGE', zColumnName: 'HEIGHT'});
  });

  await softStep('S1: Set X back to AGE, Y=HEIGHT, Z=WEIGHT', async () => {
    await set3dOptions({xColumnName: 'AGE', yColumnName: 'HEIGHT', zColumnName: 'WEIGHT'});
  });

  // ── Axis types ──────────────────────────────────────────────────────────────
  await softStep('S2: Set X axis type to logarithmic', async () => {
    await set3dOptions({xAxisType: 'logarithmic'});
  });

  await softStep('S2: Set Y axis type to logarithmic', async () => {
    await set3dOptions({yAxisType: 'logarithmic'});
  });

  await softStep('S2: Set Z axis type to logarithmic', async () => {
    await set3dOptions({zAxisType: 'logarithmic'});
  });

  await softStep('S2: Set X axis type back to linear', async () => {
    await set3dOptions({xAxisType: 'linear'});
  });

  await softStep('S2: Set Y axis type back to linear', async () => {
    await set3dOptions({yAxisType: 'linear'});
  });

  await softStep('S2: Set Z axis type back to linear', async () => {
    await set3dOptions({zAxisType: 'linear'});
  });

  // ── Color coding — categorical ──────────────────────────────────────────────
  await softStep('S3: Set Color to SEX (categorical)', async () => {
    await set3dOptions({colorColumnName: 'SEX'});
  });

  await softStep('S3: Set Color to RACE', async () => {
    await set3dOptions({colorColumnName: 'RACE'});
  });

  await softStep('S3: Clear Color', async () => {
    await set3dOptions({colorColumnName: ''});
  });

  // ── Color coding — numerical ────────────────────────────────────────────────
  await softStep('S4: Set Color to AGE (numerical)', async () => {
    await set3dOptions({colorColumnName: 'AGE'});
  });

  await softStep('S4: Clear Color', async () => {
    await set3dOptions({colorColumnName: ''});
  });

  // ── Size coding ─────────────────────────────────────────────────────────────
  await softStep('S5: Set Size to WEIGHT', async () => {
    await set3dOptions({sizeColumnName: 'WEIGHT'});
  });

  await softStep('S5: Set Size to AGE', async () => {
    await set3dOptions({sizeColumnName: 'AGE'});
  });

  await softStep('S5: Clear Size', async () => {
    await set3dOptions({sizeColumnName: ''});
  });

  // ── Labels ──────────────────────────────────────────────────────────────────
  await softStep('S6: Set Label to SEX', async () => {
    await set3dOptions({labelColumnName: 'SEX'});
  });

  await softStep('S6: Clear Label', async () => {
    await set3dOptions({labelColumnName: ''});
  });

  // ── Marker type ─────────────────────────────────────────────────────────────
  await softStep('S7: Set Marker Type to sphere', async () => {
    await set3dOptions({markerType: 'sphere'});
  });

  await softStep('S7: Set Marker Type to box', async () => {
    await set3dOptions({markerType: 'box'});
  });

  await softStep('S7: Set Marker Type to cylinder', async () => {
    await set3dOptions({markerType: 'cylinder'});
  });

  await softStep('S7: Set Marker Type to tetrahedron', async () => {
    await set3dOptions({markerType: 'tetrahedron'});
  });

  await softStep('S7: Set Marker Type to dodecahedron', async () => {
    await set3dOptions({markerType: 'dodecahedron'});
  });

  await softStep('S7: Set Marker Type back to octahedron', async () => {
    await set3dOptions({markerType: 'octahedron'});
  });

  // ── Marker opacity and rotation ─────────────────────────────────────────────
  await softStep('S8: Set Marker Opacity to 20', async () => {
    await set3dOptions({markerOpacity: 20});
  });

  await softStep('S8: Set Marker Opacity to 100', async () => {
    await set3dOptions({markerOpacity: 100});
  });

  await softStep('S8: Set Marker Opacity back to 69', async () => {
    await set3dOptions({markerOpacity: 69});
  });

  await softStep('S8: Enable Marker Random Rotation', async () => {
    await set3dOptions({markerRandomRotation: true});
  });

  await softStep('S8: Disable Marker Random Rotation', async () => {
    await set3dOptions({markerRandomRotation: false});
  });

  // ── Filtered out points ─────────────────────────────────────────────────────
  await softStep('S9: Open filter panel and filter AGE to 20–40', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const fg = tv.getFiltersGroup();
      fg.updateOrAdd({type: 'histogram', column: 'AGE', min: 20, max: 40});
      await new Promise(r => setTimeout(r, 1000));
    });
  });

  await softStep('S9: Enable Show Filtered Out Points', async () => {
    await set3dOptions({showFilteredOutPoints: true});
  });

  await softStep('S9: Disable Show Filtered Out Points', async () => {
    await set3dOptions({showFilteredOutPoints: false});
  });

  await softStep('S9: Clear AGE filter', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      tv.dataFrame.filter.setAll(true);
      await new Promise(r => setTimeout(r, 500));
    });
  });

  // ── Axes visibility and grid lines ──────────────────────────────────────────
  await softStep('S10: Disable Show Axes', async () => {
    await set3dOptions({showAxes: false});
  });

  await softStep('S10: Enable Show Axes', async () => {
    await set3dOptions({showAxes: true});
  });

  await softStep('S10: Disable Vertical Grid Lines', async () => {
    await set3dOptions({verticalGridLines: false});
  });

  await softStep('S10: Disable Horizontal Grid Lines', async () => {
    await set3dOptions({horizontalGridLines: false});
  });

  await softStep('S10: Enable Vertical Grid Lines', async () => {
    await set3dOptions({verticalGridLines: true});
  });

  await softStep('S10: Enable Horizontal Grid Lines', async () => {
    await set3dOptions({horizontalGridLines: true});
  });

  // ── Background and colors (JS API only) ─────────────────────────────────────
  await softStep('S11: Set Back Color to black', async () => {
    await set3dOptions({backColor: 0xFF000000});
  });

  await softStep('S11: Set Axis Line Color to white', async () => {
    await set3dOptions({axisLineColor: 0xFFFFFFFF});
  });

  await softStep('S11: Restore Back Color and Axis Line Color to default', async () => {
    await set3dOptions({backColor: 0xFFFFFFFF, axisLineColor: 0xFF808080});
  });

  // ── Dynamic camera movement ──────────────────────────────────────────────────
  await softStep('S12: Enable Dynamic Camera Movement', async () => {
    await set3dOptions({dynamicCameraMovement: true});
    await page.waitForTimeout(1500);
  });

  await softStep('S12: Disable Dynamic Camera Movement', async () => {
    await set3dOptions({dynamicCameraMovement: false});
  });

  // ── Zoom and navigation ──────────────────────────────────────────────────────
  await softStep('S13: Scroll wheel up ×5 (zoom in)', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const v3d = Array.from(tv.viewers).find((v: any) => v.type === '3d scatter plot') as any;
      const canvas = v3d.root.querySelector('canvas');
      const rect = canvas.getBoundingClientRect();
      const cx = rect.left + rect.width / 2;
      const cy = rect.top + rect.height / 2;
      for (let i = 0; i < 5; i++) {
        canvas.dispatchEvent(new WheelEvent('wheel', {bubbles: true, cancelable: true, clientX: cx, clientY: cy, deltaY: -120}));
        await new Promise(r => setTimeout(r, 100));
      }
      await new Promise(r => setTimeout(r, 400));
    });
  });

  await softStep('S13: Scroll wheel down ×5 (zoom out)', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const v3d = Array.from(tv.viewers).find((v: any) => v.type === '3d scatter plot') as any;
      const canvas = v3d.root.querySelector('canvas');
      const rect = canvas.getBoundingClientRect();
      const cx = rect.left + rect.width / 2;
      const cy = rect.top + rect.height / 2;
      for (let i = 0; i < 5; i++) {
        canvas.dispatchEvent(new WheelEvent('wheel', {bubbles: true, cancelable: true, clientX: cx, clientY: cy, deltaY: 120}));
        await new Promise(r => setTimeout(r, 100));
      }
      await new Promise(r => setTimeout(r, 400));
    });
  });

  await softStep('S13: Right-click plot → Reset View', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const v3d = Array.from(tv.viewers).find((v: any) => v.type === '3d scatter plot') as any;
      const canvas = v3d.root.querySelector('canvas');
      const rect = canvas.getBoundingClientRect();
      const cx = rect.left + rect.width / 2;
      const cy = rect.top + rect.height / 2;
      canvas.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, clientX: cx, clientY: cy}));
      await new Promise(r => setTimeout(r, 500));
    });
    await page.locator('.d4-menu-item:has-text("Reset View")').click();
    await page.waitForTimeout(500);
  });

  // ── Mouse-over row group highlight ───────────────────────────────────────────
  await softStep('S14: Add Bar Chart viewer', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      tv.addViewer('Bar chart');
      await new Promise(r => setTimeout(r, 1000));
    });
  });

  await softStep('S14: Enable Show Mouse Over Row Group', async () => {
    await set3dOptions({showMouseOverRowGroup: true});
  });

  await softStep('S14: Hover over 3D plot center (mouse-over highlight)', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const v3d = Array.from(tv.viewers).find((v: any) => v.type === '3d scatter plot') as any;
      const canvas = v3d.root.querySelector('canvas');
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('mousemove', {bubbles: true,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2}));
      await new Promise(r => setTimeout(r, 500));
    });
  });

  await softStep('S14: Disable Show Mouse Over Row Group', async () => {
    await set3dOptions({showMouseOverRowGroup: false});
  });

  await softStep('S14: Close Bar Chart viewer', async () => {
    await page.evaluate(async () => {
      const tv = grok.shell.tv;
      const bc = Array.from(tv.viewers).find((v: any) => v.type === 'Bar chart') as any;
      if (bc) bc.close();
      await new Promise(r => setTimeout(r, 500));
    });
  });

  // Final check
  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `[${e.step}] ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
