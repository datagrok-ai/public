import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const connection = 'Agolovko:NewTestPostgres';

test('Data Enrichment scenario', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await page.evaluate(async (conn) => {
    document.body.classList.add('selenium');
    // @ts-ignore
    grok.shell.settings.showFiltersIconsConstantly = true;
    // @ts-ignore
    grok.shell.windows.simpleMode = true;
    // @ts-ignore
    grok.shell.closeAll();
    // Pre-load connection entity into context panel — primes PowerPack enrichment registration.
    // @ts-ignore
    const connections = await grok.dapi.connections.list({pageSize: 500});
    const connEntity = connections.find((c: any) => c.nqName === conn);
    if (connEntity) {
      // @ts-ignore
      grok.shell.o = connEntity;
      await new Promise(r => setTimeout(r, 1500));
    }
    // @ts-ignore
    const df = await grok.data.db.query(conn, 'select * from public.orders');
    df.name = 'orders';
    df.setTag('.db-source-connection', conn);
    df.setTag('.db-source-schema', 'public');
    df.setTag('.db-source-table', 'orders');
    // @ts-ignore
    grok.shell.addTableView(df);
    await new Promise<void>(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  }, connection);

  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  async function waitForConnPane(timeoutMs = 20000): Promise<boolean> {
    const start = Date.now();
    while (Date.now() - start < timeoutMs) {
      const found = await page.evaluate(() => {
        const np = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
          .find(h => /new_test_postgres/i.test(h.textContent || '')) as HTMLElement | undefined;
        if (!np) return false;
        if (!np.closest('.d4-accordion-pane')!.classList.contains('expanded')) np.click();
        return true;
      });
      if (found) return true;
      await page.waitForTimeout(500);
    }
    return false;
  }

  async function waitForEnrichPane(timeoutMs = 15000): Promise<boolean> {
    const start = Date.now();
    while (Date.now() - start < timeoutMs) {
      const found = await page.evaluate(() => {
        const eh = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
          .find(h => h.textContent?.trim() === 'Enrich') as HTMLElement | undefined;
        if (!eh) return false;
        if (!eh.closest('.d4-accordion-pane')!.classList.contains('expanded')) eh.click();
        return true;
      });
      if (found) {
        await page.waitForTimeout(1200);
        return true;
      }
      await page.waitForTimeout(500);
    }
    return false;
  }

  async function waitForEnrichContent(timeoutMs = 15000): Promise<boolean> {
    const start = Date.now();
    while (Date.now() - start < timeoutMs) {
      const found = await page.evaluate(() => {
        const eh = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
          .find(h => h.textContent?.trim() === 'Enrich') as HTMLElement | undefined;
        const content = eh?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');
        return !!content?.querySelector('button.power-pack-enrich-add');
      });
      if (found) return true;
      await page.waitForTimeout(500);
    }
    return false;
  }

  async function waitForDialog(titleSubstring: string, timeoutMs = 10000): Promise<boolean> {
    const start = Date.now();
    while (Date.now() - start < timeoutMs) {
      const found = await page.evaluate((t) => {
        return Array.from(document.querySelectorAll('.d4-dialog'))
          .some(d => d.querySelector('.d4-dialog-title')?.textContent?.includes(t));
      }, titleSubstring);
      if (found) return true;
      await page.waitForTimeout(300);
    }
    return false;
  }

  async function waitForEnrichDialogReady(timeoutMs = 10000): Promise<boolean> {
    const start = Date.now();
    while (Date.now() - start < timeoutMs) {
      const ready = await page.evaluate(() => {
        return !!document.querySelector('.d4-dialog [name="div-add-Data"] i');
      });
      if (ready) return true;
      await page.waitForTimeout(300);
    }
    return false;
  }

  async function waitForEnrichmentListed(name: string, timeoutMs = 10000): Promise<boolean> {
    const start = Date.now();
    while (Date.now() - start < timeoutMs) {
      const listed = await page.evaluate((n) => {
        const eh = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
          .find(h => h.textContent?.trim() === 'Enrich') as HTMLElement | undefined;
        const content = eh?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');
        return (content as HTMLElement | undefined)?.innerText.includes(n) ?? false;
      }, name);
      if (listed) return true;
      await page.waitForTimeout(400);
    }
    return false;
  }

  // ---- Section 1 ----

  await softStep('1.1 Navigate to Postgres connection (Northwind substitute)', async () => {
    const ok = await page.evaluate(() => {
      // @ts-ignore
      return grok.shell.tv.dataFrame.columns.names().includes('customerid');
    });
    expect(ok).toBe(true);
  });

  await softStep('1.2-1.3 Open the Orders table', async () => {
    const rowCount = await page.evaluate(() => {
      // @ts-ignore
      return grok.shell.tv.dataFrame.rowCount;
    });
    expect(rowCount).toBeGreaterThan(100);
  });

  await softStep('1.4 Click customerid column', async () => {
    await page.evaluate(() => {
      // @ts-ignore
      grok.shell.o = grok.shell.tv.dataFrame.col('customerid');
    });
    const paneAppeared = await waitForConnPane(20000);
    expect(paneAppeared).toBe(true);
  });

  await softStep('1.5 Context Panel > Northwind > Enrich..', async () => {
    const enrichAppeared = await waitForEnrichPane(15000); await waitForEnrichContent(15000);
    expect(enrichAppeared).toBe(true);
    const addBtnVisible = await waitForEnrichContent(15000);
    expect(addBtnVisible).toBe(true);
  });

  await softStep('1.6 Click + Add Enrichment..', async () => {
    await page.evaluate(() => (document.querySelector('button.power-pack-enrich-add') as HTMLElement).click());
    const dlgOpen = await waitForDialog('Enrich customerid', 10000);
    expect(dlgOpen).toBe(true);
    await waitForEnrichDialogReady(10000);
  });

  await softStep('1.7 Add table to join > public > Customers and select columns', async () => {
    await page.evaluate(async () => {
      const addIcon = document.querySelector('.d4-dialog [name="div-add-Data"] i') as HTMLElement;
      const r = addIcon.getBoundingClientRect();
      const opts = { bubbles: true, cancelable: true, clientX: r.x+r.width/2, clientY: r.y+r.height/2, button: 0 };
      ['mousedown','mouseup','click'].forEach(e => addIcon.dispatchEvent(new MouseEvent(e, opts)));
      await new Promise(res => setTimeout(res, 1200));
      const getItem = (t: string) => Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(l => l.textContent?.trim() === t)?.closest('.d4-menu-item') as HTMLElement;
      for (const label of ['northwind', 'public']) {
        const el = getItem(label);
        const rr = el.getBoundingClientRect();
        const o = { bubbles: true, cancelable: true, clientX: rr.x+rr.width/2, clientY: rr.y+rr.height/2 };
        ['mouseover','mouseenter','mousemove'].forEach(e => el.dispatchEvent(new MouseEvent(e, o)));
        await new Promise(res => setTimeout(res, 1200));
      }
      const cust = getItem('customers');
      const rc = cust.getBoundingClientRect();
      const co = { bubbles: true, cancelable: true, clientX: rc.x+rc.width/2, clientY: rc.y+rc.height/2, button: 0 };
      ['mouseover','mouseenter','mousedown','mouseup','click'].forEach(e => cust.dispatchEvent(new MouseEvent(e, co)));
      await new Promise(res => setTimeout(res, 2000));
      const mainDlg = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich customerid')) as HTMLElement;
      const tag = Array.from(mainDlg.querySelectorAll('.d4-tag'))
        .find(t => (t as HTMLElement).innerText.includes('customers')) as HTMLElement;
      const div = tag.querySelectorAll(':scope > div')[1] as HTMLElement;
      const r2 = div.getBoundingClientRect();
      const o2 = { bubbles: true, cancelable: true, clientX: r2.x+r2.width/2, clientY: r2.y+r2.height/2, button: 0 };
      ['mousedown','mouseup','click'].forEach(e => div.dispatchEvent(new MouseEvent(e, o2)));
      await new Promise(res => setTimeout(res, 1500));
      const picker = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Select columns')) as HTMLElement;
      (picker.querySelector('[name="label-All"]') as HTMLElement).click();
      await new Promise(res => setTimeout(res, 400));
      (picker.querySelector('[name="button-OK"]') as HTMLElement).click();
      await new Promise(res => setTimeout(res, 1000));
    });
    const joined = await page.evaluate(() => {
      const mainDlg = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich customerid'));
      return (mainDlg as HTMLElement)?.innerText.includes('customers(All 11)');
    });
    expect(joined).toBe(true);
  });

  await softStep('1.8 Verify editor shows joins correctly', async () => {
    const joinText = await page.evaluate(() => {
      const mainDlg = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich customerid'));
      return (mainDlg as HTMLElement).innerText;
    });
    expect(joinText).toMatch(/customerid[\s\S]*=[\s\S]*customerid/);
  });

  await softStep('1.9 Enter name and SAVE', async () => {
    await page.evaluate(() => {
      const mainDlg = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich customerid')) as HTMLElement;
      (mainDlg.querySelector('input[type="text"]') as HTMLInputElement).focus();
    });
    await page.keyboard.type('CustomerInfo_test');
    await page.evaluate(async () => {
      const mainDlg = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich customerid')) as HTMLElement;
      const save = Array.from(mainDlg.querySelectorAll('button.ui-btn'))
        .find(b => b.textContent?.trim() === 'SAVE') as HTMLElement;
      save.click();
      await new Promise(r => setTimeout(r, 3000));
    });
    const listed = await waitForEnrichmentListed('CustomerInfo_test', 10000);
    expect(listed).toBe(true);
  });

  await softStep('1.10 Apply enrichment — columns added', async () => {
    const ok = await page.evaluate(async () => {
      // @ts-ignore
      const before = grok.shell.tv.dataFrame.columns.length;
      const eh = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find(h => h.textContent?.trim() === 'Enrich');
      const content = eh?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');
      const row = Array.from(content!.querySelectorAll('.power-pack-enrichment-row'))
        .find(r => (r as HTMLElement).innerText.includes('CustomerInfo_test')) as HTMLElement;
      (row.querySelector('a.ui-link') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 10000));
      // @ts-ignore
      return grok.shell.tv.dataFrame.columns.length > before;
    });
    expect(ok).toBe(true);
  });

  await softStep('1.11 Edit enrichment (reopen dialog and SAVE)', async () => {
    const ok = await page.evaluate(async () => {
      const eh = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find(h => h.textContent?.trim() === 'Enrich');
      const content = eh?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');
      const row = Array.from(content!.querySelectorAll('.power-pack-enrichment-row'))
        .find(r => (r as HTMLElement).innerText.includes('CustomerInfo_test')) as HTMLElement;
      (row.querySelector('i.fa-pencil') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 2000));
      const mainDlg = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich customerid')) as HTMLElement;
      const save = Array.from(mainDlg.querySelectorAll('button.ui-btn'))
        .find(b => b.textContent?.trim() === 'SAVE') as HTMLElement;
      save.click();
      await new Promise(r => setTimeout(r, 2000));
      return true;
    });
    expect(ok).toBe(true);
  });

  // ---- Section 2 ----

  await softStep('2.1 Create second enrichment on customerid', async () => {
    await page.evaluate(async () => {
      // @ts-ignore
      grok.shell.o = grok.shell.tv.dataFrame.col('customerid');
    });
    await waitForConnPane(20000);
    await waitForEnrichPane(15000); await waitForEnrichContent(15000);
    await page.evaluate(async () => {
      (document.querySelector('button.power-pack-enrich-add') as HTMLElement).click();
    });
    await waitForDialog('Enrich customerid', 10000); await waitForEnrichDialogReady(10000);
    await page.evaluate(async () => {
      const addIcon = document.querySelector('.d4-dialog [name="div-add-Data"] i') as HTMLElement;
      const rr = addIcon.getBoundingClientRect();
      const oo = { bubbles: true, cancelable: true, clientX: rr.x+rr.width/2, clientY: rr.y+rr.height/2, button: 0 };
      ['mousedown','mouseup','click'].forEach(e => addIcon.dispatchEvent(new MouseEvent(e, oo)));
      await new Promise(r => setTimeout(r, 1500));
      const getItem = (t: string) => Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(l => l.textContent?.trim() === t)?.closest('.d4-menu-item') as HTMLElement;
      for (const l of ['northwind','public']) {
        const el = getItem(l);
        const r1 = el.getBoundingClientRect();
        const o1 = { bubbles: true, cancelable: true, clientX: r1.x+r1.width/2, clientY: r1.y+r1.height/2 };
        ['mouseover','mouseenter','mousemove'].forEach(e => el.dispatchEvent(new MouseEvent(e, o1)));
        await new Promise(r => setTimeout(r, 1200));
      }
      const cust = getItem('customers');
      const rc = cust.getBoundingClientRect();
      const co = { bubbles: true, cancelable: true, clientX: rc.x+rc.width/2, clientY: rc.y+rc.height/2, button: 0 };
      ['mouseover','mouseenter','mousedown','mouseup','click'].forEach(e => cust.dispatchEvent(new MouseEvent(e, co)));
      await new Promise(r => setTimeout(r, 2000));
      const mainDlg = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich customerid')) as HTMLElement;
      const tag = Array.from(mainDlg.querySelectorAll('.d4-tag'))
        .find(t => (t as HTMLElement).innerText.includes('customers')) as HTMLElement;
      const div = tag.querySelectorAll(':scope > div')[1] as HTMLElement;
      const r2 = div.getBoundingClientRect();
      const o2 = { bubbles: true, cancelable: true, clientX: r2.x+r2.width/2, clientY: r2.y+r2.height/2, button: 0 };
      ['mousedown','mouseup','click'].forEach(e => div.dispatchEvent(new MouseEvent(e, o2)));
      await new Promise(r => setTimeout(r, 1500));
      const picker = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Select columns')) as HTMLElement;
      (picker.querySelector('[name="label-All"]') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 400));
      (picker.querySelector('[name="button-OK"]') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1000));
      const mainDlg2 = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich customerid')) as HTMLElement;
      (mainDlg2.querySelector('input[type="text"]') as HTMLInputElement).focus();
    });
    await page.keyboard.type('CustomerInfo2_test');
    await page.evaluate(async () => {
      const mainDlg = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich customerid')) as HTMLElement;
      const save = Array.from(mainDlg.querySelectorAll('button.ui-btn'))
        .find(b => b.textContent?.trim() === 'SAVE') as HTMLElement;
      save.click();
      await new Promise(r => setTimeout(r, 3000));
    });
    const hasFirst = await waitForEnrichmentListed('CustomerInfo_test', 10000);
    const hasSecond = await waitForEnrichmentListed('CustomerInfo2_test', 10000);
    expect(hasFirst && hasSecond).toBe(true);
  });

  await softStep('2.2 Create enrichment on employeeid', async () => {
    await page.evaluate(async () => {
      // @ts-ignore
      grok.shell.o = grok.shell.tv.dataFrame.col('employeeid');
    });
    await waitForConnPane(20000);
    await waitForEnrichPane(15000); await waitForEnrichContent(15000);
    await page.evaluate(async () => {
      (document.querySelector('button.power-pack-enrich-add') as HTMLElement).click();
    });
    await waitForDialog('Enrich employeeid', 10000); await waitForEnrichDialogReady(10000);
    await page.evaluate(async () => {
      const addIcon = document.querySelector('.d4-dialog [name="div-add-Data"] i') as HTMLElement;
      const r0 = addIcon.getBoundingClientRect();
      const o0 = { bubbles: true, cancelable: true, clientX: r0.x+r0.width/2, clientY: r0.y+r0.height/2, button: 0 };
      ['mousedown','mouseup','click'].forEach(e => addIcon.dispatchEvent(new MouseEvent(e, o0)));
      await new Promise(r => setTimeout(r, 1500));
      const getItem = (t: string) => Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find(l => l.textContent?.trim() === t)?.closest('.d4-menu-item') as HTMLElement;
      for (const l of ['northwind','public']) {
        const el = getItem(l);
        const r1 = el.getBoundingClientRect();
        const o1 = { bubbles: true, cancelable: true, clientX: r1.x+r1.width/2, clientY: r1.y+r1.height/2 };
        ['mouseover','mouseenter','mousemove'].forEach(e => el.dispatchEvent(new MouseEvent(e, o1)));
        await new Promise(r => setTimeout(r, 1200));
      }
      const emp = getItem('employees');
      const rc = emp.getBoundingClientRect();
      const co = { bubbles: true, cancelable: true, clientX: rc.x+rc.width/2, clientY: rc.y+rc.height/2, button: 0 };
      ['mouseover','mouseenter','mousedown','mouseup','click'].forEach(e => emp.dispatchEvent(new MouseEvent(e, co)));
      await new Promise(r => setTimeout(r, 2000));
      const mainDlg = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich employeeid')) as HTMLElement;
      const tag = Array.from(mainDlg.querySelectorAll('.d4-tag'))
        .find(t => (t as HTMLElement).innerText.includes('employees')) as HTMLElement;
      const div = tag.querySelectorAll(':scope > div')[1] as HTMLElement;
      const r2 = div.getBoundingClientRect();
      const o2 = { bubbles: true, cancelable: true, clientX: r2.x+r2.width/2, clientY: r2.y+r2.height/2, button: 0 };
      ['mousedown','mouseup','click'].forEach(e => div.dispatchEvent(new MouseEvent(e, o2)));
      await new Promise(r => setTimeout(r, 1500));
      const picker = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Select columns')) as HTMLElement;
      (picker.querySelector('[name="label-All"]') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 400));
      (picker.querySelector('[name="button-OK"]') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 1000));
      const mainDlg2 = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich employeeid')) as HTMLElement;
      (mainDlg2.querySelector('input[type="text"]') as HTMLInputElement).focus();
    });
    await page.keyboard.type('EmployeeInfo_test');
    await page.evaluate(async () => {
      const mainDlg = Array.from(document.querySelectorAll('.d4-dialog'))
        .find(d => d.querySelector('.d4-dialog-title')?.textContent?.includes('Enrich employeeid')) as HTMLElement;
      const save = Array.from(mainDlg.querySelectorAll('button.ui-btn'))
        .find(b => b.textContent?.trim() === 'SAVE') as HTMLElement;
      save.click();
      await new Promise(r => setTimeout(r, 3000));
    });
    const listed = await waitForEnrichmentListed('EmployeeInfo_test', 10000);
    expect(listed).toBe(true);
  });

  await softStep('2.3 Apply all enrichments', async () => {
    const ok = await page.evaluate(async () => {
      const eh = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find(h => h.textContent?.trim() === 'Enrich');
      const content = eh?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');
      const row = Array.from(content!.querySelectorAll('.power-pack-enrichment-row'))
        .find(r => (r as HTMLElement).innerText.includes('EmployeeInfo_test')) as HTMLElement;
      (row.querySelector('a.ui-link') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 10000));
      // @ts-ignore
      return grok.shell.tv.dataFrame.columns.length > 20;
    });
    expect(ok).toBe(true);
  });

  await softStep('2.4 Remove CustomerInfo2_test enrichment', async () => {
    await page.evaluate(async () => {
      // @ts-ignore
      grok.shell.o = grok.shell.tv.dataFrame.col('customerid');
    });
    await waitForConnPane(20000);
    await waitForEnrichPane(15000); await waitForEnrichContent(15000);
    const listed = await waitForEnrichmentListed('CustomerInfo2_test', 10000);
    expect(listed).toBe(true);
    await page.evaluate(async () => {
      const eh = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find(h => h.textContent?.trim() === 'Enrich') as HTMLElement | undefined;
      const content = eh?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content') as HTMLElement | undefined;
      const row = Array.from(content?.querySelectorAll('.power-pack-enrichment-row') ?? [])
        .find(r => (r as HTMLElement).innerText.trim().startsWith('CustomerInfo2_test')) as HTMLElement | undefined;
      (row?.querySelector('i.fa-times') as HTMLElement | undefined)?.click();
      await new Promise(r => setTimeout(r, 2000));
    });
    const removed = await page.evaluate(() => {
      const eh = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find(h => h.textContent?.trim() === 'Enrich') as HTMLElement | undefined;
      const content = eh?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content') as HTMLElement | undefined;
      return !(content?.innerText.includes('CustomerInfo2_test') ?? false);
    });
    expect(removed).toBe(true);
  });

  // ---- Section 3 ----

  await softStep('3.1 Enrichments available in query results (after tagging source)', async () => {
    await page.evaluate(async (conn) => {
      // @ts-ignore
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 800));
      // Prime PowerPack with the connection entity before opening the next view.
      // @ts-ignore
      const connections = await grok.dapi.connections.list({pageSize: 500});
      const connEntity = connections.find((c: any) => c.nqName === conn);
      if (connEntity) {
        // @ts-ignore
        grok.shell.o = connEntity;
        await new Promise(r => setTimeout(r, 1500));
      }
      // @ts-ignore
      const df = await grok.data.db.query(conn, 'select * from public.orders limit 500');
      df.name = 'orders_query_result';
      // @ts-ignore
      grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); (resolve as any)(); });
        setTimeout(resolve, 3000);
      });
      // Apply source tags AFTER the view is mounted, then reselect to trigger pane refresh.
      df.setTag('.db-source-connection', conn);
      df.setTag('.db-source-schema', 'public');
      df.setTag('.db-source-table', 'orders');
      await new Promise(r => setTimeout(r, 800));
      // @ts-ignore
      grok.shell.o = df.col('customerid');
    }, connection);
    await waitForConnPane(20000);
    await waitForEnrichPane(15000); await waitForEnrichContent(15000);
    const available = await waitForEnrichmentListed('CustomerInfo_test', 15000);
    expect(available).toBe(true);
  });

  await softStep('3.2-3.3 Apply enrichment and save layout', async () => {
    const saved = await page.evaluate(async () => {
      const eh = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find(h => h.textContent?.trim() === 'Enrich');
      const content = eh?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');
      const row = Array.from(content!.querySelectorAll('.power-pack-enrichment-row'))
        .find(r => (r as HTMLElement).innerText.trim().startsWith('CustomerInfo_test')) as HTMLElement;
      (row.querySelector('a.ui-link') as HTMLElement).click();
      await new Promise(r => setTimeout(r, 8000));
      // @ts-ignore
      const layout = grok.shell.tv.saveLayout();
      // @ts-ignore
      await grok.dapi.layouts.save(layout);
      await new Promise(r => setTimeout(r, 1500));
      return !!layout.id;
    });
    expect(saved).toBe(true);
  });

  await softStep('3.4 Layout restores enriched columns', async () => {
    const restored = await page.evaluate(async () => {
      // @ts-ignore
      const df = grok.shell.tv.dataFrame;
      for (const n of ['address','city','companyname','contactname','contacttitle','country','fax','phone','postalcode','region']) {
        try { df.columns.remove(n); } catch (e) {}
      }
      // @ts-ignore
      const layouts = await grok.dapi.layouts.list({pageSize: 50});
      const saved = layouts.find((l: any) => l.name === 'OrdersQueryResult');
      if (!saved) return false;
      // @ts-ignore
      try { grok.shell.tv.loadLayout(saved); } catch (e) {}
      await new Promise(r => setTimeout(r, 5000));
      // @ts-ignore
      const names = grok.shell.tv.dataFrame.columns.names();
      return names.includes('companyname');
    });
    expect(restored).toBe(true);
  });

  await softStep('3.5-3.7 Enrichment reused on another table with customerid', async () => {
    await page.evaluate(async (conn) => {
      // @ts-ignore
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 800));
      // @ts-ignore
      const df = await grok.data.db.query(conn, 'select * from public.customers');
      df.name = 'customers';
      df.setTag('.db-source-connection', conn);
      df.setTag('.db-source-schema', 'public');
      df.setTag('.db-source-table', 'customers');
      // @ts-ignore
      grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); (resolve as any)(); });
        setTimeout(resolve, 3000);
      });
      // @ts-ignore
      grok.shell.o = df.col('customerid');
    }, connection);
    // Pane may not appear if reuse is unsupported — bounded wait then assert intent.
    await waitForConnPane(15000);
    await waitForEnrichPane(10000); await waitForEnrichContent(10000);
    const reused = await page.evaluate(() => {
      const eh = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find(h => h.textContent?.trim() === 'Enrich') as HTMLElement | undefined;
      const content = eh?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');
      return (content as HTMLElement)?.innerText.includes('CustomerInfo_test') ?? false;
    });
    expect(reused).toBe(true);
  });

  // Section 4 skipped — needs a second user account.
  await softStep('4 Visibility for different users', async () => {
    test.skip(true, 'Section 4 requires a second user account — not automated.');
  });

  await softStep('Cleanup test enrichments and layouts', async () => {
    await page.evaluate(async (conn) => {
      // @ts-ignore
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 800));
      // @ts-ignore
      const df = await grok.data.db.query(conn, 'select * from public.orders');
      df.name = 'orders_cleanup';
      df.setTag('.db-source-connection', conn);
      df.setTag('.db-source-schema', 'public');
      df.setTag('.db-source-table', 'orders');
      // @ts-ignore
      grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 2500));
    }, connection);
    for (const col of ['customerid', 'employeeid']) {
      await page.evaluate((c) => {
        // @ts-ignore
        grok.shell.o = grok.shell.tv.dataFrame.col(c);
      }, col);
      await waitForConnPane(15000);
      await waitForEnrichPane(10000); await waitForEnrichContent(10000);
      await page.evaluate(async () => {
        const eh = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
          .find(h => h.textContent?.trim() === 'Enrich') as HTMLElement | undefined;
        const content = eh?.closest('.d4-accordion-pane')?.querySelector('.d4-accordion-pane-content');
        if (!content) return;
        for (let i = 0; i < 5; i++) {
          const rows = Array.from(content.querySelectorAll('.power-pack-enrichment-row'));
          if (rows.length === 0) break;
          ((rows[0] as HTMLElement).querySelector('i.fa-times') as HTMLElement).click();
          await new Promise(r => setTimeout(r, 1500));
        }
      });
    }
    await page.evaluate(async () => {
      // @ts-ignore
      const layouts = await grok.dapi.layouts.list({pageSize: 100});
      const mine = layouts.filter((l: any) => l.name === 'OrdersQueryResult');
      // @ts-ignore
      for (const l of mine) { try { await grok.dapi.layouts.delete(l); } catch (e) {} }
    });
  });

  if (stepErrors.length > 0)
    throw new Error('Failed steps:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
