/* ---
sub_features_covered: [sharing.api.check, sharing.api.get, sharing.api.grant, sharing.api.revoke, sharing.server.privileges-router]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';

test.use(specTestOptions);

test('Sharing — API permissions (grant / get / check / revoke lifecycle via grok.dapi.permissions)', async ({page}) => {
  // Pure JS-API grant/get/check/revoke round-trips, single actor, no UI lifecycle.
  test.setTimeout(120_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  
  
  
  
  
  await softStep('Setup: create throwaway Script entity and resolve stable test group', async () => {
    const res = await page.evaluate(async () => {
      const g = (window as any).grok, DG = (window as any).DG;
      try {
        const script = await g.dapi.scripts.save(DG.Script.create('sharing_api_test_' + Date.now()));
        
        
        
        const candidates = await g.dapi.groups.filter('All users').list({limit: 50});
        const group = candidates.find((x: any) =>
          x && (x.name === 'AllUsers' || x.friendlyName === 'All users'));
        (window as any).__sharingApiTest = {scriptId: script.id, groupId: group ? group.id : null};
        return {
          ok: true,
          scriptId: script.id,
          groupId: group ? group.id : null,
          groupName: group ? (group.friendlyName || group.name) : null,
        };
      } catch (e) { return {ok: false, reason: String(e).slice(0, 300)}; }
    });
    expect(res.ok, `setup must succeed (${(res as any).reason ?? ''})`).toBe(true);
    expect(res.scriptId, 'script entity must be created').toBeTruthy();
    expect(res.groupId, 'all-users test group must resolve').toBeTruthy();
    console.log(`[sharing-api] Setup: script=${res.scriptId}, group='${res.groupName}' (${res.groupId})`);
  });

  
  await softStep('Scenario 1: grant view → get shows group in view (not edit) → revoke → get shows group absent', async () => {
    const res = await page.evaluate(async () => {
      const g = (window as any).grok;
      try {
        const {scriptId, groupId} = (window as any).__sharingApiTest;
        const script = await g.dapi.scripts.find(scriptId);
        const group = await g.dapi.groups.find(groupId);
        
        await g.dapi.permissions.grant(script, group, false);
        
        let perms = await g.dapi.permissions.get(script);
        
        
        const afterGrant = {
          viewHasGroup: (perms.view || []).some((x: any) => x && x.id === groupId),
          editHasGroup: (perms.edit || []).some((x: any) => x && x.id === groupId),
        };
        
        await g.dapi.permissions.revoke(group, script);
        
        perms = await g.dapi.permissions.get(script);
        const afterRevoke = {
          viewHasGroup: (perms.view || []).some((x: any) => x && x.id === groupId),
          editHasGroup: (perms.edit || []).some((x: any) => x && x.id === groupId),
        };
        return {ok: true, afterGrant, afterRevoke};
      } catch (e) { return {ok: false, reason: String(e).slice(0, 300)}; }
    });
    expect(res.ok, `lifecycle calls must run (${(res as any).reason ?? ''})`).toBe(true);
    
    expect((res as any).afterGrant.viewHasGroup, 'group must appear in view array after view grant').toBe(true);
    expect((res as any).afterGrant.editHasGroup, 'group must NOT appear in edit array after view grant').toBe(false);
    
    expect((res as any).afterRevoke.viewHasGroup, 'group must be absent from view array after revoke').toBe(false);
    expect((res as any).afterRevoke.editHasGroup, 'group must be absent from edit array after revoke').toBe(false);
    console.log('[sharing-api] Scenario 1: view grant placed group in view; revoke cleared it');
  });

  
  await softStep('Scenario 2: grant edit → get shows group in edit array → revoke cleans up', async () => {
    const res = await page.evaluate(async () => {
      const g = (window as any).grok;
      try {
        const {scriptId, groupId} = (window as any).__sharingApiTest;
        const script = await g.dapi.scripts.find(scriptId);
        const group = await g.dapi.groups.find(groupId);
        
        await g.dapi.permissions.grant(script, group, true);
        
        let perms = await g.dapi.permissions.get(script);
        const afterGrant = {
          editHasGroup: (perms.edit || []).some((x: any) => x && x.id === groupId),
        };
        
        await g.dapi.permissions.revoke(group, script);
        perms = await g.dapi.permissions.get(script);
        const afterRevoke = {
          editHasGroup: (perms.edit || []).some((x: any) => x && x.id === groupId),
          viewHasGroup: (perms.view || []).some((x: any) => x && x.id === groupId),
        };
        return {ok: true, afterGrant, afterRevoke};
      } catch (e) { return {ok: false, reason: String(e).slice(0, 300)}; }
    });
    expect(res.ok, `edit-grant calls must run (${(res as any).reason ?? ''})`).toBe(true);
    expect((res as any).afterGrant.editHasGroup, 'group must appear in edit array after edit grant').toBe(true);
    expect((res as any).afterRevoke.editHasGroup, 'group must be absent from edit array after revoke').toBe(false);
    expect((res as any).afterRevoke.viewHasGroup, 'group must be absent from view array after revoke').toBe(false);
    console.log('[sharing-api] Scenario 2: edit grant placed group in edit; revoke cleared it');
  });

  
  
  
  
  
  await softStep('Scenario 3: check() returns owner-perspective boolean; view/edit boundary verified via get()', async () => {
    const res = await page.evaluate(async () => {
      const g = (window as any).grok;
      try {
        const {scriptId, groupId} = (window as any).__sharingApiTest;
        const script = await g.dapi.scripts.find(scriptId);
        const group = await g.dapi.groups.find(groupId);

        
        const checkEdit = await g.dapi.permissions.check(script, 'Edit');
        const checkView = await g.dapi.permissions.check(script, 'View');
        const checkShare = await g.dapi.permissions.check(script, 'Share');

        
        await g.dapi.permissions.grant(script, group, false);
        let perms = await g.dapi.permissions.get(script);
        const viewGrant = {
          inView: (perms.view || []).some((x: any) => x && x.id === groupId),
          inEdit: (perms.edit || []).some((x: any) => x && x.id === groupId),
        };
        await g.dapi.permissions.revoke(group, script);

        
        await g.dapi.permissions.grant(script, group, true);
        perms = await g.dapi.permissions.get(script);
        const editGrant = {
          inView: (perms.view || []).some((x: any) => x && x.id === groupId),
          inEdit: (perms.edit || []).some((x: any) => x && x.id === groupId),
        };
        await g.dapi.permissions.revoke(group, script);

        return {
          ok: true,
          checkTypes: {
            edit: typeof checkEdit, view: typeof checkView, share: typeof checkShare,
          },
          checkValues: {edit: checkEdit, view: checkView, share: checkShare},
          viewGrant, editGrant,
        };
      } catch (e) { return {ok: false, reason: String(e).slice(0, 300)}; }
    });
    expect(res.ok, `check() / boundary calls must run (${(res as any).reason ?? ''})`).toBe(true);
    
    expect((res as any).checkTypes.edit, 'check(Edit) must return a boolean').toBe('boolean');
    expect((res as any).checkTypes.view, 'check(View) must return a boolean').toBe('boolean');
    
    expect((res as any).checkValues.edit, 'owner must have Edit on own script').toBe(true);
    expect((res as any).checkValues.view, 'owner must have View on own script').toBe(true);
    expect((res as any).checkValues.share, 'owner must have Share on own script').toBe(true);
    
    expect((res as any).viewGrant.inView, 'view grant → group in view array').toBe(true);
    expect((res as any).viewGrant.inEdit, 'view grant → group NOT in edit array').toBe(false);
    expect((res as any).editGrant.inEdit, 'edit grant → group in edit array').toBe(true);
    expect((res as any).editGrant.inView, 'edit grant → group NOT in view array').toBe(false);
    console.log('[sharing-api] Scenario 3: check() is owner-perspective boolean; view/edit boundary confirmed via get()');
  });

  
  
  
  
  await softStep('Scenario 4: grant two groups, get reflects both, revoke both, get reflects none', async () => {
    const res = await page.evaluate(async () => {
      const g = (window as any).grok;
      try {
        const {scriptId, groupId} = (window as any).__sharingApiTest;
        const script = await g.dapi.scripts.find(scriptId);
        const groupA = await g.dapi.groups.find(groupId); 
        
        const all = await g.dapi.groups.list({limit: 50});
        const groupB = all.find((x: any) => x && x.id && x.id !== groupId);
        if (!groupB) return {ok: false, reason: 'could not resolve a second distinct group'};

        
        await g.dapi.permissions.grant(script, groupA, false);
        await g.dapi.permissions.grant(script, groupB, false);
        
        let perms = await g.dapi.permissions.get(script);
        let viewIds = (perms.view || []).map((x: any) => x && x.id);
        const bothPresent = viewIds.includes(groupA.id) && viewIds.includes(groupB.id);
        
        await g.dapi.permissions.revoke(groupA, script);
        await g.dapi.permissions.revoke(groupB, script);
        
        perms = await g.dapi.permissions.get(script);
        viewIds = (perms.view || []).map((x: any) => x && x.id);
        const neitherPresent = !viewIds.includes(groupA.id) && !viewIds.includes(groupB.id);
        return {ok: true, bothPresent, neitherPresent, groupBId: groupB.id};
      } catch (e) { return {ok: false, reason: String(e).slice(0, 300)}; }
    });
    expect(res.ok, `batch router calls must run (${(res as any).reason ?? ''})`).toBe(true);
    expect((res as any).bothPresent, 'both groups must appear in view array after two grants').toBe(true);
    expect((res as any).neitherPresent, 'both groups must be gone after two independent revokes').toBe(true);
    console.log(`[sharing-api] Scenario 4: two-group batch grant/revoke tracked independently (groupB=${(res as any).groupBId})`);
  });

  
  await page.evaluate(async () => {
    const g = (window as any).grok;
    try {
      const holder = (window as any).__sharingApiTest;
      if (holder && holder.scriptId) {
        const script = await g.dapi.scripts.find(holder.scriptId).catch(() => null);
        if (script) {
          
          if (holder.groupId) {
            const group = await g.dapi.groups.find(holder.groupId).catch(() => null);
            if (group) await g.dapi.permissions.revoke(group, script).catch(() => {});
          }
          await g.dapi.scripts.delete(script).catch(() => {});
        }
      }
    } catch (_) {  }
  }).catch(() => {});

  if (stepErrors.length > 0)
    throw new Error('Spec step failures:\n' + stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
});
