import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {before, category, expect, test, after} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';

category('Packages: project isolation', () => {

  let packageData: Uint8Array;
  let key: string;

  before(async () => {
    packageData = await _package.files.readAsBytes('package/packageWithProject.zip');
    const response = await fetch(`${grok.dapi.root}/users/current/dev_key?getNew=false`,
      {method: 'GET', credentials: 'include'});
    key = await response.json();
  });

  after(async () => {
    await safeDeletePackage('pkgisolationa');
    await safeDeletePackage('pkgisolationb');
    await cleanupProjects();
    await cleanupDashboards();
  });

  async function publish(name: string, debug: boolean = false) {
    const url = `${grok.dapi.root}/packages/dev/${name}?debug=${debug}&rebuild=false`;
    const r = await fetch(url, {method: 'POST', body: packageData as BodyInit, headers: {Authorization: `Dev ${key}`}});
    if (r.status !== 200) {
      const body = await r.text();
      throw new Error(`Publish ${name} failed with status ${r.status}: ${body.substring(0, 200)}`);
    }
    const text = await r.text();
    if (text.indexOf('ApiError') !== -1)
      throw new Error(`Publish ${name} returned ApiError: ${text.substring(0, 200)}`);
  }

  async function safeDeletePackage(name: string) {
    try {
      await fetch(`${grok.dapi.root}/packages/dev/${name}`, {method: 'DELETE', headers: {Authorization: `Dev ${key}`}});
    }
    catch (_) {}
  }

  async function cleanupProjects() {
    try {
      const projects = await grok.dapi.projects.filter('IsolationTestProject').list();
      for (const p of projects)
        await grok.dapi.projects.delete(p);
    }
    catch (_) {}
  }

  async function cleanupDashboards() {
    try {
      const projects = await grok.dapi.projects.filter('DashboardRelationTest').list();
      for (const p of projects)
        await grok.dapi.projects.delete(p);
    }
    catch (_) {}
  }

  async function getTestProjects(): Promise<DG.Project[]> {
    return await grok.dapi.projects.filter('IsolationTestProject').list();
  }

  // When two different packages contain projects with the same hardcoded IDs,
  // each package should get its own project entity. The second publish should
  // NOT steal the entity from the first package.
  test('Hardcoded IDs: cross-package isolation', async () => {
    await safeDeletePackage('pkgisolationa');
    await safeDeletePackage('pkgisolationb');
    await cleanupProjects();

    // Publish package A with a project containing hardcoded IDs
    await publish('pkgisolationa');

    const projectsAfterA = await getTestProjects();
    expect(projectsAfterA.length >= 1, true,
      `Package A should have created a project entity, got ${projectsAfterA.length}`);
    const projectAId = projectsAfterA[0].id;

    // Publish package B with the SAME content (same hardcoded project/table IDs)
    await publish('pkgisolationb');

    // Each package should have its own project entity
    const projectsAfterB = await getTestProjects();
    expect(projectsAfterB.length >= 2, true,
      `Publishing a different package with same hardcoded project IDs should create a new entity, not steal the existing one. Got ${projectsAfterB.length} projects`);

    // Original project from package A should still be accessible
    const origProject = await grok.dapi.projects.find(projectAId);
    expect(origProject != null, true,
      'Original project from package A must still exist after publishing package B');
  }, {timeout: 120000, owner: 'aparamonov@datagrok.ai'});

  // When a new version of the same package is published, it should create its own
  // project entities. The hardcoded IDs from the project JSON should be remapped,
  // not reused to steal entities from the previous version.
  test('Hardcoded IDs: version isolation', async () => {
    await safeDeletePackage('pkgisolationa');
    await cleanupProjects();

    // Publish release version
    await publish('pkgisolationa', false);

    const projectsV1 = await getTestProjects();
    expect(projectsV1.length >= 1, true,
      `Release version should create a project entity, got ${projectsV1.length}`);
    const projectV1Id = projectsV1[0].id;

    // Publish debug version of the same package (creates a second version)
    await publish('pkgisolationa', true);

    const projectsAfterV2 = await getTestProjects();
    expect(projectsAfterV2.length >= 1, true,
      `After debug publish, should have a project entity, got ${projectsAfterV2.length}`);

    // The debug version should have created its own entity, not reused the release version's
    const v2Ids = new Set(projectsAfterV2.map((p) => p.id));
    const v1EntityStolen = v2Ids.size === 1 && v2Ids.has(projectV1Id) && projectsAfterV2.length === 1;
    expect(v1EntityStolen, false,
      'New version should create its own project entity with remapped IDs, not steal from the previous version');
  }, {timeout: 120000, owner: 'aparamonov@datagrok.ai'});

  // A debug version is one developer's: publishing it leaves the dashboards that link the release
  // version's entities alone (core#1244). The links move to the current version when a version is
  // deleted or a release is published, which datlas package_version_relations_test covers.
  test('Dashboard: a debug publish does not link the debug version', async () => {
    await safeDeletePackage('pkgisolationa');
    await cleanupProjects();
    await cleanupDashboards();
    let saved: DG.Project | null = null;
    try {
      await publish('pkgisolationa', false);
      const funcs = await grok.dapi.functions.filter('package.shortName = "Pkgisolationa"').list();
      expect(funcs.length >= 1, true, `Expected at least 1 function from package, got ${funcs.length}`);
      const released = funcs[0];

      const dashboard = DG.Project.create();
      dashboard.name = 'DashboardRelationTest';
      dashboard.addLink(released);
      saved = await grok.dapi.projects.save(dashboard);
      const linksOf = async () => (await grok.dapi.projects.find(saved!.id)).links.map((e: DG.Entity) => e.id);
      const linkedBefore = await linksOf();
      expect(linkedBefore.includes(released.id), true,
        `The saved dashboard does not link the function it was given: ${linkedBefore.join()}`);

      await publish('pkgisolationa', true);
      const debugIds = (await grok.dapi.functions.filter('package.shortName = "Pkgisolationa"')
        .allPackageVersions().list()).map((f: DG.Func) => f.id).filter((id: string) => id !== released.id);
      expect(debugIds.length >= 1, true, `Expected functions from the debug version, got ${debugIds.length}`);

      // Read as the developer who published the debug version, the release function is hidden behind
      // it, so the dashboard may show no link here; what must not appear is a debug-version entity.
      const linked = await linksOf();
      expect(linked.some((id: string) => debugIds.includes(id)), false,
        `The dashboard links a function of the debug version: ${linked.join()}`);
      return `links after the debug publish: ${linked.join() || 'none'}`;
    } finally {
      if (saved != null)
        await grok.dapi.projects.delete(saved).catch(() => {});
    }
  }, {timeout: 120000, owner: 'aparamonov@datagrok.ai'});
}, {
  owner: 'aparamonov@datagrok.ai'
});
