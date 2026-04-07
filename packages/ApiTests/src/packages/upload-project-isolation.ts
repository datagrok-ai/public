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
    const url = `${grok.dapi.root}/packages/dev/${key}/${name}?debug=${debug}&rebuild=false`;
    const r = await fetch(url, {method: 'POST', body: packageData as BodyInit});
    if (r.status !== 200) {
      const body = await r.text();
      throw new Error(`Publish ${name} failed with status ${r.status}: ${body.substring(0, 200)}`);
    }
    const text = await r.text();
    if (text.indexOf('ApiError') !== -1)
      throw new Error(`Publish ${name} returned ApiError: ${text.substring(0, 200)}`);
    await DG.delay(5000);
  }

  async function safeDeletePackage(name: string) {
    try {
      await fetch(`${grok.dapi.root}/packages/dev/${key}/${name}`, {method: 'DELETE'});
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
    await DG.delay(2000);

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
    await DG.delay(2000);

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

  // When a new package version is published, dashboards that reference entities
  // from the old version should automatically get relations to the new version's entities.
  test('Dashboard: entity propagation across package versions', async () => {
    await safeDeletePackage('pkgisolationa');
    await cleanupProjects();
    await cleanupDashboards();
    await DG.delay(2000);

    // Publish release version
    await publish('pkgisolationa', false);

    // Get a function entity from the package
    const funcs = await grok.dapi.functions
      .filter('package.shortName = "Pkgisolationa"').list();
    expect(funcs.length >= 1, true,
      `Expected at least 1 function from package, got ${funcs.length}`);
    const funcV1 = funcs[0];

    // Create a dashboard and add function as a link
    const dashboard = DG.Project.create();
    dashboard.name = 'DashboardRelationTest';
    dashboard.addLink(funcV1);
    const saved = await grok.dapi.projects.save(dashboard);

    // Publish debug version (triggers _addProjectRelationsForPackageVersion)
    await publish('pkgisolationa', true);

    // Get all functions across all versions (including debug)
    const allFuncs = await grok.dapi.functions
      .filter('package.shortName = "Pkgisolationa"')
      .allPackageVersions().list();

    // Reload dashboard and check its links
    const reloaded = await grok.dapi.projects.find(saved.id);
    const linkedIds = new Set(reloaded.links.map((e: DG.Entity) => e.id));

    // There should be a function from v2 with same bind_id in the dashboard
    const v2Funcs = allFuncs.filter((f: DG.Func) => f.id !== funcV1.id);
    expect(v2Funcs.length >= 1, true,
      `Expected functions from debug version, got ${v2Funcs.length}`);
    const hasV2Entity = v2Funcs.some((f: DG.Func) => linkedIds.has(f.id));
    expect(hasV2Entity, true,
      'Dashboard should contain entity from the new package version');

    // Cleanup
    await grok.dapi.projects.delete(reloaded);
  }, {timeout: 120000, owner: 'aparamonov@datagrok.ai'});
}, {
  owner: 'aparamonov@datagrok.ai'
});
