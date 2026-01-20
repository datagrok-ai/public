// Demonstrates working with Spaces - hierarchical containers for organizing entities and files.
(async () => {
  const spaceName = 'DemoSpace_' + Math.random().toString(36).substring(7);

  // Check if a root space with this name already exists
  const exists = await grok.dapi.spaces.rootSpaceExists(spaceName);

  if (exists) {
    grok.shell.info(`Space "${spaceName}" already exists`);
    return;
  }

  // Create a new root space
  const space = await grok.dapi.spaces.createRootSpace(spaceName);
  console.log(`Created space: ${space.friendlyName} (id: ${space.id})`);

  // Get a client for the space by its ID
  const spaceClient = grok.dapi.spaces.id(space.id);

  // === Working with files ===

  // Upload a file to the space
  const fileName = 'hello.txt';
  await spaceClient.files.writeString(fileName, 'Hello from Spaces API!');
  console.log(`Uploaded file: ${fileName}`);

  // Check if file exists
  const fileExists = await spaceClient.files.exists(fileName);
  console.log(`File exists: ${fileExists}`);

  // Read file content
  const content = await spaceClient.files.readAsString(fileName);
  console.log(`File content: ${content}`);

  // === Working with subspaces ===

  // Create a subspace by creating a directory (directories become subspaces)
  const subspaceName = 'SubFolder';
  await spaceClient.files.createDirectory(subspaceName);
  console.log(`Created subspace via directory: ${subspaceName}`);

  // Check if subspace exists
  const subspaceExists = await spaceClient.subspaceExists(subspaceName);
  console.log(`Subspace exists: ${subspaceExists}`);

  // Alternatively, create a subspace directly
  const childSpace = await spaceClient.addSubspace('ChildSpace');
  console.log(`Created child space: ${childSpace.friendlyName}`);

  // Move file to subspace
  await spaceClient.files.move([fileName], subspaceName);
  console.log(`Moved file to ${subspaceName}`);

  // Verify file was moved
  const movedFileExists = await spaceClient.files.exists(`${subspaceName}/${fileName}`);
  console.log(`File in subspace: ${movedFileExists}`);

  // === Working with children ===

  // List all children of the space
  let children = await spaceClient.children.list();
  console.log(`Space has ${children.length} children`);
  children.forEach(child => console.log(`  - ${child.friendlyName} (${child.constructor.name})`));

  // Filter children by type (only subspaces/projects)
  const subspaces = await spaceClient.children.filter('Project', false).list();
  console.log(`Space has ${subspaces.length} subspaces`);

  // === Working with entities ===

  // Create a script and add it to the space
  // Note: In a real scenario, you would use an existing entity
  let script = DG.Script.create(`#name: Template
#description: Calculates number of cells in the table
#language: python
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]

count = table.shape[0] * table.shape[1]`);
  script = await grok.dapi.scripts.save(script);
  await spaceClient.addEntity(script.id);
  console.log(`Added script to space`);

  // === Cleanup ===

  // Delete the space (cascades to subspaces and files)
  await grok.dapi.spaces.delete(space);
  console.log(`Deleted space: ${spaceName}`);

  grok.shell.info('Spaces API demo completed! Check console for details.');
})();
