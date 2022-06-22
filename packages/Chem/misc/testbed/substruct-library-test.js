function testSubstructLibrary(inputSet, inputQueries, n = 10) {
  time('substructLibrary init', 1, () => { buildLibrary(inputSet); });
  time('substructLibrary searches', n, () => {
    for (let query of inputQueries) {
      searchInLibrary(query);
    }
  });
  closeLibrary();
}