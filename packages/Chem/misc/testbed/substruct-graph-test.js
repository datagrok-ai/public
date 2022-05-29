function testSubstructGraph(inputSet, inputQueries, n = 10) {
  time('substructGraph searches', n, () => {
    for (let query of inputQueries) {
      searchByGraph(inputSet, query);
    }
  });
}