async function time(name, n, foo) {
  let start = new Date();
  for (let k = 0; k < n; ++k) {
    foo();
  } 
  let stop = new Date();
  console.log(`${name}: ${(stop - start)/n} ms`);
}