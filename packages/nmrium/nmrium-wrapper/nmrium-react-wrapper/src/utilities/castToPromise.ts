export function castToPromise(func, ...props) {
  return new Promise((resolve) => {
    const spectra = [];
    func(spectra, ...props);
    resolve(spectra);
  });
}
