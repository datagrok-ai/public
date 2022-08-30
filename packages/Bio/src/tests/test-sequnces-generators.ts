export function generateManySequences(): string {
  let csvData = `MSA,Activity
meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me,5.30751`;
  for (let i = 0; i < 10 ** 6; i++) {
    csvData += `\n meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr//Phe_4Me,5.30751`;
  }
  return csvData;
}

export function generateLongSequence(): string {
  let longSequence = `meI/hHis/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/E/N/D-Orn/D-aThr`;
  for (let i = 0; i < 10 ** 6; i++) {
    longSequence += `/Aca/N/T/dE/Thr_PO3H2/Aca/D-Tyr_Et/Tyr_ab-dehydroMe/dV/dv`;
  }
  longSequence += `//Phe_4Me,5.30751`;
  let csvData = `MSA,Activity
  \n ${longSequence}`;
  for (let i = 0; i <= 10 ** 2; i++) {
    csvData += `\n ${longSequence}`;
  }
  return csvData;
}
