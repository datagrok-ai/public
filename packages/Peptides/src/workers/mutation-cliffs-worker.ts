
onmessage = async (event): Promise<void> => {
  const {startIdx, endIdx, activityArray, monomerInfoArray, settings, currentTargetIdx, targetOptions} = event.data;
  // const monomers1: string[] = [];
  // const monomers2: string[] = [];
  const pos: string[] = [];
  const seq1Idxs: number[] = [];
  const seq2Idxs: number[] = [];
  const chunkSize = endIdx - startIdx;
  //const mi = startRow;
  //const mj = startCol;
  let cnt = 0;
  const startRow = activityArray.length - 2 - Math.floor(
    Math.sqrt(-8 * startIdx + 4 * activityArray.length * (activityArray.length - 1) - 7) / 2 - 0.5);
  const startCol = startIdx - activityArray.length * startRow + Math.floor((startRow + 1) * (startRow + 2) / 2);
  let seq1Idx = startRow;
  let seq2Idx = startCol;
  const tempData = new Array(monomerInfoArray.length);
  while (cnt < chunkSize) {
    if (!(currentTargetIdx !== -1 && (targetOptions.targetCol?.rawData[seq1Idx] !== currentTargetIdx ||
        targetOptions.targetCol?.rawData[seq2Idx] !== currentTargetIdx))) {
      let substCounter = 0;
      const activityValSeq1 = activityArray[seq1Idx];
      const activityValSeq2 = activityArray[seq2Idx];
      const delta = activityValSeq1 - activityValSeq2;
      if (Math.abs(delta) >= settings.minActivityDelta) {
        let substCounterFlag = false;
        let tempDataIdx = 0;
        for (const monomerInfo of monomerInfoArray) {
          const seq1categoryIdx = monomerInfo.rawData[seq1Idx];
          const seq2categoryIdx = monomerInfo.rawData[seq2Idx];
          if (seq1categoryIdx === seq2categoryIdx)
            continue;

          substCounter++;
          substCounterFlag = substCounter > settings.maxMutations;
          if (substCounterFlag)
            break;

          tempData[tempDataIdx++] = {
            pos: monomerInfo.name,
            // seq1monomer: monomerInfo.cat![seq1categoryIdx],
            // seq2monomer: monomerInfo.cat![seq2categoryIdx],
            seq1Idx: seq1Idx,
            seq2Idx: seq2Idx,
          };
        }
        if (!(substCounterFlag || substCounter === 0)) {
          for (let i = 0; i < tempDataIdx; i++) {
            const tempDataElement = tempData[i];
            const position = tempDataElement.pos;
            // const seq1monomer = tempDataElement.seq1monomer;
            // const seq2monomer = tempDataElement.seq2monomer;
            // monomers1.push(seq1monomer);
            // monomers2.push(seq2monomer);
            pos.push(position);
            seq1Idxs.push(seq1Idx);
            seq2Idxs.push(seq2Idx);
          }
        }
      }
    }
    cnt++;
    seq2Idx++;
    if (seq2Idx === activityArray.length) {
      seq1Idx++;
      seq2Idx = seq1Idx + 1;
    }
  }
  postMessage({
    // monomers1: monomers1,
    // monomers2: monomers2,
    pos: pos,
    seq1Idxs: new Uint32Array(seq1Idxs),
    seq2Idxs: new Uint32Array(seq2Idxs),
  });
};
