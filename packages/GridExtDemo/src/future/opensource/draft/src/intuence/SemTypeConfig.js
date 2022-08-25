class SemTypeConfigEntry
{
    constructor(classSemType, strURLPathToIcon)
    {
      this.m_classSemType = classSemType;
      this.m_strURLPathToIcon = strURLPathToIcon === undefined ? null : strURLPathToIcon;
    }
}

export class SemTypeConfig
{
 constructor()
 {
  this.m_arEntries = new Array();
 }

 size()
 {
  return this.m_arEntries.length;
 }

 add(classSemType, strURLPathToIcon)
 {
     this.m_arEntries.push(new SemTypeConfigEntry(classSemType, strURLPathToIcon));
 }

 getSemTypeclass(nEntry)
 {
     return this.m_arEntries[nEntry].m_classSemType;
 }


 findSemType(arCols)
 {
     let fScore = 0;
     let fScoreMax = -1;
     let nEScoreMax = -1;
     let classSemType = null;
     let type = null;
     let typeBest = null;
     let entry = null;
     const nEntryCount = this.m_arEntries.length;
     for(var nE=0; nE<nEntryCount;  ++nE)
     {
        entry = this.m_arEntries[nE];
        classSemType = entry.m_classSemType;
        type = new classSemType(entry.m_strURLPathToIcon);

        fScore = type.getCompatibilityLevel(arCols);
        if(fScore > 0 && (fScoreMax === 0  || fScore > fScoreMax))
        {
            fScoreMax = fScore;
            nEScoreMax = nE;
            typeBest = type;
        }
     }

     return nEScoreMax < 0 ? null :  typeBest;
 }
}
