export class TanimotoResponseMessage
{
 constructor(arScores)
 {
     this.m_arScores = arScores;
 }

 getScores()
 {
     return this.m_arScores;
 }
}