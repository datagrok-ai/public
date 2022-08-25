export class Insets
{
  constructor(nT, nL, nB, nR)
  {
      this.m_nT = nT;
      this.m_nL = nL;
      this.m_nB = nB;
      this.m_nR = nR;
  }


  getT()
  {
      return this.m_nT;
  }

  getL()
  {
      return this.m_nL;
  }

  getB()
  {
      return this.m_nB;
  }


  getR()
  {
      return this.m_nR;
  }

  clone()
  {
      return new Insets(this.m_nT, this.m_nL, this.m_nB, this.m_nR);
  }
}