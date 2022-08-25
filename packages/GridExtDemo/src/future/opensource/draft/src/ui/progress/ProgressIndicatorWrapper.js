export class ProgressIndicatorWrapper
{
    constructor(progress) {
        this.m_progress = progress;

        this.m_strUpper = "";
    }

    setUpperText(str)
    {
      this.m_strUpper = str;
    }

    setBottomText(str)
    {
        this.m_progress.update(this.m_progress.percent, this.m_strUpper + " " + str);
    }

    setProgress(nPct)
    {
        this.m_progress.update(nPct, this.m_progress.description);
    }

    close() {}
}

