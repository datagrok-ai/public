export class InitRequestMessage
{
    constructor(eCanvas)
    {
        this.m_eCanvas = eCanvas;

    }

    getCanvas()
    {
        return this.m_eCanvas;
    }

}