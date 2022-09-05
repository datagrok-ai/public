import {render} from "datagrok-api/ui";

export class ClassMap
{
    constructor()
    {
        this.m_map = new Map();
    }

    set(type, renderer)
    {
        this.m_map.set(type, renderer);
    }

    get(clazz)
    {

        if(clazz.prototype === undefined)
            throw new Error("The argument is not an instance of class.");

        let renderer = this.m_map.get(clazz);
        while((renderer === undefined || renderer === null) && clazz !== Object)
        {
            clazz = clazz.prototype.__proto__.constructor;
            renderer = this.m_map.get(clazz);
            let asd = 0;
        }

        return renderer;
    }
}