/** Grok User */
class User {
    constructor(d) { this.d = d; }
    static fromId(id) { return new User(grok_User_From_Id(id)); }
    static current() { return new User(grok_User()); }

    get firstName() { return grok_User_Get_FirstName(this.d); }
    get lastName() { return grok_User_Get_LastName(this.d); }
    get email() { return grok_User_Get_Email(this.d); }
    get picture() { return grok_User_Get_Picture(this.d); }
    get login() { return grok_User_Get_Login(this.d); }

    toMarkup() { return grok_User_ToMarkup(this.d); }
}


class Entity {
    constructor(d) { this.d = d;}

    get id() { return grok_Entity_Get_Id([this.d]); }
    set id(x) { return grok_Entity_Set_Id([this.d], x); }
    get name() { return grok_Entity_Get_Name([this.d]); }
    set name(x) { return grok_Entity_Set_Name([this.d], x); }
    get path() { return grok_Entity_Path([this.d]); }
}


/** Represents a project */
class Project extends Entity {
    constructor(d) { super(d); }

    get description() { return grok_Project_Description(this.d); }
    get isDirty() { return grok_Project_IsDirty(this.d); }
    get isEmpty() { return grok_Project_IsEmpty(this.d); }

    toMarkup() { return grok_Project_ToMarkup(this.d); }
}


/** Represents a data query */
class DataQuery extends Entity {
    constructor(d) { super(d); }

    get query() { return grok_Query_Query([this.d]); }
}


/** Represents a data connection */
class DataConnection extends Entity {
    constructor(d) { super(d); }

    get parameters() { return grok_DataConnection_Parameters([this.d]); }
}


/** Represents a credentials */
class Credentials extends Entity {
    constructor(d) { super(d); }

    get parameters() { return grok_Credentials_Parameters([this.d]); }
}


class GrokPackage {
    constructor(webRoot) {
        this.webRoot = webRoot;
    }

    // Override init() method to provide package-specific initialization.
    // It is guaranteed to get called exactly once before the execution of any function below.
    /*async*/ init() { return Promise.resolve(null); }
}


class JsEntityMeta {
    get type() { throw 'Not defined.'; }

    isApplicable(x) { throw 'Not defined.'; }

    /** String representation of the [item], by default item.toString(). */
    getCaption(x) { return `${x}`; };

    renderIcon(x) { return ui.divText(this.getCaption(x)); }

    renderMarkup(x) { return ui.divText(this.getCaption(x)); }
    renderTooltip(x) { return ui.divText(this.getCaption(x)); }
    renderCard(x) { return ui.divText(this.getCaption(x)); }
    renderProperties(x) { return ui.divText(this.getCaption(x)); }
    renderView(x) {return this.renderProperties(x); }

    /** Gets called once upon the registration of meta class. */
    init() {}

    static register(meta) { grok_Meta_Register(meta); }

    registerParamFunc(name, run) {
        grok.functions.registerParamFunc(name, this.type, run, this.isApplicable);
    }
}


/** Strongly-typed property associated with an object.
 *  Used for reflection, serialization, UI generation, and other introspection-dependent tasks. */
class Property {
    constructor(d) { this.d = d; }

    get get() { return grok_Property_Get_Get(this.d); }
    set get(x) { return grok_Property_Set_Get(this.d, x); }

    get set() { return grok_Property_Get_Set(this.d); }
    set set(x) { return grok_Property_Set_Set(this.d, x); }

    get name() { return grok_Property_Get_Name(this.d); }
    set name(s) { grok_Property_Set_Name(this.d, s); }

    /** Property type */
    get propertyType() { return grok_Property_Get_PropertyType(this.d); }
    set propertyType(s) { grok_Property_Set_PropertyType(this.d, s); }

    /** Semantic type */
    get semType() { return grok_Property_Get_SemType(this.d); }
    set semType(s) { grok_Property_Set_SemType(this.d, s); }

    /** Description */
    get description() { return grok_Property_Get_Description(this.d); }
    set description(s) { grok_Property_Set_Description(this.d, s); }

    /** Default value */
    get defaultValue() { return grok_Property_Get_DefaultValue(this.d); }
    set defaultValue(s) { grok_Property_Set_DefaultValue(this.d, s); }

    static create(name, type, getter, setter, defaultValue = null) { return new Property(grok_Property(name, type, getter, setter, defaultValue)); }
    static int(name, getter, setter, defaultValue) { return Property.create(name, TYPE_INT, getter, setter, defaultValue); }
    static float(name, getter, setter, defaultValue) { return Property.create(name, TYPE_FLOAT, getter, setter, defaultValue); }
    static string(name, getter, setter, defaultValue) { return Property.create(name, TYPE_STRING, getter, setter, defaultValue); }
}

class SemanticValue {
    constructor(d) { this.d = d; }

    static fromValueType(value, semType) { return new SemanticValue(grok_SemanticValue(value, semType)); }

    get value() { return grok_SemanticValue_Get_Value(this.d); }
    set value(x) { grok_SemanticValue_Set_Value(this.d, x); }

    get semType() { return grok_SemanticValue_Get_SemType(this.d); }
    set semType(x) { grok_SemanticValue_Set_SemType(this.d, x); }
}

class DateTime {
    constructor(d) { this.d = d; }

    static fromDate(date) { return DateTime.fromMillisecondsSinceEpoch(date.getTime()); }
    static fromMillisecondsSinceEpoch(millisecondsSinceEpoch) {
        return new DateTime(grok_DateTime_FromMillisecondsSinceEpoch(millisecondsSinceEpoch));
    }
}
