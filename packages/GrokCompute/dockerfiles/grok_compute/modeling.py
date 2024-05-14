from chemprop import ChemProp


engines = [ChemProp]


def get_engine_by_type(type):
    """
    Returns model engine by type.

    :param type: Modeling engine type.
    :return: Modeling engine instance.
    """
    engine = next((eng for eng in engines if eng.type == type), None)
    if engine is None:
        raise Exception('Unknown modeling engine name')
    return engine()


def get_all_engines():
    """
    Returns dictionary with all registered modeling engines with parameters and options.

    :return: Dictionary with all registered modeling engines with parameters and options.
    """
    engs = {}
    for eng in engines:
        engs[eng.type] = {
            'description': eng.description,
            'options': eng.options,
            'parameters': eng.parameters
        }
    return engs
