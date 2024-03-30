import multiprocessing
import numpy as np

def set_type(value, value_type):
    """
    Returns dictionary with input value and value type.

    :param value:
    :param value_type:
    :return: Dictionary with input value and value type.
    """
    return {'type': value_type, 'value': value.tolist()}

def np_none(shape):
    """
    Returns np.array initialized with None elements.

    :param shape: see np.array
    :return: np.array initialized with NaNs.
    """
    return np.full(shape, None, dtype=object)


def _parallel_worker(id, results, target, args):
    results[id] = target(*args)


def parallelize(target, args, splittable_args, num_cores=4, split_threshold=64):
    """
    Gets descriptors for each input molecule.

    :param target: Parallelizable function
    :param args: List of all function arguments
    :param splittable_args: List of function arguments that will be splitted between cores
    :param num_cores: Number of cores for task
    :param split_threshold: threshold to start splitting task to multiple cores
    :return: "target" function output
    """
    len_arg = len(splittable_args[0])
    if len_arg >= split_threshold:
        if len_arg < num_cores:
            num_cores = len_arg
        segment_length = len_arg // num_cores
        manager = multiprocessing.Manager()
        results = manager.dict()
        processes = []
        for n in range(0, num_cores):
            _args = ()
            for arg in args:
                if arg in splittable_args:
                    _args += (arg[(n * segment_length):((n + 1) * segment_length if n != num_cores - 1 else len_arg)],)
                else:
                    _args += (arg,)
            process = multiprocessing.Process(target=_parallel_worker, args=(n, results, target, _args))
            processes.append(process)
            process.start()
        for process in processes:
            process.join()
        if isinstance(results[0], list):
            result = results[0]
            for n in range(1, len(results)):
                result.extend(results[n])
            return result
        elif isinstance(results[0], dict):
            result = {}
            keys = set()
            for value in results.values():
                keys.update(value.keys())
            for n in range(0, len(results)):
                for key in keys:
                    if key not in result:
                        result[key] = {'value': [], 'type': None}
                    if key in results[n]:
                        result[key]['value'].extend(results[n][key]['value'])
                        result[key]['type'] = results[n][key]['type']
                    else:
                        result[key]['value'].extend([None] * (segment_length if n != num_cores - 1 else
                                                              len_arg - n * segment_length))
            return result
    else:
        return target(*args)
