from copy import copy

example_dict = { 'key1' : 'value1',
                 'key2' : 'value2',
                 'key3' : { 'key3a': 'value3a' },
                 'key4' : { 'key4a': { 'key4aa': 'value4aa',
                                       'key4ab': 'value4ab',
                                       'key4ac': 'value1'},
                            'key4b': 'value4b'}
                }


def get_keys(d, target):
    result = []
    path = []
    for k, v in d.iteritems():
        path.append(k)
        if isinstance(v, dict):
            get_keys(v, target)
        if v == target:
            result.append(copy(path))
        path.pop()

    return result


print get_keys(example_dict, 'value1')        