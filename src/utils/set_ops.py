def shared_words(*args):

    """
    :param args:    an arbitrary number of sets (other structures are coerced to set)
    :return:        the intersection between the input arguments
    """

    shared = set()
    for i, arg in enumerate(args):
        if arg:
            if not isinstance(arg, set):
                print("Input argument {}, originally of type {}, is coerced to set".format(i, type(arg)))
                arg = set(arg)
            if not shared:
                shared = arg
            else:
                shared = shared.intersection(arg)

    return shared


def union_words(*args):

    """
    :param args:    an arbitrary number of sets (other structures are coerced to set)
    :return:        the union of the input arguments
    """

    union_set = set()
    for i, arg in enumerate(args):
        if arg:
            if not isinstance(arg, set):
                print("Input argument {}, originally of type {}, is coerced to set".format(i, type(arg)))
                arg = set(arg)
            if not union_set:
                union_set = arg
            else:
                union_set = union_set.union(arg)

    return union_set