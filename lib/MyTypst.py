from janim.imports import *

def flatten(iterable):
    if not isinstance(iterable, (list, tuple)):
        return iterable,
    return list(it.chain.from_iterable(
        flatten(x)
        for x in iterable))


def same_shape(vitem1: VItem, vitem2: VItem, threshold: float = 1/20) -> bool:
    '''
    Determines whether two VItem objects have the same shape within a specified tolerance.

    Parameters:
        vitem1 (VItem): The first VItem object to compare.
        vitem2 (VItem): The second VItem object to compare.
        threshold (float, optional): The tolerance level for considering two shapes as the same.
                                     Defaults to 1/20.

    Returns:
        bool: True if the two shapes are considered the same within the given threshold, False otherwise.
    '''

    # Extract subpaths from both items
    paths1, paths2 = (j.points.get_subpaths() for j in (vitem1, vitem2))

    # Check if the number of subpaths is the same
    if len(paths1) != len(paths2):
        return False
    for path1, path2 in zip(paths1, paths2):
        if path1.shape != path2.shape:
            return False

    # Normalize both point sets by centering and scaling to height 1
    (center1, height1), (center2, height2) = (
        (item.points.box.center, item.points.box.height)
        for item in (vitem1, vitem2)
    )

    min_height = 1e-3
    if height1 < min_height or height2 < min_height:
        height1, height2 = (item.points.box.width
                            for item in (vitem1, vitem2))

    # Compare normalized paths
    for path1, path2 in zip(paths1, paths2):
        normalized_path1 = (path1 - center1) / height1
        normalized_path2 = (path2 - center2) / height2
        if np.sum(np.isclose(normalized_path1, normalized_path2)) / np.size(normalized_path1) < 1 - threshold:
            return False

    return True


class MyTypst(Typst):

    def __init__(self, arg: str | TypstDoc, **kwargs):
        match arg:
            case str():
                # 如果arg是字符串，按照常规方式初始化
                super().__init__(arg, **kwargs)
            case TypstDoc():
                # 如果arg是MyTypst实例，复制其属性
                self.__dict__.update(arg.__dict__)
            case _:
                raise TypeError(
                    'Argument must be a string or an instance of TypstDoc.')

    def indices(formula, pattern: 'MyTypst', cmp=same_shape) -> list[int]:
        '''Find all occurrences of a pattern within this object.

        Args:
            pattern (MyTypst): The pattern to search for.

        Returns:
            list[int]: A list of starting indices where the pattern is found.
        '''
        if not isinstance(pattern, MyTypst):
            pattern = MyTypst(pattern)

        lps = pattern._compute_lps(cmp)
        indices, p = [], 0

        for index, shape in enumerate(formula):
            while not (same := cmp(shape, pattern[p])) and p != 0:
                p = lps[p - 1]
            if same:
                p += 1
            if p == len(pattern):
                indices.append(index - (p - 1))
                p = lps[p - 1]

        return indices

    lps = {}

    def _compute_lps(self, cmp):
        # 原则上不提供对分数线等很难有 MyTypst.text 成员的比较
        if self.text not in MyTypst.lps:
            lps = [0] * len(self)
            for index, shape in enumerate(self):
                p, same = index, False
                while p > 0 and not same:
                    p = lps[p - 1]
                    same = cmp(shape, self[p])
                if same:
                    p += 1
                lps[index] = p
            MyTypst.lps[self.text] = lps
        return MyTypst.lps[self.text]

    def __getitem__(self, key: int | slice | Iterable):
        # 可能不会处理脱靶情形，因为空的 Typst 很难定义
        match key:
            case int() | slice():
                return super().__getitem__(key)
            case list() if all(isinstance(x, int) for x in key):
                return Group(*map(super().__getitem__, key))
            case list() if all(isinstance(x, bool) for x in key):
                return Group(*[i for i, j in zip(self, key) if j])

            case str(pattern):
                return self.get(self.slices(pattern, 0))
            case str(), int():
                return self.get(self.slices(*key))
            case str(), list():
                return Group(*flatten(self.get(self.slices(*key))))
            case _ if self._with_empty(key):
                return self
            case dict():
                return Group(*flatten(self.get(self.multi_slices(key))))
            case tuple() | list():
                return Group(*flatten(self.get(self.multi_slices(*key))))

    @staticmethod
    def _with_empty(arg):
        empty_signs = None, [], (), "", ...
        if arg in empty_signs:
            return True
        for i in flatten(arg):
            if i in empty_signs:
                return True
        return False

    def get(self, slices, gapless: bool = False):
        fragments = []
        if not gapless:
            if isinstance(slices, slice):
                return self[slices]
            else:
                for i in slices:
                    fragment = self.get(i)
                    fragments.append(fragment)
        else:
            indices = {0, len(self)}
            for i in flatten(slices):
                assert isinstance(i, slice)
                indices.update({i.start, i.stop})
            for start, stop in it.pairwise(sorted(indices)):
                fragments.append(self[start:stop])
        return fragments

    def slices(self, pattern: "MyTypst", ordinal: list[int] = None):
        match pattern:
            case _ if self._with_empty(pattern):
                return slice(0, len(self))
            case str():
                pat = MyTypst(pattern)
            case TypstDoc():
                pat = pattern

        indices = self.indices(pat)

        slices = []
        match ordinal:
            case int(i) if i < len(indices):
                slices = slice(indices[i], indices[i] + len(pat))
            case _ if self._with_empty(ordinal):
                if len(indices) == 1:
                    slices = slice(indices[0], indices[0] + len(pat))
                else:
                    slices = [slice(i, i + len(pat)) for i in indices]
            case list():
                slices = [slice(indices[i], indices[i] + len(pat))
                          for i in ordinal]
        return slices

    def multi_slices(self, *args):
        # args = list(args)
        while len(args) == 1 and isinstance(args, (tuple, dict)):
            if isinstance(args, dict):
                args = tuple(args.items())
            args, = args

        if isinstance(args, str) or \
                all(isinstance(i, str) for i in args):
            args = args, None

        if all(isinstance(i, (tuple, dict)) for i in args):
            return [self.multi_slices(i) for i in args]

        pattern, *ordinals = args
        # assert isinstance(pattern, (TypstDoc, str, list))
        # assert all(isinstance(i, (int, list)) for i in ordinals)
        ordinals = list(flatten(ordinals))

        if isinstance(pattern, (tuple, list)):
            return [self.slices(i, ordinals) for i in pattern]
        else:
            return self.slices(pattern, ordinals)

