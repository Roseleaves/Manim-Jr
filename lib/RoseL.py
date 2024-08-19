from janim.imports import *


def get_between(string: str, first: str = None, last: str = None,
                include_first: bool = False, include_last: bool = False) -> str:
    '''
    Extracts a substring from 'string' that is between the substrings 'first' and 'last'.

    Parameters:
        string (str): The main string to search within.
        first (str, optional): The starting delimiter string. Defaults to None.
        last (str, optional): The ending delimiter string. Defaults to None.
        include_first (bool, optional): If True, include 'first' in the result. Defaults to False.
        include_last (bool, optional): If True, include 'last' in the result. Defaults to False.

    Returns:
        str: The extracted substring between 'first' and 'last', or an empty string if delimiters are not found.
    '''

    # Determine the start position
    start = string.find(first) if first else 0
    if start == -1:
        return ''  # Return empty string if 'first' is not found

    if not include_first:
        start += len(first)

    # Determine the end position
    end = string.find(last, start) if last else len(string)
    if end == -1:
        return ''  # Return empty string if 'last' is not found

    if include_last:
        end += len(last)

    return string[start:end]


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
        ordinals = list(flatten(ordinals))

        if isinstance(pattern, (tuple, list)):
            return [self.slices(i, ordinals) for i in pattern]
        else:
            return self.slices(pattern, ordinals)


class TransformInParts(AnimGroup):
    def __init__(
        self,
        source: Iterable[Item],
        target: Iterable[Item],
        from_moved: bool = False,
        durations: Iterable[int] | None = None,
        trs_keywords: Iterable[dict] = ({},),
        **kwargs
    ):
        t = trs_keywords
        t = t if type(t) is not dict else (t,)

        if durations:
            cd = it.cycle(durations) if type(durations) \
                is not int else it.cycle((durations,))
            for trs_kw in t:
                trs_kw['duration'] = next(cd)

        cycle_keywords = it.cycle(t)

        anims = []
        if not from_moved:
            for s, t in zip(source, target):
                anims.append(Transform(
                    s, t, **next(cycle_keywords)))
        else:
            anims = self.from_moved(
                source, target, trs_keywords, **kwargs)
        super().__init__(*anims, **kwargs)

    @staticmethod
    def from_moved(
        source, target,
        trs_keywords: Iterable[dict] = ({},),
        move_keywords: Iterable[dict] = ({},),
        joint_keywords: dict = {},
        **kwargs
    ):
        m = move_keywords
        m = m if type(m) is not dict else (m,)
        m_keywords = it.cycle(m)

        t = trs_keywords
        t = t if type(t) is not dict else (t,)
        t_keywords = it.cycle(t)

        anims = []
        for s, t in zip(source, target):
            assert isinstance(s, Group)
            assert isinstance(t, Group)
            # 只对宽度较窄的进行移动
            if s.points.box.width < t.points.box.width:
                st = s.copy().points.move_to(t).r
                sk, tk = m_keywords, t_keywords
            else:
                if all(same_shape(ss, tt) for ss, tt in zip(s, t)):
                    anims.append(Transform(s, t, **next(t_keywords)))
                    continue
                st = t.copy().points.move_to(s).r
                sk, tk = t_keywords, m_keywords
            anims.append(Succession(
                Transform(s, st, **next(sk), show_target=False),
                Transform(st, t, **next(tk)),
                **joint_keywords
            ))
        return AnimGroup(*anims, **kwargs)

    @classmethod
    def from_segments(
        cls,
        source: Item,
        source_segments: Iterable[Iterable[int]] | Iterable[int],
        target: Item | types.EllipsisType,
        target_segments: Iterable[Iterable[int]] | Iterable[int] | types.EllipsisType,
        **kwargs
    ):
        if target_segments is ...:
            target_segments = source_segments

        tuples = []
        for segs in (source_segments, target_segments):
            standard = []
            assert len(segs) > 0
            if not isinstance(segs[0], Iterable):
                segs = [segs]

            for seg in segs:
                for l, r in it.pairwise(seg):
                    standard.append((l, r) if l >= 0 > r or l < r else (r, l))

            tuples.append(standard)

        s = [source[l:r] for l, r in tuples[0]]
        t = [target[l:r] for l, r in tuples[1]]
        return cls(s, t, **kwargs)

    @classmethod
    def matching_patterns(cls, source, target, gapless: bool = False,
                          from_moved: bool = False, **kwargs):
        s, s_patterns = source
        t, t_patterns = target
        assert isinstance(s, MyTypst)
        assert isinstance(t, MyTypst)

        # 每段动画的整理不影响下一段动画
        s_pats, t_pats = cls._assimilation(s_patterns, t_patterns)

        s_parts = s.get(flatten(s.multi_slices(*s_pats)), gapless)
        t_parts = t.get(flatten(t.multi_slices(*t_pats)), gapless)

        if from_moved:
            return cls.from_moved(s_parts, t_parts, **kwargs)
        else:
            return cls(s_parts, t_parts, **kwargs)

    @classmethod
    def _assimilation(cls, pattern0, pattern1):
        # Handle cases where one of the patterns is None, empty list or Ellipsis
        if pattern0 in (None, [], ...):
            pattern0 = (pattern1 := tuple(pattern1))
        if pattern1 in (None, [], ...):
            pattern1 = (pattern0 := tuple(pattern0))

        merged_patterns = [
            cls._merge_items(*cls._merge_items(item0, item1))
            for item0, item1 in zip(pattern0, pattern1)
        ]

        # Unzip the merged patterns, but into two tuples
        return zip(*merged_patterns)

    @staticmethod
    def _merge_items(item0: tuple, item1: tuple):
        '''merge item0 from item1, and then swap them.
        use it twice to merge them each other.'''

        k1, v1, *_ = item1 + ([],) if item1 else (None, [])
        k0, v0, *rest = item0 + (v1,) if item0 else (k1, v1)
        # tuple 逗号优先级低于 if 三目运算, 用 tuple 加法直接提供缺省值

        # Replace None or empty strings with counterpart's values
        # MyTypst.slice() 代入 None, [], ... 都是全选,
        # 这里填 None, [] 和不填也是全选, 但填 ... 就是选对家
        k0 = k0 if k0 is not ... else k1
        v0 = v0 if v0 is not ... else v1

        # Ensure v0 is a list
        v0 = v0 if isinstance(v0, list) else list(item0[1:]) if rest else [v0]

        return (k1, v1), (k0, v0)

    @classmethod
    def from_transposed(
        cls,
        source: MyTypst,
        target: MyTypst,
        source_patterns: Iterable[str],
        target_patterns: Iterable[str] | types.EllipsisType = ...,
        source_indices: Iterable[int] | types.EllipsisType = [0],
        target_indices: Iterable[int] | types.EllipsisType = [...],
        **anim_kwargs
    ):
        sp = source_patterns
        sp = (sp,) if isinstance(sp, str) else tuple(sp)
        tp = target_patterns
        tp = sp if tp is ... else (tp,) if type(tp) is str else tuple(tp)
        if all(isinstance(i, int) for i in tp):
            # 忘了打target_patterns的省略号. target_indices就不手动顺延了
            source_indices = tp
            tp = sp
        return cls.matching_patterns(
            (source, zip(sp, it.cycle(source_indices))),
            (target, zip(tp, it.cycle(target_indices))),
            **anim_kwargs)


class TransformMatchingPatterns(Succession):
    def __init__(
        self,
        *Typst_fragments,
        trs_keywords: dict = {},
        step_keywords: list[dict] = [],
        **succ_keywords
    ):
        anims = []
        cycled_keywords = it.cycle(step_keywords)
        for source, target in it.pairwise(Typst_fragments):
            anims.append(TransformInParts(
                source, target,
                trs_keywords=trs_keywords,
                **next(cycled_keywords))
            )
        super().__init__(*anims, **succ_keywords)
