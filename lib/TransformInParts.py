from janim.imports import *
from MyTypst import *


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
        # assert all(isinstance(i, Item) for i in source)
        # assert all(isinstance(i, Item) for i in target)
        # assert len(source) == len(target)
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
            # params = {'move_durations', 'move_keywords', 'movtrs_keywords'}
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
        # print(kwargs)
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
        s_pats, t_pats = cls.assimilation(s_patterns, t_patterns)

        s_parts = s.get(flatten(s.multi_slices(*s_pats)), gapless)
        t_parts = t.get(flatten(t.multi_slices(*t_pats)), gapless)

        if from_moved:
            return cls.from_moved(s_parts, t_parts, **kwargs)
        else:
            return cls(s_parts, t_parts, **kwargs)

    @classmethod
    def assimilation(cls, pattern0, pattern1):
        # Handle cases where one of the patterns is None, empty list or Ellipsis
        if pattern0 in (None, [], ...):
            pattern0 = (pattern1 := tuple(pattern1))
        if pattern1 in (None, [], ...):
            pattern1 = (pattern0 := tuple(pattern0))

        merged_patterns = [
            cls.merge_items(*cls.merge_items(item0, item1))
            for item0, item1 in zip(pattern0, pattern1)
        ]

        # Unzip the merged patterns, but into two tuples
        return zip(*merged_patterns)

    @staticmethod
    def merge_items(item0: tuple, item1: tuple):
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
        # durations: Iterable[int] = None,
        # trs_keywords: Iterable[dict] = None,
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
            # durations=durations, trs_keywords=trs_keywords,
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
