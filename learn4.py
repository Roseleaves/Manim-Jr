from janim.imports import *
from RoseL import *


class MyHelloWorld(Timeline):
    def construct(self):

        hello = MyTypst('a^2+b^2=c^2')
        typst = MyTypst('b^2 = c^2 - a^2')
        manim = MyTypst('a^2 = c^2 - b^2')
        latex = MyTypst('a^2 = (c+b)(c-b)')

        sqsum = MyTypst('a^2 + 2a b + b^2 = c^2 + 2a b')
        sumsq = MyTypst('(a+b)^2 = c^2 + 2a b')
        sqsum.points.shift(DOWN)
        sumsq.points.shift(DOWN)

        self.show((hello))
        self.play(TransformInParts(
            (hello[0:2], hello[2], hello[3:5], hello[5], hello[6:8]),
            (typst[6:8], typst[5], typst[0:2], typst[2], typst[3:5]),
        ))
        self.play(TransformInParts.matching_patterns(
            (typst, (('a^2',), ('b^2',), ('c^2',), ('=',), ('-',))),
            (manim, ...),
        ))
        print(manim.multi_slices('""^2'), manim.multi_slices('-'))
        print(manim[('""^2', 0)])
        self.play(
            TransformInParts.matching_patterns(
                (manim, zip('=-abc', [[0]] * 3 + [[0, -1]] * 2)),
                (latex, ...),
                offset=1/6
            ),
            TransformInParts.matching_patterns(
                (manim, [('""^2', [1, 2])] * 2 + [('-', ), ('""^2', [0]), ]),
                (latex, zip([*'()+'] + ['""^2'], [[]] * 4)),
                from_moved=True,
                offset=1/6
            ),
        )
        self.play(Uncreate(latex))
        self.play(Create(hello))
        self.play(
            TransformInParts.from_transposed(
                hello, sqsum,
                ['a^2', '+b^2=c^2'],
                trs_keywords={"hide_src": False},
                from_moved=True
            )
        )
        self.forward()


if __name__ == '__main__':
    from janim.gui.anim_viewer import AnimViewer
    AnimViewer.views(MyHelloWorld().build())
