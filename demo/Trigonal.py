from functools import reduce
from janim.imports import *
from RoseL import *

with open('Trigonal.typ', 'r', encoding='utf-8') as file:
    content = file.read()


class MyTrigonal(Timeline):

    def complete(self, *items: Iterable[Item], anim=Write):
        self.forward(1/60)
        self.play(anim(Group(*filter(
            lambda x: not self.is_displaying(x), it.chain(*items)
        ))))
        for i in items:
            if isinstance(i, TypstDoc):
                self.show(i)
        if isinstance(items, TypstDoc):
            self.show(items)
        self.forward(1/60)
        return

    def construct(self):
        TrsP = TransformInParts.from_transposed

        def TransformSegmentWithTrace(s1, s2):
            anims = [Transform(segment[s1], segment[s2], hide_src=False)]
            for i, j in zip(s1, s2):
                if i + j in segment:
                    anims.append(Create(segment[i + j]))

            return self.play(*anims)
        TrST = TransformSegmentWithTrace

        def get_equal(formula, ith=0):
            return formula['=', ith].points.box.center

        self.formula = TypstDoc(content)
        form2 = self.formula.copy()
        form2.color.set("#CC99CC")
        form2.fill.set(alpha=0.5)

        self.show(form2)
        play_parts = [1, 2, 3, 4]

############ FORMULA I, III(ii) ############
        if 1 in play_parts:
            coord_A = np.array([3, 0, 0])
            len_AB = 3

            angle = np.arctan(3/4)
            arc = Arc(0, angle, len_AB, arc_center=coord_A)

            def unit(theta):
                return np.array([np.cos(theta), np.sin(theta), 0])

            vec = {
                "AB": len_AB * unit(0),
                "AC": len_AB * unit(angle),
                "AD": len_AB * np.cos(angle) * unit(0),
                "BE": len_AB * np.tan(angle) * unit(PI/2),
            }
            vec["AE"] = vec["AB"] + vec["BE"]

            coord = {P: coord_A + vec[f"A{P}"] for P in "BCDE"}
            coord["A"] = coord_A
            segment = {P1+P2: Line(coord[P1], coord[P2], alpha=1/2) for (P1, P2) in
                       ["AB", "BC", "CA", "CD", "BE", "EA",]}
            seg_title = {"CA": Typst("1"),
                         "BC": Typst("x"),
                         "CD": Typst("sin x"),
                         "BE": Typst("tan x")}
            for i in seg_title:
                seg_title[i].points.scale(2/3).\
                    move_to(segment[i]).shift(1/3 * LEFT)  # .rotate(-PI/2)
            seg_title["BC"].points.scale(3/2).shift(1/2 * UR)
            seg_title["BE"].points.shift(2/3 * RIGHT)

            tan_trig = Polygon(*map(coord.get, "ABE"),
                               color=TEAL, fill_alpha=0.2)
            rad_sect = Sector(angle=angle, radius=len_AB,
                              arc_center=coord_A, color=YELLOW, fill_alpha=0.2)
            sin_trig = Polygon(*map(coord.get, "ABC"),
                               color=GOLD, fill_alpha=0.2)

            self.play(Create(arc))
            self.play(Write(seg_title["BC"]))
            self.complete(segment.values(), anim=lambda x: AnimGroup(
                *[Create(i) for i in x], offset=1/2
            ))
            self.complete(seg_title.values())
            self.forward()

            areas = MyTypst("1/2 r^2 sin x < 1/2 r^2 x < 1/2 r^2 tan x")
            less_index = areas.indices("<")
            self.play(FadeIn(tan_trig),
                      Write(areas[less_index[1] + 1:]))
            self.forward()
            self.play(FadeIn(rad_sect),
                      Write(areas[less_index[0] + 1:less_index[1] + 1]))
            self.forward()
            self.play(FadeIn(sin_trig),
                      Write(areas[:less_index[0] + 1]))
            self.forward()

            self.play(*[FadeOut(areas[i:i+5]) for i in [
                0, less_index[0] + 1, less_index[1] + 1
            ]], offset=1/5)

            order = MyTypst("sin x < x < tan x quad (exists a>0, 0<x<a)")
            self.play(TrsP(
                areas, order, ["sin", "x", "<", "x", "<", "tan", "x"],
                [0, 0, 0, 1, 1, 0, 2],  # 它们在areas里面是分离的
                offset=1/5))
            self.complete(order)

            self.play(Uncreate(Group(
                *segment.values(), *seg_title.values(),
                arc, sin_trig, rad_sect, tan_trig,
            )), duration=2)
            self.forward()

            sin_special = MyTypst("sin 0 = 0")
            sin_special.points.align_to(order, LEFT).shift(DOWN)

            self.play(
                Indicate(order["0<x"], duration=2),
                TrsP(order, sin_special,
                     "sin x < x", '',
                     trs_keywords={"hide_src": False}),
                offset=1/2)
            self.forward()

            self.play(
                Transform(order, self.formula[:23]),
                Write(self.formula[23:26]),
                Transform(sin_special, self.formula[87:93]),
                Write(self.formula[93:96]),
                offset=1/6
            )

            self.forward()

############ FORMULA II(ii), III(i), IV(i), V ############
        if 2 in play_parts:
            # 作图场景，用欧拉直角三角形和两次投影导入
            alpha = 1
            beta = (np.sqrt(5)-1)/2
            coord_A = np.array([3, 0, 0])
            len_AB = 3

            def unit(theta):
                return np.array([np.cos(theta), np.sin(theta), 0])

            vec = {
                "AB": len_AB * unit(alpha-beta),
                "AC": len_AB * np.cos(beta) * unit(alpha),
                "AD": len_AB * np.cos(beta) * np.cos(alpha) * unit(0),
                "BE": len_AB * np.sin(beta) * np.sin(alpha) * unit(PI),
                "AF": len_AB * np.cos(alpha-beta) * unit(0),
            }
            vec["AE"] = vec["AB"] + vec["BE"]

            coord = {P: coord_A + vec[f"A{P}"] for P in "BCDEF"}
            coord["A"] = coord_A

            segment = {P1+P2: Line(coord[P1], coord[P2], color=color)
                       for (P1, P2), color in zip(
                ["AB", "AC", "BC", "AD", "CD", "EC", "BE", "FD", "BF"],
                [YELLOW, TEAL, GOLD, BLUE, GREY, GREY, RED, RED, GREY]
            )}
            seg_title = {"AB": Typst("1"),
                         "AC": Typst("cos beta"),
                         "BC": Typst("sin beta"),
                         "AD": Typst("cos beta cos alpha"),
                         "BE": Typst("sin beta sin alpha"),
                         "FD": Typst("cos (alpha - beta)")}
            for i in seg_title:
                seg_title[i].points.scale(2/3).\
                    move_to(segment[i]).shift(1/5 * DOWN)
            seg_title["AB"].points.shift(1/3 * UP).rotate(alpha-beta)
            seg_title["AC"].points.shift(1/5 * UL).rotate(alpha)
            seg_title["BC"].points.shift(2/5 * UP).rotate(alpha - PI/2)
            seg_title["BE"].points.shift(2/5 * UP)
            seg_title["FD"].points.shift(1/4 * LEFT)

            angle = {"alpha": Arc(0, alpha, 1/3, arc_center=coord_A),
                     "beta": Arc(alpha-beta, beta, 1/2, arc_center=coord_A),
                     "alpha-beta": Arc(0, alpha-beta, 5/6, arc_center=coord_A)}
            for i in angle:
                title = Typst(i)
                title.points.scale(2/3).next_to(angle[i], RIGHT)
                angle[i] = Group(angle[i], title)
            angle["beta"][1].points.shift(1/6 * UL)

            self.play(Create(segment["AB"]))
            TrST("AB", "AC")
            self.play(Create(angle["beta"]), *[
                Write(seg_title[i]) for i in ["AB", "AC", "BC"]])
            TrST("AC", "AD")
            self.play(Create(angle["alpha"]),
                      Write(seg_title["AD"]),
                      )
            TrST("BC", "EC")
            self.play(Write(seg_title["BE"]))
            TrST("BE", "FD")
            self.forward()

            # 余弦差角公式, cosine of the difference of angles formula
            cd_str = "cos (alpha - beta) = cos alpha cos beta + sin alpha sin beta"
            cos_diff = MyTypst(cd_str)
            cos_diff.points.shift(1.2 * UP)

            self.play(Create(angle["alpha-beta"]))
            self.play(Write(seg_title["FD"]))
            self.play(Write(cos_diff))
            self.play(Uncreate(Group(*segment.values(), *
                                     seg_title.values(), *angle.values())), duration=2)
            self.forward()

################ 以上是作图环节. 作图是用来引入. 后面是代数. ################

            cd0_str = [cd_str.replace(i, "0") for i in ["alpha", "beta"]]
            cos_diff_0 = MyTypst(cd0_str[1])
            cos_diff_0.points.shift(0.6 * UP)
            cos_0_diff = MyTypst(cd0_str[0])
            self.play(Transform(cos_diff, cos_diff_0, hide_src=False),
                      Transform(cos_diff, cos_0_diff, hide_src=False,
                                duration=2, rate_func=double_smooth))

            self.play(
                Indicate(self.formula[87:93], duration=2),
                FadeOut(Group(
                    cos_diff_0[["(", "- 0 )", get_between(
                        cd0_str[1], "+", "", True)]],
                    cos_0_diff[["0", get_between(
                        cd0_str[0], "+", "", True)], 0]
                )),
                offset=1
            )

            cos_0 = MyTypst("cos 0 = 1")
            cos_0.points.shift(0.6 * UP).align_to(
                cos_diff_0["cos 0"], LEFT)
            self.play(TrsP(
                cos_diff_0, cos_0,
                ["cos 0", "="] + ["cos", "alpha"]*2,
                ["cos 0", "="] + ["1"] * 4,
                [0] * 4 + [1] * 2, [0],
                offset=1/6
            ))
            self.play(FadeOut(cos_0_diff["cos 0"]))

            cos_negative = MyTypst("cos (-beta) = cos beta")
            cos_negative.points.shift(0.85*RIGHT)
            self.play(TrsP(
                cos_0_diff, cos_negative,
                ["beta", "cos", ] * 2 + [*"=)-("],
                [1, -1] + [0] * 6,
                offset=1/6
            ))

            self.forward()

            cos_same = MyTypst(cd_str.replace("beta", "alpha"))
            cos_same.points.shift(0.6 * DOWN)
            self.play(Transform(cos_diff, cos_same, hide_src=False))

            # 正弦-余弦的勾-股关系. 勾 (ko) 股 (ka) 是上古拟音.
            koka = MyTypst("cos^2 alpha + sin^2 alpha=1")
            koka.points.shift(0.6 * DOWN).align_to(
                cos_same["cos", 1], LEFT)

            self.play(TrsP(
                cos_same, koka,
                ["cos", "alpha", "sin", "alpha"] * 2 + ["+"], ...,
                [1, 2, 0, 4, 2, 3, 1, 5, 0],
                ([0] * 3 + [1]) * 2 + [0],
                offset=1/6
            ))
            self.complete(koka[:-2])

            self.play(TrsP(
                cos_same, cos_0,
                ["cos", "(alpha - alpha)", "="],
                ["cos", "0", "="],
            ))
            self.play(Transform(cos_0[-2:], koka[-2:], hide_src=False))
            # 尽量只讲主干. 除了一个tangent是sine/cosine之外其他同角三角函数就不多说了
            # 零角和负角只是一个抽象示意, 是函数的拓展. 但平方和就是不简单的关系, 可以补上一张完全不同的图

            st = np.sin(alpha)
            ct = np.cos(alpha)
            coord_O = np.array([4, 0, 0])
            gien = 2
            vec = {
                "OC": st * ct * UP,
                "OA": ct * ct * LEFT,
                "AB": RIGHT,
                "AD": DR,
                "AF": ct * UR,
                "BE": st * UR
            }
            vec["OB"] = vec["OA"] + vec["AB"]
            vec["OD"] = vec["OA"] + vec["AD"]
            vec["OE"] = vec["OB"] + vec["BE"]
            vec["OF"] = vec["OA"] + vec["AF"]
            coord = {i: coord_O + gien * vec["O"+i] for i in "ABCDEF"}

            seg_title = {"AB": Typst("1"),
                         "AC": Typst("cos alpha"),
                         "BC": Typst("sin alpha")}
            for (i, j), title in seg_title.items():
                title.points.scale(2/3).move_to(
                    (coord[i]+coord[j])/2 + 1/8 * UP
                )

            right_triangle = Polygon(*map(coord.get, "ABC"))
            square = {i+j: Rect(coord[i], coord[j], fill_alpha=1/3, color=k)
                      for (i, j), k in zip(["AD", "BE", "AF"], [YELLOW, TEAL, GOLD])}
            square["AD"].points.rotate(PI/2)
            square["BE"].points.rotate(PI)\
                .rotate(alpha, about_point=coord["B"])
            square["AF"].points.rotate(3 * PI/2)\
                .rotate(alpha, about_point=coord["A"])
            # 矩形的默认顶点顺序是 (UR, UL, DL, DR)
            # 自旋用来把首边贴在三角形上
            self.play(
                Create(right_triangle),
                *map(Create, seg_title.values()),
                *map(Create, square.values()),
                offset=0.5)
            self.forward()

            # 沿着BC边对AC边的正方形做切变
            sqAF = square["AF"].copy()
            sqBE = square["BE"].copy()

            self.play(
                sqAF.anim.points.apply_matrix(
                    np.array([[1/ct, -st], [0, ct]]) @
                    np.array([[ct, st], [-st, ct]]),
                    about_point=coord["A"]),
                sqBE.anim.points.apply_matrix(
                    np.array([[ct, -1/st], [st, 0]]) @
                    np.array([[ct, st], [-st, ct]]),
                    about_point=coord["B"])
            )
            self.play(Rotate(sqAF, -PI/2, about_point=coord["A"]),
                      Rotate(sqBE, PI/2, about_point=coord["B"]))
            self.play(
                sqAF.anim.points.apply_matrix(
                    np.array([[1, 0], [-st/ct, 1]]),
                    about_point=coord["A"]),
                sqBE.anim.points.apply_matrix(
                    np.array([[1, 0], [ct/st, 1]]),
                    about_point=coord["B"])
            )

            self.forward()

            self.play(Uncreate(Group(
                right_triangle, sqAF, sqBE, *seg_title.values(), *square.values(),
            )))

            self.play(
                Transform(cos_diff, self.formula[52:78]),
                Transform(cos_0, self.formula[81:87]),
                Transform(cos_negative, self.formula[96:108]),
                Transform(koka, self.formula[124:137]),
                Write(self.formula[78:81]),
                Write(self.formula[93:96]),
                Write(self.formula[121:124]),
                Write(self.formula[137:140]),
                offset=1/6
            )
            self.forward()

############ FORMULA II(i, ii) ############
        if 3 in play_parts:

            cd_str = "cos (alpha - beta) = cos alpha cos beta + sin alpha sin beta"
            cos_diff = MyTypst(cd_str)
            cos_diff.points.shift(1.2 * UP)

            # 前面的变换都相当容易, 至少没有需要拉长或者缩短公式的紧急需求
            # 但这里开始就变得复杂, 不得不祭出来 "先移动再变换", 有的时候还得要 copy().
            # inversed difference, 差角的差角, 反转的差距
            cos_inv_diff = MyTypst(cd_str.replace("beta", "(alpha - beta)"))
            cos_inv_diff.points.shift(np.array([
                (get_equal(cos_diff) - get_equal(cos_inv_diff))[0], 0.6, 0]))

            self.play(FadeIn(cos_diff))
            self.play(TrsP(
                cos_diff.copy(), cos_inv_diff,
                "beta", "(alpha - beta)", [None],
                gapless=True,  # offset=1/12,
                from_moved=True,
            ))

            de_str = "cos beta = cos alpha (cos alpha cos beta + sin alpha sin beta) + sin alpha sin (alpha - beta)"
            double_expanded = MyTypst(de_str)
            double_expanded.points.shift(np.array([
                (get_equal(cos_diff) - get_equal(double_expanded))[0], 0.6, 0,]))

            bra = cos_inv_diff["( alpha -", [0, 1]]
            ket = cos_inv_diff[")", [0, 1]]

            self.play(
                AnimGroup(Transform(bra[1], bra[0]),
                          Transform(ket[0], ket[1])),
                AnimGroup(FadeOut(bra[0]), FadeOut(ket[1])),
                offset=1
            )
            self.play(
                TrsP(cos_inv_diff, cos_diff, "cos (alpha - beta)"),
                TrsP(cos_inv_diff, double_expanded,
                     ["+ sin alpha sin (alpha - beta)",
                      "= cos alpha", "beta", "cos", ],
                     ),
                TrsP(cos_diff, double_expanded,
                     "cos alpha cos beta + sin alpha sin beta",
                     trs_keywords={"hide_src": False}
                     ),
                offset=0.5
            )
            self.complete(double_expanded)

            sd1_str = "sin alpha sin (alpha - beta) = \
            cos beta - cos alpha cos alpha cos beta -\
            cos alpha sin alpha sin beta"
            sine_diff_1 = MyTypst(sd1_str)
            sine_diff_1.points.shift(np.array([
                (get_equal(cos_diff) - get_equal(sine_diff_1))[0], 0, 0,]))
            sd1_pats = ["cos beta", "cos alpha", "cos alpha cos beta",
                        "cos alpha", "sin alpha sin beta",
                        "sin alpha sin (alpha - beta)", "="]

            self.play(TrsP(double_expanded, sine_diff_1,
                           sd1_pats, sd1_pats,
                           [0], [0, 0, 0, 2, 0, 0, 0],
                           trs_keywords={"hide_src": False},
                           offset=1/3
                           ))
            self.complete(sine_diff_1)

            sine_diff_2 = MyTypst("cos beta(1 - cos alpha cos alpha)")
            sine_diff_2.points.move_to(
                sine_diff_1["cos beta"], aligned_edge=LEFT)

            self.play(TrsP(
                sine_diff_1, sine_diff_2,
                ["- cos alpha cos alpha", "cos beta"], [0, -1],
            ))
            self.complete(sine_diff_2)
            sd2_cc = sine_diff_2["cos alpha cos alpha"]

            sine_diff_3 = MyTypst("sin alpha sin alpha")
            sine_diff_3.points.move_to(sd2_cc, aligned_edge=DOWN)
            self.play(
                Uncreate(sine_diff_2[")"]),
                Uncreate(sine_diff_2["(1-"]),
                Transform(sd2_cc, sine_diff_3),
                offset=1/3
            )

            sd_str = "sin (alpha - beta) = sin alpha cos beta - cos alpha sin beta"
            sine_diff = MyTypst(sd_str)
            sine_diff.points.shift(0.6 * DOWN)
            self.play(TrsP(
                sine_diff_1, sine_diff,
                ["sin(alpha - beta) =", "cos beta", "-cos alpha", "sin beta"],
                [0, 0, -1, 0],
                trs_keywords={"hide_src": False}, offset=1/3
            ))
            self.complete(sine_diff)
            self.forward()
            self.play(FadeOut(Group(
                cos_diff, double_expanded,
                sine_diff_1, sine_diff_2[:4], sine_diff_3)))

################ 以上是代数. 后面是作图, 用来形象展示. ################
            alpha = 1
            beta = (np.sqrt(5)-1)/2
            coord_A = np.array([3, 0, 0])
            len_AB = 3

            def unit(theta):
                return np.array([np.cos(theta), np.sin(theta), 0])

            vec = {
                "AB": len_AB * unit(alpha-beta),
                "AC": len_AB * np.cos(beta) * unit(alpha),
                "AD": len_AB * np.cos(beta) * np.cos(alpha) * unit(0),
                "BE": len_AB * np.sin(beta) * np.sin(alpha) * unit(PI),
                "AF": len_AB * np.cos(alpha-beta) * unit(0),
            }
            vec["AE"] = vec["AB"] + vec["BE"]

            coord = {P: coord_A + vec[f"A{P}"] for P in "BCDEF"}
            coord["A"] = coord_A

            segment = {P1+P2: Line(coord[P1], coord[P2], color=color)
                       for (P1, P2), color in zip(
                ["AB", "AC", "BC", "AD", "CD", "EC", "BE", "DF", "ED", "BF"],
                [YELLOW, TEAL, GOLD, GREY, BLUE, RED, GREY, GREY, BLUE, BLUE]
            )}
            seg_title = {"AB": Typst("1"),
                         "AC": Typst("cos beta"),
                         "BC": Typst("sin beta"),
                         "CD": Typst("cos beta sin alpha"),
                         "EC": Typst("sin beta cos alpha"),
                         "BF": Typst("sin (alpha - beta)")}
            for i in seg_title:
                seg_title[i].points.scale(2/3).\
                    move_to(segment[i]).shift(1/5 * DOWN)
            seg_title["AB"].points.shift(1/10 * UR).rotate(alpha-beta)
            seg_title["AC"].points.shift(1/5 * UL).rotate(alpha)
            seg_title["BC"].points.shift(2/5 * UP).rotate(alpha - PI/2)

            angle = {"alpha": Arc(0, alpha, 1/3, arc_center=coord_A),
                     "beta": Arc(alpha-beta, beta, 1/2, arc_center=coord_A),
                     "alpha-beta": Arc(0, alpha-beta, 5/6, arc_center=coord_A)}
            for i in angle:
                title = Typst(i)
                title.points.scale(2/3).next_to(angle[i], RIGHT)
                angle[i] = Group(angle[i], title)
            angle["beta"][1].points.shift(1/6 * UL)

            self.play(Create(segment["AB"]))
            TrST("AB", "AC")
            self.play(Create(angle["beta"]), *[
                Write(seg_title[i]) for i in ["AB", "AC", "BC"]])
            TrST("AC", "AD")
            self.play(Create(angle["alpha"]),
                      Write(seg_title["CD"]),
                      )
            TrST("BC", "EC")
            self.play(Write(seg_title["EC"]))
            TrST("ED", "BF")
            self.play(Create(angle["alpha-beta"]))
            self.play(Write(seg_title["BF"]))
            self.forward()

            self.play(Uncreate(Group(
                *segment.values(), *seg_title.values(), *angle.values()
            )), duration=2)
            self.forward()

################ 以上是作图环节. 作图是用来引入. 后面是代数. ################

            sine_0_diff = MyTypst(sd_str.replace("alpha", "0"))
            sine_0_diff.points.shift(np.array([
                (get_equal(sine_diff) - get_equal(sine_0_diff))[0], -1.2, 0,]))
            sine_negative = MyTypst("sin (-beta) = - sin beta")
            sine_negative.points.shift(np.array([
                (get_equal(sine_diff) - get_equal(sine_negative))[0], -1.2, 0,]))

            self.play(TrsP(
                sine_diff.copy(), sine_0_diff, "alpha", "0", [None],
                gapless=True, from_moved=True, offset=1/6
            ))
            self.play(
                Indicate(self.formula[81:93], duration=3),
                FadeOut(sine_0_diff[(["0", "sin 0 cos beta", "cos 0"], 0)]),
                TrsP(
                    sine_0_diff, sine_negative,
                    ["sin (", "- beta)=", "-", "sin beta"],  [-1],
                    offset=1/6
                ),
                offset=0.5
            )

            sine_sum_1 = MyTypst(sd_str.replace("beta", "( - beta)"))
            sine_sum_1.points.shift(np.array([
                (get_equal(sine_diff) - get_equal(sine_sum_1))[0], -1.8, 0,]))
            sine_sum = MyTypst(sd_str.replace("-", "+"))
            sine_sum.points.shift(np.array([
                (get_equal(sine_diff) - get_equal(sine_sum))[0], -1.8, 0,]))

            self.play(TrsP(
                sine_diff.copy(), sine_sum_1, "beta", "( - beta)", [None],
                gapless=True, from_moved=True, offset=1/6
            ))
            self.play(
                Indicate(self.formula[96:108]),
                FadeOut(sine_sum_1[["( -", ")"], [0, -1, -2]]),
                FadeOut(sine_sum_1["-", [0, 3]]),
                TrsP(
                    sine_sum_1, sine_sum,
                    ["sin ( alpha ", ") = sin alpha cos",
                     "cos alpha sin ", "beta"], ..., [0] * 3 + [None],
                    offset=1/6
                ),
                FadeIn(sine_sum["+", [0, 1]]),
                offset=1/2
            )
            self.complete(sine_sum)

            cos_sum_1 = MyTypst(cd_str.replace("beta", "( - beta)"))
            cos_sum_1.points.shift(np.array([
                (get_equal(cos_diff) - get_equal(cos_sum_1))[0], 0.6, 0,]))
            cs_str = "cos (alpha + beta) = cos alpha cos beta - sin alpha sin beta"
            cos_sum = MyTypst(cs_str)
            cos_sum.points.shift(np.array([
                (get_equal(cos_diff) - get_equal(cos_sum))[0], 0.6, 0,]))

            self.play(FadeIn(cos_diff))
            self.play(TrsP(
                cos_diff.copy(), cos_sum_1,  "beta", "( - beta)", [None],
                gapless=True, from_moved=True, offset=1/6
            ))

            self.play(
                Indicate(self.formula[96:108]),
                FadeOut(cos_sum_1[["( -", ")"], [0, -1, -2]]),
                FadeOut(cos_sum_1[[*"+-"], 0]),
                TrsP(
                    cos_sum_1, cos_sum,
                    ["cos ( alpha ", ") = cos alpha cos", "sin alpha sin "]
                    + ["beta"] * 3, [0]*4 + [1, 2],
                    offset=1/6
                ),
                FadeIn(cos_sum[*"+-"]),
                offset=1/2
            )
            self.forward()
            self.play(*[Transform(i, self.formula[j:k]) for i, (j, k) in zip(
                [sine_sum, sine_diff, cos_diff, cos_sum, sine_negative],
                [(26, 52)] * 2 + [(52, 78)] * 2 + [(108, 121)]
            )], offset=1/3)
            self.forward()

############ FORMULA VI, VII ############
        if 4 in play_parts:

            cs_sd_str = [
                "sin (alpha + beta) = sin alpha cos beta + cos alpha sin beta",
                "sin (alpha - beta) = sin alpha cos beta - cos alpha sin beta",
                "cos (alpha + beta) = cos alpha cos beta - sin alpha sin beta",
                "cos (alpha - beta) = cos alpha cos beta + sin alpha sin beta",
            ]
            cs_sd = Group(*map(MyTypst, cs_sd_str))
            cs_sd.points.arrange(5.4 * DOWN).shift(
                0.6 * UP - get_equal(cs_sd[0])[0] * RIGHT
            )
            self.play(FadeIn(cs_sd))

            cs_sd_1_str = [
                "sin (alpha + beta) + sin (alpha - beta) = 2 sin alpha cos beta",
                "sin (alpha + beta) - sin (alpha - beta) = 2 cos alpha sin beta",
                "cos (alpha + beta) + cos (alpha - beta) = 2 cos alpha cos beta",
                "cos (alpha + beta) - cos (alpha - beta) = -2 sin alpha sin beta",
            ]
            cs_sd_1 = Group(*map(MyTypst, cs_sd_1_str))
            cs_sd_1.points.arrange(5.4 * DOWN).shift(
                0 * UP - get_equal(cs_sd_1[0])[0] * RIGHT
            )

            big_plus = Typst("+")
            big_plus.points.scale(5).next_to(cs_sd_1, LEFT)
            big_minus = Typst("-")
            big_minus.points.scale(5).next_to(cs_sd_1, LEFT)

            for sign, fi, fj in ((big_plus, 0, 2), (big_minus, 1, 3)):
                # fi = formula index
                self.play(FadeIn(sign))
                self.play(AnimGroup(
                    *[TrsP(
                        cs_sd[i], cs_sd_1[j],
                        [get_between(cs_sd_str[i], "", "="),
                         get_between(cs_sd_1_str[j], "2", ""), ],
                        trs_keywords={"hide_src": False}
                    )for i, j in zip(range(4), [fi, fi, fj, fj])],
                    offset=1/6,
                ))
                self.forward(1/60)
                self.complete(cs_sd_1[fi], cs_sd_1[fj])
                self.play(FadeOut(sign))

            self.forward()

            substitution = {
                "(alpha + beta)": "theta",
                "(alpha - beta)": "phi",
                "alpha": "(theta + phi)/2",
                "beta": "(theta - phi)/2",
            }

            cs_factor_str = [reduce(
                lambda s, t: s.replace(t[0], t[1]),
                substitution.items(), s
            ) for s in cs_sd_1_str]
            cs_factor = Group(*(MyTypst(i, fill_color=GOLD)
                                for i in cs_factor_str))
            cs_factor.points.arrange(3.8 * DOWN).shift(
                0.52 * DOWN - get_equal(cs_factor[0])[0] * RIGHT
            )
            self.play(*[TrsP(
                cs_sd_1[i], cs_factor[i], substitution.keys(),
                substitution.values(), [0, 0, -1, -1], [0],
                gapless=True,
                trs_keywords={"hide_src": False},
                offset=1/6,
            )for i in range(4)], offset=1/2)

            cs_split_strs = [np.array([
                get_between(s, "2", ""),
                get_between(s, "(alpha - beta)", "2"), " 1/2 [",
                get_between(s, "", "="), "]"
            ]) for s in cs_sd_1_str]
            cs_split = Group(*[MyTypst("".join(s), fill_color=TEAL)
                               for s in cs_split_strs])
            cs_split.points.arrange(3.8 * DOWN).shift(
                0.08 * UP - get_equal(cs_split[0])[0] * RIGHT
            )

            self.forward()
            self.play(
                FadeOut(Group(*[i["2"] for i in cs_sd_1])))
            self.play(*[TrsP(
                cs_sd_1[i], cs_split[i], cs_split_strs[i][[0, 1, 3]],
                offset=1/6,
            )for i in range(4)], offset=1/2)

            self.forward(1/60)
            self.complete(*cs_split)
            self.forward()

            self.play(FadeOut(cs_sd))
            self.play(
                *[Transform(i, self.formula[140:171]) for i in cs_split],
                Write(self.formula[171:174]),
                offset=1/6
            )

            th_str = '''tan (theta + phi)/2 = 
(sin theta + sin phi)/(cos theta + cos phi) = 
-(cos theta - cos phi)/(sin theta - sin phi)'''
            tan_half = MyTypst(th_str)
            tan_half.points.shift(0.5 * DOWN)

            self.play(Write(tan_half["tan"]))
            self.play(*[TrsP(
                cs_factor[i], tan_half,
                ["(theta + phi)/2",
                 get_between(cs_factor_str[i], "", "=")],
                trs_keywords={"hide_src": False},
            )for i in [0, 2, 1, 3]], offset=1)

            self.complete(tan_half)

            self.forward()

            self.play(
                *[Transform(i, self.formula[174:201]) for i in cs_factor],
                Transform(tan_half, self.formula[204:253]),
                Write(self.formula[201:204]),
                Write(self.formula[253:256]),
                offset=1/6
            )
        return

    def others(self):
        return

############ FORMULA IX & X ############

        Cnt1 = Typst("cos n theta cos theta = 1/2 ( \
                     cos (n theta + theta) + \
                     cos (n theta - theta))")
        Cnt1.points.shift(3*UP+RIGHT)
        Cnt2 = Typst("cos (n+1) theta = 2 cos n theta cos theta \
                     - cos (n-1) theta")
        Cnt2.points.shift(2*UP+RIGHT)
        Cnt3 = Typst("cos 2 theta = 2 cos^2 theta - cos 0 theta")
        Cnt3.points.shift(1*UP+RIGHT)
        Cnt4 = Typst("T_(n+1) = 2 T_n id - T_(n-1)")
        Cnt4.points.shift(0.5*UP+RIGHT)
        Cnt5 = Typst("cos n theta = T_n (cos theta)")

        self.play(TrsS(
            Cp2s, Cnt1,
            [0, 3, 4, 17, 18, 26, 27, 31],
            [0, 3, 5, 18, 20, 28, 30, 34],
        ))
        self.play(TrsS(
            Cnt1.copy(), Cnt2,
            [[14, 19], [20, 23], [19, 20], [9, 10], [10, 13],
             [0, 9], [23, 24], [24, 29], [30, 33], [29, 30]],
            [0, 5, 8, 9, 10, 11,
             20, 21, 26, 29, 30],
            offset=1/5
        ))
        self.play(TrsS(
            Cnt2.copy(), Cnt3,
            [0, 3, 8, 14, 19, 21, 24, 29, 30],
            [0, 3, 4, 10, 11, 13, 16, 17, 18],
            offset=1/5
        ))
        self.play(Indicate(formula[58:64]))
        One = Typst("1")
        One.points.move_to(Cnt3[13:])
        self.play(Transform(Cnt3[13:], One))
        self.play(Transform(Cnt3[:13].add(One[0]), formula[215: 229]))
        self.play(TrsS(
            Cnt2.copy(), Cnt4,
            [0, 4, 7, 9,
             11, 14, 15, 16,
             20, 21, 25, 28, 30],
            [[0, 1], [1, 4], [0, 1], [4, 6],
             [6, 7], [7, 8], [6, 7], [8, 10],
             [10, 11], [11, 12], [12, 15], [11, 12]],
            offset=1/5
        ))
        Tedges = Typst("T_1=id, quad T_0=1")
        Tedges.points.shift(RIGHT)

        self.play(Succession(
            Write(Tedges[:6]),
            Indicate(formula[58:64]),
            Write(Tedges[6:])
        ))

        self.forward()
        self.show(Cnt1, Cnt2, Cnt3, Cnt4)
        self.forward(10)


if __name__ == "__main__":
    from janim.gui.anim_viewer import AnimViewer
    AnimViewer.views(MyTrigonal().build(), interact=True)
