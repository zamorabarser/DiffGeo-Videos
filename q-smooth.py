from fnmatch import translate
from this import d
from tkinter import E
from manim import *
from numpy import sqrt




class tism(Scene):
    def construct(self):
        t1 = Text("Smooth Curves", font_size=60)
        self.play(Write(t1))
        self.wait()        


class deriv(Scene):
    def construct(self):
        ti = Tex(r"Derivatives", font_size = 50, color = RED).to_edge(UL)
        v0 = ti.get_left()[0]
        exp = Tex(r"For a curve $\gamma : [a,b] \to \mathbb{R}^2$, its derivative is given by", font_size = 44).next_to(ti,DOWN)
        exp.shift((v0 - exp.get_left()[0])*RIGHT)
        f0 = Tex(r"$\gamma^{\prime}(t) : = \lim \dfrac{\gamma(t+h)-\gamma(t)}{h}$", font_size = 44).next_to(exp,DOWN).shift((0.1)*DOWN)
        vv = f0.get_center()[0]
        f0.shift(vv*LEFT)
        box = SurroundingRectangle(f0, buff = .2, color = YELLOW)
        f0h = Tex(r"$ h \to 0$", font_size = 22).next_to(exp,DOWN).shift((0.8)*DOWN+(0.675)*RIGHT)
        sp = ParametricFunction(lambda t : [(1.2)*t,(t-2)*t*(t+2)/4+t/2-1,0], t_range = [-2,2], color = BLUE)
        self.play(Write(ti), Write(exp), Write(f0), Write(f0h), Create(sp), Create(box))
        self.wait()
        dot0 = Dot(sp.point_from_proportion(0.31), color = ORANGE)
        dot1 = Dot(sp.point_from_proportion(0.75), color=GREEN)
        newt = Arrow(dot0.get_center() , dot1.get_center(), color= YELLOW, stroke_width=6, max_stroke_width_to_length_ratio=10, max_tip_length_to_length_ratio=0.15)
        newt.scale(1.1)
        self.play(Create(dot0), Create(dot1), Create(newt))
        self.wait(2)
        para = ValueTracker(0.75)
        dot1.add_updater(  lambda x: x.become(Dot(sp.point_from_proportion(para.get_value()), color=GREEN)))
        newt.add_updater(  lambda x: x.become(Arrow( dot0.get_center(), dot0.get_center() +  (0.47)*(sp.point_from_proportion(para.get_value())-dot0.get_center())/(para.get_value()-0.31)  , color= YELLOW, stroke_width=6, max_stroke_width_to_length_ratio=10, max_tip_length_to_length_ratio=0.15)))
        self.play(para.animate.set_value(0.32), run_time = 4)
        self.wait()
        


class chain(Scene):
    def construct(self):
        ti = Tex(r"Derivatives", font_size = 50, color = RED).to_edge(UL)
        v0 = ti.get_left()[0]
        exp = Tex(r"For a curve $\gamma : [a,b] \to \mathbb{R}^2$, its derivative is given by", font_size = 44).next_to(ti,DOWN)
        exp.shift((v0 - exp.get_left()[0])*RIGHT)
        f0 = Tex(r"$\gamma^{\prime}(t) : = \lim \dfrac{\gamma(t+h)-\gamma(t)}{h}$", font_size = 44).next_to(exp,DOWN).shift((0.1)*DOWN)
        vv = f0.get_center()[0]
        f0.shift(vv*LEFT)
        box = SurroundingRectangle(f0, buff = .2, color = YELLOW)
        f0h = Tex(r"$ h \to 0$", font_size = 22).next_to(exp,DOWN).shift((0.8)*DOWN+(0.675)*RIGHT)
        self.play(Write(ti), Write(exp), Write(f0), Write(f0h), Create(box))
        self.wait()
        d0 = Tex(r"$[c,d]$", font_size = 44).move_to(2*LEFT )
        d1 = Tex(r"$[a,b]$", font_size = 44).move_to(2*LEFT +2* DOWN)
        tar = Tex(r"$\mathbb{R}^2$", font_size = 44).move_to(2*RIGHT+DOWN)
        phiar = Arrow([-2,-0.2,0], [-2,-1.8,0],  tip_length = 0.2)
        phi = Tex(r"$\phi$", font_size = 44, ).next_to(phiar,LEFT)
        g1ar = Arrow([-1.5, - 1.85, 0], [1.5,-1.15, 0], tip_length = 0.2)
        g1 = Tex(r"$\gamma_1$", font_size = 44).next_to(g1ar, DOWN)
        g2ar = Arrow([-1.5,  -0.15, 0], [1.5,-0.85, 0],  tip_length = 0.2)
        g2 = Tex(r"$\gamma_2$", font_size = 44).next_to(g2ar, UP)
        self.play(Write(d0), Write(d1), Write(tar),Write(g1) , Write(g2), Create(g1ar), Create(g2ar), Create(phiar), Write(phi))
        cr = Tex(r"$\gamma_2^{\prime}(t) =$", r" $ (\gamma_1 \circ \phi ) ^{\prime} (t) = $", r"$\gamma_1^{\prime}(\phi (t) ) (\phi^{\prime}(t))$.", font_size = 44).shift((3.2)*DOWN)
        self.play(Write(cr[0]))
        self.wait()
        self.play(Write(cr[1]))
        self.wait()
        self.play(Write(cr[2]))
        self.wait()




class indep(Scene):
    def construct(self):
        le = Tex(r"Lemma", font_size = 50, color = RED).to_edge(UL)
        v0 = le.get_left()[0]
        lest = Tex(r"For a $C^1$ curve $\gamma : [a,b] \to \mathbb{R}^2$,", font_size = 40).next_to(le,DOWN)
        lest.shift((v0 - lest.get_left()[0])*RIGHT)
        f0 = Tex(r"length$(\gamma ) = \int_a^b \vert \gamma^{\prime} (t) \vert $d$t$", font_size = 40).next_to(lest,DOWN).shift((0.1)*DOWN)
        vv = f0.get_center()[0]
        f0.shift(vv*LEFT)
        box = SurroundingRectangle(f0, buff = .2, color = YELLOW)
        self.play(Write(le), Write(lest), Write(f0), Create(box))
        self.wait()
        d0 = Tex(r"$[c,d]$", font_size = 44).move_to(2*LEFT )
        d1 = Tex(r"$[a,b]$", font_size = 44).move_to(2*LEFT +2* DOWN)
        tar = Tex(r"$\mathbb{R}^2$", font_size = 44).move_to((1.9)*RIGHT+DOWN)
        phiar = Arrow([-2,-0.2,0], [-2,-1.8,0],  tip_length = 0.2)
        phi = Tex(r"$\phi$", font_size = 44, ).next_to(phiar,LEFT)
        g1ar = Arrow([-1.5, - 1.85, 0], [1.5,-1.15, 0], tip_length = 0.2)
        g1 = Tex(r"$\gamma_1$", font_size = 44).next_to(g1ar, DOWN)
        g2ar = Arrow([-1.5,  -0.15, 0], [1.5,-0.85, 0],  tip_length = 0.2)
        g2 = Tex(r"$\gamma_2$", font_size = 44).next_to(g2ar, UP)
        self.play(Write(d0), Write(d1), Write(tar),Write(g1) , Write(g2), Create(g1ar), Create(g2ar), Create(phiar), Write(phi))
        comd = Group(d0, d1, tar, phiar, phi, g1, g1ar, g2, g2ar)
        self.play(comd.animate.shift(3.7*LEFT))
        opt1 = Tex(r"$\int_c^d \vert \gamma_2^{\prime} (t) \vert $d$t$ = ", font_size = 40).next_to(f0,DOWN).shift(DOWN+(0.6)*RIGHT)
        h0 = opt1.get_left()[0]
        self.play(Write(opt1))
        opt2 = Tex(r"$\int_c^d  \vert   \gamma_1 ^{\prime} ( \phi (t)) \phi ^{\prime} (t) \vert $d$t$ =", font_size = 40).next_to(opt1,DOWN)
        opt2.shift((opt2.get_left()[0] - h0 )*LEFT )
        self.play(Write(opt2))
        u = Tex(r"$u = \phi (t) $", font_size = 40).next_to(opt2,RIGHT).shift((0.6)*RIGHT+(0.2)*UP)
        ul = u.get_left()[0]
        du = Tex(r"$du = \phi ^{\prime} (t) dt $", font_size = 40).next_to(u,DOWN)
        du.shift((du.get_left()[0]-ul)*LEFT + (0.15)*UP)
        self.play(Write(u), Write(du))
        opt3 = Tex(r"$\int_a^b  \vert   \gamma_1 ^{\prime} (u)  \vert $d$u$", font_size = 40).next_to(opt2,DOWN)
        opt3.shift((opt3.get_left()[0] - h0 )*LEFT )
        self.play(Write(opt3))
        self.wait()









class lem(Scene):
    def construct(self):
        le = Tex(r"Lemma", font_size = 50, color = RED).to_edge(UL)
        v0 = le.get_left()[0]
        lest = Tex(r"For a $C^1$ curve $\gamma : [a,b] \to \mathbb{R}^2$,", font_size = 40).next_to(le,DOWN)
        lest.shift((v0 - lest.get_left()[0])*RIGHT)
        f0 = Tex(r"length$(\gamma ) = \int_a^b \vert \gamma^{\prime} (t) \vert $d$t$", font_size = 40).next_to(lest,DOWN).shift((0.1)*DOWN)
        vv = f0.get_center()[0]
        f0.shift(vv*LEFT)
        box = SurroundingRectangle(f0, buff = .2, color = YELLOW)
        self.play(FadeIn(le), FadeIn(lest), FadeIn(f0), FadeIn(box))
        self.wait()
        coro = Tex(r"Corollary", font_size = 50, color = RED).next_to(f0, DOWN)
        coro.shift((v0 - coro.get_left()[0])*RIGHT)
        v0 = le.get_left()[0]
        alp = Tex(r"A $C^1$ curve $\gamma : [a,b] \to \mathbb{R}^2$ is a ", r" parametrization by arc length", font_size = 40).next_to(coro,DOWN)
        alp.shift((v0 - alp.get_left()[0])*RIGHT)
        ifandif = Tex(r"if and only if", font_size = 40).next_to(alp,DOWN)
        ifandif.shift((v0 - ifandif.get_left()[0])*RIGHT)
        f1 = Tex(r"$ \vert \gamma ^{\prime} (t)  \vert  = 1 $  ", r" for all $t \in [a,b]$.", font_size = 40).next_to(ifandif,DOWN).shift((0.1)*DOWN)
        vv = f1.get_center()[0]
        f1.shift(vv*LEFT)
        boxc0 = SurroundingRectangle(alp[1], buff = .1, color = GREEN_D)
        boxc1 = SurroundingRectangle(f1[0], buff = .1, color = GREEN_D)
        self.play(FadeIn(coro), FadeIn(alp), FadeIn(ifandif), FadeIn(f1))
        self.wait()
        self.play( Create(boxc0), Create(boxc1))
        self.wait()
        coro2 = Tex(r"Corollary", font_size = 50, color = RED).next_to(f1, DOWN)
        coro2.shift((v0 - coro2.get_left()[0])*RIGHT)
        salp = Tex(r"The arc length parametrization of a  smooth regular curve is", font_size = 40).next_to(coro2,DOWN)
        salp2 = Tex(r"smooth and regular", font_size = 40).next_to(salp,DOWN)
        salp2.shift((v0 - salp2.get_left()[0])*RIGHT)
        salp.shift((v0 - salp.get_left()[0])*RIGHT)
        self.play(Write(coro2), Write(salp), Write(salp2))
        self.wait()
        cor2 = Group(coro2, salp, salp2)
        self.play(FadeOut(le, lest, f0,box, coro, alp, ifandif, f1, boxc0, boxc1, cor2))
        #self.wait()
        #pro = Tex(r"Proof", font_size = 40, color = RED).next_to(salp2, DOWN)
        #pro.shift((v0 - pro.get_left()[0])*RIGHT)
        #self.play(Write(pro))
        #self.wait()
        #len = Tex(r"The length function of $\gamma$ is given by", font_size = 40).next_to(pro, DOWN)
        #len.shift((v0 - len.get_left()[0])*RIGHT)
        #for1 =  Tex(r"$\ell (t)  = \int_a ^t \vert \gamma^{\prime}(u) \vert du$", font_size = 40).next_to(len, DOWN)
        #for1.shift(for1.get_center()[0]*LEFT)
        #para = Tex(r"The arc length parametrization is given by", font_size = 40).next_to(for1, DOWN)
        #para.shift((v0 - para.get_left()[0])*RIGHT)
        #for2 =  Tex(r"$\alpha (s)  = \gamma ( \ell^{-1} (s) ) $", font_size = 40).next_to(para, DOWN)
        #for2.shift(for2.get_center()[0]*LEFT)
        #deri = Tex(r"Its derivative is given by", font_size = 40).next_to(for2, DOWN)
        #deri.shift((v0 - deri.get_left()[0])*RIGHT)
        #for3 =  Tex(r"$\alpha ^{\prime} (s)  = \gamma ^{\prime}( \ell^{-1} (s) ) / \ell^{\prime} ( \ell^{-1} (s) ) =  \gamma ^{\prime}( \ell^{-1} (s) ) / \vert \gamma^{\prime } (  \ell^{-1} (s) ) \vert $", font_size = 40).next_to(deri, DOWN)
        #for3.shift(for3.get_center()[0]*LEFT)


class ugex(Scene):
    def construct(self):
        ex = Tex(r"Exercise", font_size = 36, color = RED).to_edge(UL)
        v0 = ex.get_left()[0]
        lest = Tex(r"If $\gamma : [0,2 \sqrt{2} + \pi /\sqrt{2} ] \to \mathbb{R}^2$ is the arc length parametrization, then", font_size = 36).next_to(ex,DOWN)
        lest.shift((v0 - lest.get_left()[0])*RIGHT)
        axes = Axes(x_range = [-2,2,1], y_range = [-0.2,1.7,1], x_axis_config={ "numbers_to_include": np.arange(-2, 2.1, 1)},y_axis_config={ "numbers_to_include": np.arange(0, 2, 1)}, tips=False)
        axes.scale(0.335*RIGHT+0.325*UP).shift((1.15)*UP+0.08*RIGHT)
        l1 = ParametricFunction(lambda t : [-1 + t + PI/4, -3.1+t+PI/4,0], t_range = [-PI/4-1,-PI/4 ])#Line([-2,-4,0], [-1,-3,0])
        l2 = ParametricFunction(lambda t : [sqrt(2)* np.sin(t), sqrt(2)* np.cos(t)-4.1,0], t_range = [-PI/4, PI/4])
        l3 = ParametricFunction(lambda t : [ t + 1  - PI/4, - t + PI/4 -3.1,0], t_range = [PI/4,PI/4 +1 ])# Line([1,-3, 0], [2,-4,0])
        l = l1.copy().append_points(l2.points).append_points(l3.points)
        l.set_color(BLUE).shift(4.5*UP)
        self.play(Write(ex), Write(lest))#        self.play(Create(l1), Create(l2), Create(l3))
        self.wait()
        self.play(Create(l), Create(axes))#        self.play(Create(l1), Create(l2), Create(l3))
        self.wait()
        gam = Tex(r"$\gamma ^{\prime} ( t) = $", font_size = 36).to_edge(UL).shift(DOWN + (0.25)*RIGHT)
        gam1 = Tex(r"$ \left( 1/\sqrt{2}, 1/ \sqrt{2} \right) , $",  font_size = 36 ).next_to(gam, (0.85)*RIGHT).shift((0.9)*RIGHT+(0.7)*UP)
        pos = Tex(r"$\left( \sin \left( \frac{t }{\sqrt{2}} + \frac{\pi - 4 }{4}  \right) , \cos  \left( \frac{t }{\sqrt{2}} + \frac{\pi - 4 }{4}   \right), $" ,  font_size = 36).next_to(gam1, DOWN).shift(0.15*LEFT+0.1*DOWN)
        zer = Tex(r"$(1 / \sqrt{2}, -1 / \sqrt{2} ), $",  font_size = 36).next_to(pos, DOWN).shift(0.15*LEFT+0.1*DOWN)
        tp =  Tex(r"$ t \in [0,\sqrt{2}]  $", font_size = 36).next_to(gam1, RIGHT).shift((2.3)*RIGHT)
        tn = Tex(r"$ t \in [\sqrt{2}, \sqrt{2} + \pi /\sqrt{2} ] $", font_size = 36).next_to(tp,  DOWN).shift(0.15*DOWN +(0.1)*RIGHT)
        tz = Tex(r"$ t \in [\sqrt{2}  + \pi /\sqrt{2} , 2 \sqrt{2} + \pi /\sqrt{2} ] $", font_size = 36).next_to(tn, DOWN).shift(0.15*DOWN +(0.1)*RIGHT)
        key = Tex(r"$\{$", font_size = 180).shift((4.35)*LEFT + (2.1)*UP).scale(0.5*RIGHT + (1.4)*UP)
        key.shift(0.2*LEFT)
        gam.shift(0.2*DOWN+1.2*LEFT)
        key.shift(1.5*LEFT)
        eq = Group(gam, pos, tp, tn, key, gam1, zer, tz).scale(0.8).shift(3.6*DOWN + 1.7*RIGHT)
        self.play(FadeIn(eq))
        self.wait()
        box = SurroundingRectangle(eq, buff = .25, color = BLUE)
        self.play(Create(box))
        self.wait()


class length(Scene):
    def construct(self):
        le = Tex(r"Lemma", font_size = 50, color = RED).to_edge(UL)
        v0 = le.get_left()[0]
        lest = Tex(r"For a $C^1$ curve $\gamma : [a,b] \to \mathbb{R}^2$,", font_size = 40).next_to(le,DOWN)
        lest.shift((v0 - lest.get_left()[0])*RIGHT)
        f0 = Tex(r"length$(\gamma ) = \int_a^b \vert \gamma^{\prime} (t) \vert $d$t$", font_size = 40).next_to(lest,DOWN).shift((0.1)*DOWN)
        vv = f0.get_center()[0]
        f0.shift(vv*LEFT)
        box = SurroundingRectangle(f0, buff = .2, color = YELLOW)
        self.play(Write(le), Write(lest), Write(f0), Create(box))
        self.wait()
        pr = Tex(r"Proof:", font_size = 50, color = RED).next_to(f0, DOWN)
        pr.shift((v0 - pr.get_left()[0])*RIGHT)
        self.play(Write(pr))
        self.wait()
        begin = Tex(r"First prove:", font_size = 40).next_to(pr, DOWN)
        begin.shift((v0 - begin.get_left()[0])*RIGHT)
        f00 =  Tex(r"$ \int_a^b \vert \gamma^{\prime} (t) \vert $d$t \geq \vert \gamma ( b ) - \gamma (a ) \vert $", font_size = 40).next_to(begin,DOWN)
        f00.shift(f00.get_center()[0]*LEFT)
        self.play(Write(begin), Write(f00))
        self.wait()
        let = Tex(r"Let $\eta : [a,b] \to \mathbb{R}$ be given by", font_size = 40).next_to(f00, DOWN)
        let.shift((v0 - let.get_left()[0])*RIGHT)
        etad =  Tex(r"$ \eta ( t ) : = (\gamma (t) -\gamma (a) ) \cdot (\gamma ( b) - \gamma (a)) $", font_size = 40).next_to(let,DOWN)
        etad.shift(etad.get_center()[0]*LEFT)
        etadp = Tex(r"$ \eta ^{\prime } ( t )  = (\gamma^{\prime} (t) ) \cdot (\gamma ( b) - \gamma (a))$", r"$ \leq \vert \gamma^{\prime}(t) \vert \vert \gamma ( b ) - \gamma ( a) \vert  $", font_size = 40).next_to(etad,DOWN)
        etadp.shift(etadp.get_center()[0]*LEFT)
        etaf = MathTex(r"{{\vert \gamma (b) - \gamma (a ) \vert ^2}} {{= \eta ( b ) - \eta ( a )}} {{=  \int _a ^b \eta^{\prime}(t) \mathrm{d} t }} {{\leq}} {{  \vert \gamma ( b ) - \gamma ( a)  \vert }} {{ \int _a ^b \vert \gamma^{\prime}(t) \vert \mathrm{d}t}}", font_size = 40).next_to(etadp,DOWN)
        #etaf = MathTex(r"gbga", r"eta", r"inteta", r"intgt", font_size = 40).next_to(etadp,DOWN)
        etaf.shift(etaf.get_center()[0]*LEFT)
        simply0 = MathTex(r"{{\vert \gamma (b) - \gamma (a ) \vert ^2}} {{\leq}} {{  \vert \gamma ( b ) - \gamma ( a)  \vert }} {{ \int _a ^b \vert \gamma^{\prime}(t) \vert \mathrm{d}t}} ", font_size = 40)
        simply1 = MathTex(r"{{  \vert \gamma ( b ) - \gamma ( a)  \vert }} {{\leq}}   {{ \int _a ^b \vert \gamma^{\prime}(t) \vert \mathrm{d}t}} ", font_size = 40)
        etadp.shift(etadp.get_center()[0]*LEFT)
        self.play(Write(let), Write(etad))
        self.wait()
        self.play(Write(etadp))
        self.wait()
        self.play(Write(etaf[0]), Write(etaf[2]))
        self.wait()
        self.play(Write(etaf[4]))
        self.wait()
        self.play(Write(etaf[6]), Write(etaf[8]), Write(etaf[10]))
        self.wait()
        bulk = Group(begin, f00, let, etad, etadp)
        self.play(FadeOut(bulk), FadeOut(etaf[2]), FadeOut(etaf[4]), Transform(etaf[0], simply0[0]), Transform(etaf[6], simply0[2]), Transform(etaf[8], simply0[4]), Transform(etaf[10], simply0[6]))
        self.wait()
        self.play(Transform(etaf[8], simply1[0]), Transform(etaf[6], simply1[2]), Transform(etaf[10], simply1[4]), FadeOut(etaf[0]))
        self.wait()
        bf = SurroundingRectangle(simply1, buff = .15, color = RED)
        self.play(Create(bf))
        self.wait()
        self.remove(etaf[8], etaf[0], etaf[6], etaf[10])
        self.add(simply1)
        self.play( FadeOut(pr), FadeOut(box), FadeOut(bf), FadeOut(le), FadeOut(lest), FadeOut(f0), FadeOut(simply1))



class p1(Scene):
    def construct(self):
        ps = Tex(r"If $\gamma : [a,b] \to \mathbb{R}^2$ is $C^1$, then for any partition:", font_size = 40).to_edge(UL)
        v0 = ps.get_left()[0]
        eq1 = Tex(r"$ \int_a^b \vert \gamma ^{\prime } (t) \vert \mathrm{d}t = \sum_{j=1}^k \int_{t_{j-1}}^{t_j} \vert \gamma ^{\prime } (t) \vert \mathrm{d}t \geq \sum_{j=1}^k \vert \gamma (t_j ) - \gamma (t_{j-1}) \vert $", font_size = 40).next_to(ps, DOWN).shift((0.1)*DOWN)
        eq1.shift(eq1.get_center()[0]*LEFT)
        in1 = Tex(r"$\int_a^b \vert \gamma ^{\prime } (t) \vert \mathrm{d}t  \geq $ length$(\gamma)$.", font_size = 40).next_to(eq1, DOWN).shift((0.1)*DOWN)
        in1.shift(in1.get_center()[0]*LEFT)
        b1 = SurroundingRectangle(in1, buff = .15, color = LIGHT_PINK)
        pc = Tex(r"Fix $\varepsilon > 0$ and choose $a = t_0 < \ldots < t_k =b$ such that for each $j$,", font_size = 40).next_to(in1,DOWN).shift((0.1)*DOWN)
        pc.shift((v0 - pc.get_left()[0])*RIGHT)
        opt1 = Tex(r"$\bullet$ ", r"$ \vert \gamma^{\prime} (t) \vert \leq  \varepsilon / (b-a) $ for all $t \in [t_{j-1}, t_j ]$.", font_size=40).next_to(pc,DOWN).shift((0.1)*DOWN)
        opt1.shift((v0 - opt1.get_left()[0])*RIGHT)
        bo1 = SurroundingRectangle(opt1, buff = .1, color = BLUE)
        ortext = Tex(r"or", font_size=40).next_to(opt1,DOWN)
        ortext.shift((- ortext.get_center()[0])*RIGHT)
        opt2 = Tex(r"$\bullet$ ", r"$ \vert \gamma^{\prime} (t) \vert \geq \varepsilon /2 (b-a) $ for all $t \in [t_{j-1}, t_j ]$.", font_size=40).next_to(ortext,DOWN).shift((0.1)*DOWN)
        opt2.shift((v0 - opt2.get_left()[0])*RIGHT)
        bo2 = SurroundingRectangle(opt2, buff = .1, color = ORANGE)
        exersu = Tex(r"$  \int_a^b \vert \gamma ^{\prime } (t) \vert \mathrm{d}t  \leq $", r"$ \varepsilon $", r"$ + $", r"$  \sum_{j } \int_{t_{j-1}}^{t_j} \vert \gamma ^{\prime } (t) \vert \mathrm{d}t $", font_size = 40).next_to(opt2,DOWN).shift((0.2)*DOWN)
        exersu.shift(exersu.get_center()[0]*LEFT)
        ex11 = Group(exersu[1], exersu[2], exersu[3])
        ex11.shift(0.1*RIGHT)
        exersu[2].shift(0.1*RIGHT)
        exersu[3].shift(0.2*RIGHT)
        exersu.shift(0.2*LEFT)
        bb1 = SurroundingRectangle(exersu[1], buff = 2/30, color = BLUE)
        bb2 = SurroundingRectangle(exersu[3], buff = .1, color = ORANGE)
        self.play(Write(ps), Write(eq1))
        self.wait()
        self.play(Write(in1), Create(b1))
        self.wait()
        self.play(Write(pc), Write(opt1), Write(opt2), Write(ortext))
        self.wait()
        self.play(Write(exersu))
        self.wait()
        self.play(Create(bo1), Create(bo2), Create(bb1), Create(bb2))
        self.wait()
        all = Group(ps, eq1,  in1, b1, pc, opt1, bo1, opt2, bo2, ortext)
        eqs2 = Group(exersu, bb1, bb2)
        self.play(FadeOut(all), eqs2.animate.shift(5*UP))
        exer = Tex(r"Exercise", font_size=40, color = RED).next_to(exersu,DOWN).shift(0.1*DOWN)
        exer.shift((v0 - exer.get_left()[0])*RIGHT)
        extx = Tex(r"If $\alpha : [c,d] \to \mathbb{R}^2$ is such that:", font_size=40).next_to(exer,DOWN).shift(0.1*DOWN)
        extx.shift((v0 - extx.get_left()[0])*RIGHT)
        c1 = Tex(r"$\bullet$ $ \vert \alpha (t) \vert \geq \varepsilon$ for all $t \in [c,d]$.", font_size=40).next_to(extx,DOWN).shift(0.1*DOWN)
        c2 = Tex(r"$\bullet$ $ \vert \alpha (t) - \alpha (s) \vert \leq \delta < \varepsilon$ for all $s,t \in [c,d]$.", font_size=40).next_to(c1,DOWN).shift(0.1*DOWN)
        c1.shift((v0 - c1.get_left()[0])*RIGHT)
        c2.shift((v0 - c2.get_left()[0])*RIGHT)
        thtx = Tex(r"Then", font_size=40).next_to(c2,DOWN).shift(0.1*DOWN)
        thtx.shift((v0 - thtx.get_left()[0])*RIGHT)
        con = Tex(r"$ \int_c^d \vert \alpha (t) \vert \mathrm{d}t \leq \dfrac{\varepsilon}{\varepsilon - \delta } \vert \int_c^d  \alpha (t)  \mathrm{d}t  \vert  $.", font_size=40).next_to(thtx,DOWN).shift(0.1*DOWN)
        con.shift((con.get_center()[0])*LEFT)
        self.play(Write(exer), Write(extx), Write(c1), Write(c2), Write(thtx), Write(con))
        self.wait()
        self.play( extx.animate.shift((1/15)*UP), c1.animate.shift((2/15)*UP), c2.animate.shift((1/5)*UP), thtx.animate.shift((4/15)*UP), con.animate.shift((1/3)*UP))
        use = Tex(r"Use this to prove that for each $j$,", font_size=40).next_to(con,DOWN)
        use.shift((v0 - use.get_left()[0])*RIGHT)
        final = Tex(r"$ \int_{t_{j-1}}^{t_j} \vert \gamma ^{\prime } (t) \vert \mathrm{d}t  = $ length$(\gamma _{[t_{j-1}, t_j]})$", font_size=40).next_to(use,DOWN)
        final.shift(final.get_center()[0]*LEFT)
        somb = SurroundingRectangle(final, buff = 0.1, color = ORANGE)
        self.play(Write(use), Write(final), Write(somb))
        self.wait()
        allagain = Group(exer, extx, c1, c2, thtx, con, use, final, somb)
        self.play(FadeOut(allagain), FadeOut(eqs2))
        

        














class pplan(Scene):
    def construct(self):
        prop = Tex(r"Theorem", r": If $\lambda : \{ Curves \} \to \mathbb{R}$ satisfies:", font_size=44).to_edge(UL)
        prop[0].set_color(LIGHT_PINK)
        v0 = prop.get_left()[0]
        p0 = Tex(r"$\bullet$ ", r" $ \lambda (\gamma ) \geq \vert \gamma (b) - \gamma (a) \vert$ ", r", with equality if $\gamma$ is a line.", font_size=44).next_to(prop,DOWN).shift((0.5)*DOWN)
        p0.shift((v0 - p0.get_left()[0])*RIGHT)
        b0 = SurroundingRectangle(p0[1], buff = .1, color = RED)
        p1 = Tex(r"$\bullet$ ", r"$\lambda (\gamma_1 ) =  \lambda  (\gamma_2)$", r" when $\gamma_1, \gamma_2$ are reparametrizations of each other.", font_size=44).next_to(p0,DOWN).shift((0.5)*DOWN)
        p1.shift((v0 - p1.get_left()[0])*RIGHT)
        b1 = SurroundingRectangle(p1[1], buff = .1, color = YELLOW)
        p2 = Tex(r"$\bullet$ ", r" $ \lambda (\gamma_1 \ast \gamma_2 ) = \lambda  (\gamma_1) + \lambda  (\gamma _2)$ ", r" for concatenations.", font_size=44).next_to(p1,DOWN).shift((0.5)*DOWN)
        p2.shift((v0 - p2.get_left()[0])*RIGHT)
        p2[1].shift(0.1*RIGHT)
        p2[2].shift(0.2*RIGHT)
        b2 = SurroundingRectangle(p2[1], buff = .1, color = GREEN)
        b21 = SurroundingRectangle(p2[1], buff = .2, color = BLUE)
        p3 = Tex(r"$\bullet$ ", r"$ \lambda  (T \circ \gamma ) = \lambda  (\gamma) $ ", r" for $T$ any isometry.", font_size=44).next_to(p2,DOWN).shift((0.5)*DOWN)
        p3.shift((v0 - p3.get_left()[0])*RIGHT)
        b3 = SurroundingRectangle(p3[1], buff = .1, color = BLUE)
        p4 = Tex(r"$\bullet$ If $\gamma_n \to \gamma$ pointwise, ", r" $ \lambda  ( \gamma ) \leq  \liminf \lambda  (\gamma_n) $ ", r" .", font_size=44).next_to(p3,DOWN).shift((0.5)*DOWN)
        p4.shift((v0 - p4.get_left()[0])*RIGHT)
        b4 = SurroundingRectangle(p4[1], buff = .1, color = ORANGE)
        conc = Tex(r"Then ", r"$\lambda (\gamma ) = $ length$(\gamma)$ ", r" for all curves $\gamma \in$", r"$ \{ Curve \} $.")
        box = SurroundingRectangle(conc[1], buff = .1, color = LIGHT_PINK)
        self.play(Write(prop), Write(p0), Create(b0), Write(p1), Create(b1), Write(p2), Create(b2), Create(b21), Write(p3), Create(b3), Write(p4), Create(b4), Write(conc), Create(box) )
        self.wait()

class eqtest(Scene):
    def construct(self):
        e1 = MathTex(r"{{part0}} {{part1}} {{part2}} {{part3}}")
        self.play(Write(e1[0]), Write(e1[2]), Write(e1[4]))
        self.wait()