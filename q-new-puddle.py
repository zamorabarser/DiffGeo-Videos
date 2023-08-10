from this import d
from tkinter import E
from manim import *
from numpy import sqrt
import math



class title(Scene):
    def construct(self):
        t1 = Text("Moon in a puddle", font_size=60)
        self.play(Write(t1))
        self.wait()        


class moon(Scene):
    def construct(self):
        thm = Tex(r"Theorem", font_size=40, color= BLUE).to_edge(UL)
        thm.shift(0.3*DOWN)
        l = thm.get_left()[0]
        thm1 = Tex(r"Let $\gamma : [a,b] \to \mathbb{R}^2$ be a simple smooth regular closed curve.", font_size = 40).next_to(thm, DOWN)
        thm1.shift((thm1.get_left()[0]-l)*LEFT)
        thm2 = Tex(r"If $\kappa (t) \leq 1$ for all $t$, then $\gamma$ contains a unit disk in its interior.", font_size = 40).next_to(thm1, DOWN)
        thm2.shift((thm2.get_left()[0]-l)*LEFT)
        self.play(Write(thm), Write(thm1), Write(thm2))
        self.wait()
        pud = ParametricFunction(lambda t : [(0.7*np.cos( 3*t) + 1.7)*np.cos(t+0.1*np.cos( 3*t)) , (0.7*np.cos( 3*t) + 1.6)*np.sin(t+0.1*np.cos( 3*t))  -1.2 ,0], t_range = [0,TAU], color = BLUE)
        moo = Circle(radius = 0.7, color = YELLOW_A, fill_opacity = 0.3)
        moo.shift(1.25*DOWN+1.1*RIGHT)
        self.play(Create(pud), Create(moo))
        self.wait()




class tan(Scene):
    def construct(self):
        ti = Tex(r"Tangent curves", font_size=40, color= BLUE).shift(2.5*UP)
        coo =  Tex(r"Co-oriented", font_size=40).shift(2.8*DOWN+ 2*RIGHT)
        cou =Tex(r"Counter-oriented", font_size=40).shift(2.8*DOWN+ 2*RIGHT)
        p0 = Tex(r"$p_0 = \gamma_1(t_1)=  \gamma_2(t_2) $", font_size=40).shift(4*LEFT + UP)
        p1 = Tex(r"$ \gamma_1^{\prime}(t_1)  \vert \vert  \gamma_2^{\prime}(t_2) $", font_size=40).shift(4*LEFT +0.2*DOWN)
        cool = Line([-1.6+2,-3.1, 0], [1.6+2 ,-3.1,0], color = GREEN)
        coul = Line([-1.6+2,-3.1, 0],[1.6+2 ,-3.1,0], color = ORANGE)
        d0 = Dot([2,-0.5, 0], color = YELLOW)
        p0l = Line([-5.5,0.7, 0], [-2.5 ,0.7,0], color = YELLOW)
        self.play(Write(ti))
        g1 = ParametricFunction(lambda t : [t+2 , -0.2 * t**2 -0.5,0], t_range = [-3,3], color = BLUE)
        g2 = ParametricFunction(lambda t : [t+2 , 0.2 * t**2 + 0.03*t**3 -0.5 ,0], t_range = [-3,3], color = BLUE)
        g3 = ParametricFunction(lambda t : [-t +2, 0.2 * t**2 - 0.03*t**3 -0.5 ,0], t_range = [-3,3], color = BLUE)
        self.play(Create(g1), Create(g2), Write(p0), Create(p0l), Create(d0))
        self.wait()
        tang = DashedLine([-3+2,-0.5,0], [3+2,-0.5,0])
        self.play(Create(tang), Write(p1))
        self.wait()
        d1 = Dot([-3+2,-2.3, 0], color = GREEN)
        d2 = Dot([-3+2, 0.51, 0], color = GREEN)
        self.play(Create(d1), Create(d2), FadeOut(tang), Write(coo), Create(cool))
        self.wait()
        self.play(MoveAlongPath(d1,g1), MoveAlongPath(d2, g2), run_time = 2)
        self.wait()
        self.play(FadeOut(coo), FadeOut(d1), FadeOut(d2), FadeOut(cool))
        d1 = Dot([-3+2,-2.3, 0], color = ORANGE)
        d2 = Dot([3+2, 2.11, 0], color = ORANGE)
        self.play(Create(d1), Create(d2), Write(cou), Create(coul))
        self.wait()
        self.play(MoveAlongPath(d1,g1), MoveAlongPath(d2, g3), run_time = 2)
        self.wait()

        


class sup(Scene):
    def construct(self):
        ti = Tex(r"Supporting curves", font_size=40, color= BLUE).shift(2.5*UP)
        self.play(Write(ti))
        l1 = Tex(r"$\gamma_1$", font_size=40).shift(4.5*RIGHT + 1.5*DOWN)
        l2 = Tex(r"$\gamma_2$", font_size=40).shift(4.5*RIGHT+1.5*UP)
        p0 = Tex(r"$p_0 = \gamma_1(t_1)=  \gamma_2(t_2) $", font_size=40).shift(4*LEFT +1.6*UP)
        c1 = Tex(r"$\gamma_1(t_1-\varepsilon , t_1 + \varepsilon ) \subset \partial  R$", font_size=40).shift(1.8*DOWN)
        c2 = Tex(r"$\gamma_2(t_2-\varepsilon , t_2 + \varepsilon ) \subset R$", font_size=40).next_to(c1, DOWN).shift(0.2*DOWN)
        rl = Tex(r"$R$", font_size=40).shift(2*RIGHT+1.6*UP)
        rlb = Line([1.8, 1.3, 0], [2.2, 1.3, 0], color=PINK)
        p0b = Line([-5.3, 1.2, 0], [-2.7, 1.2, 0], color=RED)
        l1b = Line([4.2,-1.8, 0], [4.8, -1.8, 0], color = YELLOW)
        l2b = Line([4.2,1.2, 0], [4.8, 1.2, 0], color = TEAL)
        clb1 = Line([-2, -2.1, 0], [2,-2.1, 0], color=YELLOW)
        clb2 = Line([-2,-3, 0], [2, -3, 0], color=TEAL)
        self.wait()
        p = Dot([0,-0.37, 0], color = RED)
        ax = Axes(x_range=[-5, 5], y_range=[-4, 4])
        g1 = ax.plot(lambda x: -0.2 * x**2 -0.5, x_range=[-3, 3], color=YELLOW)
        g2 = ax.plot(lambda x: 0.2 * x**2 -0.5+ 0.03*x**3 , x_range=[-3, 3], color=TEAL)
        g11 = ax.plot(lambda x: -0.2 * x**2 -0.5, x_range=[-1.8, 1.8], color=YELLOW)
        g22 = ax.plot(lambda x: 0.2 * x**2 -0.5+ 0.03*x**3 , x_range=[-1.8, 1.8], color=TEAL)
        g3 = ax.plot(lambda x: 1 , x_range=[-3, 3], color=BLUE)
        self.play(Create(g1), Create(g2), Create(l1b), Create(l2b), Write(l1), Write(l2), Write(p0), Create(p), Create(p0b))
        r = ax.get_area(g1, [-2, 2], bounded_graph=g3, color=PINK, opacity=0.3)
        self.wait()
        self.play(Create(r), Create(rlb), Write(rl))
        self.wait()
        self.add(g11, g22)
        self.play(Write(c1), Write(c2), FadeOut(g1), FadeOut(g2), Create(clb1), Create(clb2))
        self.wait()
        tip1 = Line([2.2,-0.9,0], [2, -0.7,0], color = YELLOW)
        tip11 = Line([2.2,-0.9,0], [1.95, -0.85,0], color = YELLOW)
        tip2 = Line([2.2,0.25,0], [2.03, 0.1,0], color = TEAL)
        tip22 = Line([2.2,0.25,0], [1.95, 0.2,0], color = TEAL)
        self.play(Create(tip1), Create(tip11), Create(tip2), Create(tip22))
        self.wait()
        self.play(FadeOut(c1), FadeOut(c2), FadeOut(clb1), FadeOut(clb2), FadeOut(l1b), FadeOut(l2b), FadeOut(l1), FadeOut(l2), FadeOut(p0), FadeOut(p0b), FadeOut(r), FadeOut(rl), FadeOut(rlb))
        self.wait()




class gra(Scene):
    def construct(self):
        ti = Tex(r"Supporting curves", font_size=40, color= BLUE).shift(2.5*UP)
        l1 = Tex(r"$\gamma_2  = $ graph of $ f_2$", font_size=40).shift(4.2*LEFT + 0.5*UP)
        l2 = Tex(r"$\gamma_1  = $ graph of $ f_1$", font_size=40).shift(4.2*LEFT+0.5*DOWN)
        l12 = Tex(r"$ f_1 (x) \leq f_2 (x)$", font_size=40).shift(4.2*RIGHT)
        l1b = Line([-5.2,0.2, 0], [-3.2, 0.2, 0], color = TEAL)
        l2b = Line([-5.2,-0.8, 0], [-3.2, -0.8, 0], color = YELLOW)
        l12b = Line([3.5,-0.4, 0], [4.9, -0.4, 0], color = PINK)
        p = Dot([0,-0.37, 0], color = RED)
        ax = Axes(x_range=[-5, 5], y_range=[-4, 4])
        g11 = ax.plot(lambda x: -0.2 * x**2 -0.5, x_range=[-1.8, 1.8], color=YELLOW)
        g22 = ax.plot(lambda x: 0.2 * x**2 -0.5+ 0.03*x**3 , x_range=[-1.8, 1.8], color=TEAL)
        r = ax.get_area(g11, [-1.5, 1.5], bounded_graph=g22, color=PINK, opacity=0.3)
        tip1 = Line([2.2,-0.9,0], [2, -0.7,0], color = YELLOW)
        tip11 = Line([2.2,-0.9,0], [1.95, -0.85,0], color = YELLOW)
        tip2 = Line([2.2,0.25,0], [2.03, 0.1,0], color = TEAL)
        tip22 = Line([2.2,0.25,0], [1.95, 0.2,0], color = TEAL)
        self.add(g11, g22, tip1, tip11, tip2, tip22, p , ti)
        axx = Axes(x_range=[0,8], y_range=[0 ,8], x_length = 4, y_length = 2, tips=False)
        axx.shift(0.3*DOWN)
        self.wait()
        self.play(Create(axx), Write(l1), Write(l2), Create(l1b), Create(l2b))
        self.wait()
        self.play(Create(l12), Create(l12b), Create(r))
        self.wait()
        all = Group(l12, l12b, r, axx, l1, l2, l1b, l2b, g11, g22, tip1, tip11, tip2, tip22, p)
        self.play(all.animate.shift(0.7*UP))
        self.wait()     
        pro = Tex(r"Proposition", font_size=40, color= BLUE).shift(4.5*LEFT + 1.3*DOWN)
        l = pro.get_left()[0]
        pro1 = Tex(r"If $\gamma_1$ and $\gamma_2$ are tangent co-oriented at $\gamma_1(t_1) = \gamma_2(t_2)$, then", font_size=40).next_to(pro, DOWN)
        pro1.shift((pro1.get_left()[0]-l)*LEFT)
        pro2 = Tex(r"$\bullet$ if $\gamma_1$ supports $\gamma_2$ from the right, then $k_1(t_1) \leq k_2(t_2)$", font_size=40).next_to(pro1, DOWN)
        pro2.shift((pro2.get_left()[0]-l)*LEFT)
        pro3 = Tex(r"$\bullet$ if $k_1(t_1) < k_2(t_2)$, then $\gamma_1$ supports $\gamma_2$ from the right ", font_size=40).next_to(pro2, DOWN)
        pro3.shift((pro3.get_left()[0]-l)*LEFT)
        self.play(Write(pro), Write(pro1), Write(pro2), Write(pro3))
        self.wait()





class glo(Scene):
    def construct(self):
        g1 = ParametricFunction(lambda t : [ 4*np.cos(t) ,2* np.sin(t)  ,0], t_range = [0,TAU], color = BLUE)
        g2 = ParametricFunction(lambda t : [ (1-0.8*np.sin(t))* np.cos(t), (1-0.8*np.sin(t))*np.sin(t) -0.2 ,0], t_range = [0,TAU], color = ORANGE)
        g1l = Tex(r"$\gamma_1$", font_size=40).shift(2.6*RIGHT+2*UP)
        g2l = Tex(r"$\gamma_2$", font_size=40).shift(1.5*RIGHT)
        g1b = Line([2.4,1.7,0],[2.8,1.7,0], color = BLUE)
        g2b = Line([1.3,-0.3, 0], [1.7,-0.3, 0 ], color = ORANGE)
        self.play(Create(g1), Create(g2), Write(g1l), Write(g2l), Create(g1b), Create(g2b))
        self.wait()




class oscu(Scene):
    def construct(self):
        ti = Tex(r"Osculating circles", font_size = 40, color = BLUE).to_edge(UL).shift(DOWN)
        g = ParametricFunction(lambda t : [ 4*np.cos(t) ,2* np.sin(t)  ,0], t_range = [0,TAU], color = BLUE)
        kr = Tex(r"$r(t) = \frac{1}{\vert k(t) \vert }$", font_size = 40).to_edge(UR).shift(DOWN+0.8*LEFT)
        krb = SurroundingRectangle(kr,buff = .2, color= RED)
        def r(t):
            return sqrt(1 + 3* np.sin(t)**2)**3
        def n10(t):
            return -2*np.cos(t) + 6*np.cos(t)*(np.sin(t)**2)/(1 + 3* np.sin(t)**2)
        def n20(t):
            return -np.sin(t) -  3*np.sin(t)*(np.cos(t)**2)/(1 + 3* np.sin(t)**2)
        def n1(t):
            return n10(t)/sqrt( n10(t)**2 + n20(t)**2 )
        def n2(t):
            return n20(t)/sqrt( n10(t)**2 + n20(t)**2 )
        def c1(t):
            return 4*np.cos(t) + n1(t)*r(t)
        def c2(t):
            return 2*np.sin(t) + n2(t)*r(t)
        c = Dot([c1(PI/8), c2(PI/8), 0], color = PINK)  
        oc = Circle(radius = r(PI/8)).shift( c1(PI/8)*RIGHT + c2(PI/8)*UP)
        tl = DashedLine([4*np.cos(PI/8)+4*np.sin(PI/8) , 2*np.sin(PI/8) -2* np.cos(PI/8) , 0], [4*np.cos(PI/8)-4*np.sin(PI/8), 2*np.sin(PI/8)+2*np.cos(PI/8), 0])
        p = Dot([4*np.cos(PI/8), 2*np.sin(PI/8), 0 ], color = YELLOW)
        self.play(Write(ti), Create(g), Create(p))
        self.wait()
        self.play(Create(tl))
        self.wait()
        self.play(Create(oc), FadeOut(tl), Create(c), Write(kr), Create(krb))
        self.wait()
        s = ValueTracker(PI/8)
        oc.add_updater( lambda x : x.become(Circle(radius = r(s.get_value())).shift( c1(s.get_value())*RIGHT + c2(s.get_value())*UP))) 
        p.add_updater( lambda x : x.become(Dot( [4*np.cos(s.get_value()), 2*np.sin(s.get_value()), 0 ], color = YELLOW ))) 
        c.add_updater( lambda x : x.become(Dot( [c1(s.get_value()), c2(s.get_value()), 0 ], color = PINK ))) 
        self.play(s.animate.set_value(TAU), run_time = 15)

        
        
        

class lem(Scene):
    def construct(self):
        thm = Tex(r"Theorem", font_size=40, color= BLUE).to_edge(UL)
        thm.shift(RIGHT)
        l = thm.get_left()[0]
        thm1 = Tex(r"Let $\gamma : [a,b] \to \mathbb{R}^2$ be a simple smooth regular closed curve.", font_size = 40).next_to(thm, DOWN)
        thm1.shift((thm1.get_left()[0]-l)*LEFT)
        thm2 = Tex(r"If $\kappa (t) \leq 1$ for all $t$, then $\gamma$ contains a unit disk in its interior.", font_size = 40).next_to(thm1, DOWN)
        thm2.shift((thm2.get_left()[0]-l)*LEFT)
        self.play(Write(thm), Write(thm1), Write(thm2))
        self.wait()
        lem = Tex(r"Lemma", font_size=40, color= BLUE).to_edge(UL)
        lem.shift(2.2*DOWN+RIGHT)
        lem1 = Tex(r"Let $\gamma : [a,b] \to \mathbb{R}^2$ be a simple smooth regular loop.", font_size = 40).next_to(lem, DOWN)
        lem1.shift((lem1.get_left()[0]-l)*LEFT)
        lem2 = Tex(r"Then $\gamma$ supports one of its osculating circles from outside.", font_size = 40).next_to(lem1, DOWN)
        lem2.shift((lem2.get_left()[0]-l)*LEFT)
        self.play(Write(lem), Write(lem1), Write(lem2))
        self.wait()
        def l1(t):
            return 2.2*np.sin(t) - 0.3*np.cos(10*t)-1 + (0.2)*t*(PI-t)*np.cos(t)
        def l2(t):
            return -1.3*np.sin(2*t) + 0.15*np.sin(8*t)-2.5
        g = ParametricFunction(lambda t : [l1(t) ,l2(t) ,0], t_range = [0,PI], color = BLUE)
        oscc = Circle(radius = 0.1, fill_opacity = 0.4, color = YELLOW).shift( (l1(0.58)+0.14/sqrt(2) )*RIGHT + (l2(0.58)+0.1/sqrt(2))*UP )
        self.play(Create(g))
        self.wait()
        self.play(Create(oscc))
        self.wait()
        



class mproof(Scene):
    def construct(self):
        def nor(t):
            return math.exp(-2*t**2)
        def f(t):
            u = (t+2)/1000
            i = 1
            val = 0
            while i < 1001:
                val = val + nor( -2 + i*(t+2)/1000)*u
                i = i + 1
            return val
        def x0(t):
            return f(5*t-1)/3 + f(4*t-4)/3 + 2*f(5*t/2 - 5)/3 - 4*f(4*t - 11.7)/3
        def y0(t):
            return -f(6*t-1)/2 + f(4*t - 5/2)/3 - f(7*t - 7)/2 + f(7*t - 9)/3 - f(3*t - 6)/3 + 2* f(13*t/10 - 13/5)/3
        def x(t):
            return 4*x0(t)-3
        def y(t):
            return 4*y0(t) + 1
        lem = Tex(r"Lemma", font_size=40, color= BLUE).to_edge(UL).shift(RIGHT)
        l = lem.get_left()[0]
        lem1 = Tex(r"Let $\gamma : [a,b] \to \mathbb{R}^2$ be a simple smooth regular loop.", font_size = 40).next_to(lem, DOWN)
        lem1.shift((lem1.get_left()[0]-l)*LEFT)
        lem2 = Tex(r"Then $\gamma$ supports one of its osculating circles from outside.", font_size = 40).next_to(lem1, DOWN)
        lem2.shift((lem2.get_left()[0]-l)*LEFT)
        self.add(lem, lem1,lem2)
        self.wait()
        g = ParametricFunction(lambda t : [x(t) ,y(t) ,0], t_range = [0,3.5], color = BLUE)
        gl = Tex(r"$\gamma$", font_size = 40).shift(2*LEFT + 1.4*UP)
        glb = Line([-2.1, 1.2, 0], [-1.9, 1.2, 0], color = BLUE)
        self.play(Create(g), Write(gl), Create(glb))
        self.wait()
        p = Dot([x(0.21), y(0.21),0], color = YELLOW)
        pl = Tex(r"$p$", font_size = 40).shift(2.6*LEFT + (0.7)*DOWN)
        plb = Line([-2.5, -1, 0], [-2.7 , -1, 0], color = YELLOW)
        q = Dot([x(0.65), y(0.65),0], color = GREEN)
        self.play(Create(p), Write(pl), Create(plb))
        self.wait()
        ic = Circle(radius = 0)
        sl = Tex(r"$\sigma$", font_size = 40).shift(1.75*LEFT + (0.25)*UP)
        slb = Line([-1.85, 0.05, 0], [-1.65 , 0.05, 0], color = RED)
        ql = Tex(r"$q$", font_size = 40).shift(LEFT +0.5*DOWN)
        qlb = Line([-1.1, -0.8, 0], [-0.9 , -0.8, 0], color = GREEN)
        ic.shift(p.get_center())
        c = Dot(p.get_center(), radius = 0)
        self.add(ic)
        ic.add_updater( lambda x : x.become(Circle(radius = sqrt((c.get_center()[0]-p.get_center()[0] )**2 + (c.get_center()[1]-p.get_center()[1] )**2) ).shift( c.get_center() ))) 
        #c.add_updater( lambda x : x.become(Dot( [c1(s.get_value()), c2(s.get_value()), 0 ], color = PINK ))) 
        #self.play(s.animate.set_value(TAU), run_time = 15)
        self.play(c.animate.shift((0.24)*UP+(0.24)*(1.45)* RIGHT))
        self.wait()
        self.play(Create(q), Create(slb), Create(qlb), Write(ql), Write(sl))
        self.wait()
        cap = Tex(r"$\sigma =$ incircle of $p$ = largest circle tangent to $\gamma$ at $p$ supported by $\gamma$", font_size = 40).shift(2.7*DOWN)
        sll = slb.copy().shift(3.05*DOWN+4.15*LEFT)
        pll = plb.copy().shift(2*DOWN+0.6*LEFT)
        pll2 = pll.copy().shift(6.2*RIGHT)
        gll = glb.copy().shift(4.2*DOWN+4.1*RIGHT)
        gll2 = gll.copy().shift(3.8*RIGHT)
        self.play(Create(cap), Create(sll), Create(pll), Create(gll), Create(gll2), Create(pll2))
        self.wait()
        self.play(FadeOut(p), FadeOut(q), FadeOut(ic), FadeOut(ql), FadeOut(pl), FadeOut(plb), FadeOut(qlb), FadeOut(sl), FadeOut(slb) )
        self.wait()
        c = Dot([x(0.85), y(0.85),0], radius = 0)
        p1 = Dot([x(0.85), y(0.85),0], color = YELLOW)
        q1 = Dot([x(2.996), y(2.996),0], color = GREEN)
        gamma1 = ParametricFunction(lambda t : [x(t) ,y(t) ,0], t_range = [0.85,2.996], color = ORANGE)
        ic1 = Circle(radius = 0).shift(p1.get_center())
        ic1.add_updater( lambda x : x.become(Circle(radius = sqrt((c.get_center()[0]-p1.get_center()[0] )**2 + (c.get_center()[1]-p1.get_center()[1] )**2) ).shift( c.get_center() ))) 
        p1l = Tex(r"$p_1$", font_size = 40).shift((1.05)*LEFT + (0.4)*DOWN)
        gamma1l = Tex(r"$\gamma_1$", font_size = 40).next_to(p1l, RIGHT).shift(3*RIGHT+0.3*UP)
        gamma1b = Line([gamma1l.get_center()[0]-0.1,gamma1l.get_center()[1]-0.2, 0 ], [gamma1l.get_center()[0]+0.1, gamma1l.get_center()[1]-0.2,0], color = ORANGE)
        p1lb = Line([-1.15,-0.6, 0], [-0.95, -0.6, 0], color = YELLOW)
        q1l = Tex(r"$q_1$", font_size = 40).shift(LEFT +(1.4)* UP)
        q1lb = Line([-1.1,1.2, 0], [-0.9, 1.2, 0], color = GREEN)
        self.play(Create(p1), Create(ic1), Create(p1l), Create(p1lb))
        self.wait()
        self.play(c.animate.shift((0.458)*UP + 0.025*RIGHT))
        s1l = Tex(r"$\sigma_1$", font_size = 40).next_to(ic1, LEFT).shift(0.1*RIGHT)
        s1b = Line([s1l.get_center()[0]-0.1,s1l.get_center()[1]-0.2, 0 ], [s1l.get_center()[0]+0.1, s1l.get_center()[1]-0.2,0], color = RED)
        self.wait()
        self.play(Create(q1), Create(q1l), Create(q1lb), Create(s1l), Create(s1b))
        self.wait()
        self.play(Create(gamma1), Create(gamma1b), Write(gamma1l))
        self.wait()
        c2 = Dot([x(1.782), y(1.782),0], radius = 0)
        p2 = Dot([x(1.782), y(1.782),0], color = YELLOW)
        q2 = Dot([x(2.909), y(2.909),0], color = GREEN)
        ic2 = Circle(radius = 0).shift(p2.get_center())
        ic2.add_updater( lambda x : x.become(Circle(radius = sqrt((c2.get_center()[0]-p2.get_center()[0] )**2 + (c2.get_center()[1]-p2.get_center()[1] )**2) ).shift( c2.get_center() ))) 
        p2l = Tex(r"$p_2$", font_size = 40).shift((0.85)*RIGHT + (0.2)*DOWN)
        p2lb = Line([0.95,-0.4, 0], [0.75, -0.4, 0], color = YELLOW)
        q2l = Tex(r"$q_2$", font_size = 40).shift((0.8)*RIGHT +(1.4)* UP)
        q2lb = Line([0.7,1.2, 0], [0.9, 1.2, 0], color = GREEN)
        self.play(Create(p2), Create(ic2), Create(p2l), Create(p2lb))
        self.wait()
        self.play(c2.animate.shift((0.415)*UP+ 0.111*LEFT))
        s2l = Tex(r"$\sigma_2$", font_size = 40).next_to(ic2, LEFT).shift(0.1*RIGHT)
        s2b = Line([s2l.get_center()[0]-0.1,s2l.get_center()[1]-0.2, 0 ], [s2l.get_center()[0]+0.1, s2l.get_center()[1]-0.2,0], color = RED)
        self.wait()
        self.play(Create(q2), Create(q2l), Create(q2lb), Create(s2l), Create(s2b))
        self.wait()
        self.remove(gamma1, gamma1b, gamma1l)
        self.wait()
        c3 = Dot([x(2.05), y(2.05),0], radius = 0)
        p3 = Dot([x(2.05), y(2.05),0], color = YELLOW)
        q3 = Dot([x(2.818), y(2.818),0], color = GREEN)
        ic3 = Circle(radius = 0).shift(p3.get_center())
        ic3.add_updater( lambda x : x.become(Circle(radius = sqrt((c3.get_center()[0]-p3.get_center()[0] )**2 + (c3.get_center()[1]-p3.get_center()[1] )**2) ).shift( c3.get_center() ))) 
        p3l = Tex(r"$p_3$", font_size = 40).shift((2.35)*RIGHT + (0.2)*DOWN)
        p3lb = Line([2.25,-0.4, 0], [2.45, -0.4, 0], color = YELLOW)
        q3l = Tex(r"$q_3$", font_size = 40).shift((2.4)*RIGHT +(1.4)* UP)
        q3lb = Line([2.3,1.2, 0], [2.5, 1.2, 0], color = GREEN)
        self.play(Create(p3), Create(ic3), Create(p3l), Create(p3lb))
        self.wait()
        self.play(c3.animate.shift((0.4)*UP+ 0.02*RIGHT))
        s3l = Tex(r"$\sigma_3$", font_size = 40).next_to(ic3, LEFT).shift(0.1*RIGHT)
        s3b = Line([s3l.get_center()[0]-0.1,s3l.get_center()[1]-0.2, 0 ], [s3l.get_center()[0]+0.1, s3l.get_center()[1]-0.2,0], color = RED)
        self.wait()
        self.play(Create(q3), Create(q3l), Create(q3lb), Create(s3l), Create(s3b))
        self.wait()
        c4 = Dot([x(2.25), y(2.25),0], radius = 0)
        p4 = Dot([x(2.25), y(2.25),0], color = YELLOW)
        q4 = Dot([x(2.74), y(2.74),0], color = GREEN)
        ic4 = Circle(radius = 0).shift(p3.get_center())
        ic4.add_updater( lambda x : x.become(Circle(radius = sqrt((c4.get_center()[0]-p4.get_center()[0] )**2 + (c4.get_center()[1]-p4.get_center()[1] )**2) ).shift( c4.get_center() ))) 
        p4l = Tex(r"$p_4$", font_size = 40).shift((3.25)*RIGHT + (0.1)*DOWN)
        p4lb = Line([3.15,-0.3, 0], [3.35, -0.3, 0], color = YELLOW)
        q4l = Tex(r"$q_4$", font_size = 40).shift((3.2)*RIGHT +(1.4)* UP)
        q4lb = Line([3.3,1.2, 0], [3.1, 1.2, 0], color = GREEN)
        self.play(Create(p4), Create(ic4), Create(p4l), Create(p4lb))
        self.wait()
        self.play(c4.animate.shift((0.31)*UP+ (0.165)*LEFT))
        self.wait()
        self.play(Create(q4), Create(q4l), Create(q4lb))
        self.wait()
        p5 = Dot([x(2.35), y(2.35),0], color = YELLOW)
        q5 = Dot([x(2.69), y(2.69),0], color = GREEN)
        p5l = Tex(r"$p_5$", font_size = 40).shift((3.85)*RIGHT + (0.1)*UP)
        p5lb = Line([3.75,-0.1, 0], [3.95, -0.1, 0], color = YELLOW)
        q5l = Tex(r"$q_5$", font_size = 40).shift((3.8)*RIGHT +(1.3)* UP)
        q5lb = Line([3.7,1.1, 0], [3.9, 1.1, 0], color = GREEN)
        self.play(Create(p5), Create(q5), Create(p5lb), Create(q5lb), Write(p5l), Write(q5l))
        pinf = Dot([x(2.6), y(2.6),0], color = PINK)
        pq6 = Tex(r"$p_{\infty} = q_{\infty}$", font_size = 40).shift((4.7)*RIGHT + (0.8)*UP)
        p6b = Line([3.95, 0.6, 0], [5.45, 0.6, 0], color = PINK)
        self.play(Create(pinf), Write(pq6), Create(p6b))
        self.wait()
        



class mp2(Scene):
    def construct(self):
        def r(t):
            return 1.8 + np.sin(4*t) + 0.5*np.sin(2*t) - math.exp(- (2*t - PI/4)**2)-  math.exp( -(2*t - 5*PI/4)**2 ) /3 - math.exp(-(t - 3*PI/8)**2)/4
        def x(t):
            return np.cos(t-3*PI/8)*r(t)
        def y(t):
            return np.sin(t - 3*PI/8)*r(t)-0.5
        def r2(t):
            return 1.8 + np.sin(4*t)
        def x2(t):
            return np.cos(t-3*PI/8)*r2(t)
        def y2(t):
            return np.sin(t - 3*PI/8)*r2(t)-0.5
        cla = Tex(r"Claim: if $p_2$ is between $p_1$ and $q_1$, then so is $q_2$", font_size = 40).shift(2.5*UP)
        s1lab = Tex(r"$\sigma _1 =$ incircle of $p_1$", font_size = 40).shift(4.5*LEFT + UP)
        s2lab = Tex(r"$\sigma _2 =$ incircle of $p_2$", font_size = 40).next_to(s1lab, DOWN)
        s1b = Line([s1lab.get_left()[0], s1lab.get_center()[1]-0.3, 0], [s1lab.get_left()[0]+0.3, s1lab.get_center()[1]-0.3, 0], color = RED)
        s2b = Line([s2lab.get_left()[0], s2lab.get_center()[1]-0.3, 0], [s2lab.get_left()[0]+0.3, s2lab.get_center()[1]-0.3, 0], color = TEAL)
        lcp1 = Line([-2.4,2.2, 0], [-2.2, 2.2,0], color = YELLOW)
        lcp2 = Line([0.1,2.2, 0], [0.3, 2.2,0], color = YELLOW)
        lcq1 = Line([1.7,2.2, 0], [1.5, 2.2,0], color = GREEN)
        lcq2 = Line([4,2.2, 0], [4.2, 2.2,0], color = GREEN)
        claimb = Line([-3.15,2.2, 0], [-4.25, 2.2,0])
        self.play(Write(cla), Create(lcp1), Create(lcp2), Create(lcq1), Create(lcq2), Create(claimb))
        self.wait()
        c11 = Ellipse( width=0.75, height=1.6, color=RED ).shift(0.53*DOWN)
        g = ParametricFunction(lambda t : [x(t) ,y(t) ,0], t_range = [-PI,PI], color = BLUE)
        g1 = ParametricFunction(lambda t : [x(t) ,y(t) ,0], t_range = [-PI/8,7*PI/8], color = ORANGE)
        g2 = ParametricFunction(lambda t : [x2(t) ,y2(t) ,0], t_range = [-PI,PI], color = BLUE)
        g12 = ParametricFunction(lambda t : [x2(t) ,y2(t) ,0], t_range = [-PI/8,7*PI/8], color = ORANGE)
        p1 = Dot([x(-PI/8), y(-PI/8), 0], color = YELLOW)
        q1 = Dot([x(7*PI/8), y(7*PI/8), 0], color = GREEN)
        p2 = Dot([x(3*PI/8), y(3*PI/8), 0], color = YELLOW)
        c12 = Ellipse( width=1.6, height=0.75, color=TEAL ).shift(0.53*DOWN)
        q2 = Dot([x2(-5*PI/8), y2(-5*PI/8), 0], color = GREEN)
        p1l = Tex(r"$p_1$", font_size = 30).shift( (1.35)*DOWN)
        p2l = Tex(r"$p_2$", font_size = 30).next_to(p2, RIGHT)
        gam1l = Tex(r"$\gamma_1$", font_size = 30).shift(2*RIGHT + 2*DOWN)
        gam1b = Line([gam1l.get_center()[0]-0.1, gam1l.get_center()[1]-0.2, 0], [gam1l.get_center()[0]+0.1, gam1l.get_center()[1]-0.2, 0], color = ORANGE)
        p2lb = Line([0,0,0], [0.2, 0,0], color = YELLOW).next_to(p2l, DOWN).shift(0.2*UP)
        p1lb = Line([-0.1,-1.5, 0], [0.1, -1.5, 0], color = YELLOW)
        q1l = Tex(r"$q_1$", font_size = 30).shift((0.45)* UP)
        q1lb = Line([-0.1,0.3, 0], [0.1, 0.3, 0], color = GREEN)
        q2l = Tex(r"$q_2$", font_size = 30).next_to(q2, LEFT)
        q2lb = Line([0,0, 0], [0.2, 0, 0], color = GREEN).next_to(q2l, DOWN).shift(0.2*UP)
        c1 = Circle(radius = sqrt( (x(-PI/8) - x(7*PI/8) )**2 + (y(-PI/8) - y(7*PI/8) )**2 )/2).shift((x(-PI/8) + x(7*PI/8))*RIGHT/2 + (y(-PI/8) + y(7*PI/8))*UP/2 )
        self.play(Create(g))
        self.wait()
        self.play(Create(p1), Create(q1), Create(c1), Write(p1l), Write(q1l), Create(p1lb), Create(q1lb), Create(g1), Create(s1lab), Create(s1b), Create(gam1l), Create(gam1b))
        self.wait()
        self.play(Create(p2), Create(p2l), Create(p2lb))
        self.wait()
        self.play(Transform(g,g2),Transform(g1, g12) , Transform(c1, c11), p1.animate.shift(0.5*DOWN), p1l.animate.shift(0.4*DOWN), p1lb.animate.shift(0.4*DOWN), q1.animate.shift(0.4*UP), q1l.animate.shift(0.4*UP), q1lb.animate.shift(0.4*UP), gam1l.animate.shift(0.7*RIGHT), gam1b.animate.shift(0.7*RIGHT))
        self.wait()
        self.play(Create(c12), Create(q2), Create(q2lb), Write(q2l), Create(s2lab), Create(s2b))
        self.wait()







class fvt(Scene):
    def construct(self):
        defi = Tex(r"Definition", font_size=40, color= BLUE).to_edge(UL).shift(0.5*DOWN)
        l = defi.get_left()[0]
        def1 = Tex(r"For a smooth regular curve $\gamma : [a,b] \to \mathbb{R}^2$, the point $\gamma (t)$ is a ", r"vertex", r" if", font_size = 40).next_to(defi, DOWN)
        def1.shift((def1.get_left()[0]-l)*LEFT)
        def2 = Tex(r"$k ^{\prime} (t) = 0$", font_size = 40).next_to(def1, DOWN)
        def2.shift((def2.get_center()[0])*LEFT)
        ver1 = Line([0,0,0], [1, 0,0], color = PINK).next_to(def1[1], DOWN).shift(0.2*UP)
        vb = SurroundingRectangle(def2, color = PINK, buff = 0.1)
        self.play(Write(defi), Write(def1), Write(def2), Create(ver1), Create(vb))
        self.wait()
        thm = Tex(r"Theorem", font_size=40, color= BLUE).next_to(def2, DOWN)
        thm.shift((thm.get_left()[0]-l)*LEFT)
        t1 = Tex(r"A closed smooth regular curve $\gamma : \mathbb{S}^1 \to \mathbb{R}^2$ has at least four vertices", font_size = 40).next_to(thm, DOWN)
        t1.shift((t1.get_left()[0]-l)*LEFT)
        self.play(Write(thm), Write(t1))
        self.wait()
        g = ParametricFunction(lambda t : [ 2*np.cos(t) , np.sin(t) -1.5 ,0], t_range = [0,TAU], color = BLUE)
        def r(t):
            return sqrt(1 + 3* np.sin(t)**2)**3/2
        def n10(t):
            return -2*np.cos(t) + 6*np.cos(t)*(np.sin(t)**2)/(1 + 3* np.sin(t)**2)
        def n20(t):
            return -np.sin(t) -  3*np.sin(t)*(np.cos(t)**2)/(1 + 3* np.sin(t)**2)
        def n1(t):
            return n10(t)/sqrt( n10(t)**2 + n20(t)**2 )
        def n2(t):
            return n20(t)/sqrt( n10(t)**2 + n20(t)**2 )
        def c1(t):
            return 2*np.cos(t) + n1(t)*r(t)
        def c2(t):
            return np.sin(t) + n2(t)*r(t) -1.5
        oc = Circle(radius = r(PI/8)).shift( c1(PI/8)*RIGHT + c2(PI/8)*UP)
        p = Dot([2*np.cos(PI/8), np.sin(PI/8) -1.5, 0 ], color = YELLOW)
        self.play(Create(g), Create(p))
        self.wait()
        self.play(Create(oc))
        self.wait()
        s = ValueTracker(PI/8)
        oc.add_updater( lambda x : x.become(Circle(radius = r(s.get_value())).shift( c1(s.get_value())*RIGHT + c2(s.get_value())*UP))) 
        p.add_updater( lambda x : x.become(Dot( [2*np.cos(s.get_value()), np.sin(s.get_value()) - 1.5, 0 ], color = YELLOW ))) 
        self.play(s.animate.set_value(PI/2), run_time =1)
        v1 = Dot([0,-0.5,0], color = PINK)
        v2 = Dot([-2,-1.5,0], color = PINK)
        v3 = Dot([0,-2.5,0], color = PINK)
        v4 = Dot([2,-1.5,0], color = PINK)
        self.play(Create(v1))
        self.play(s.animate.set_value(PI), run_time =2)
        self.play(Create(v2))
        self.play(s.animate.set_value(3*PI/2), run_time =2)
        self.play(Create(v3))
        self.play(s.animate.set_value(2*PI), run_time =2)
        self.play(Create(v4))
        self.wait()
        




class fp(Scene):
    def construct(self):
        def x(t):
            return (7/4 + np.cos(2*t))*np.cos(t)
        def y(t):
            return (7/4 + np.cos(2*t))*np.sin(t)*2/3
        def u(t):
            return  x(t)**2 + y(t)**2
        g = ParametricFunction(lambda t : [ x(t) , y(t) ,0], t_range = [0,TAU], color = BLUE)
        g0 = ParametricFunction(lambda t : [ x(t) , y(t) ,0], t_range = [0,TAU], color = BLUE)
        c = Circle(radius = sqrt(1.3), color = YELLOW)
        ginv =  ParametricFunction(lambda t : [ (1.3)*x(t) / u(t) , (1.3)*y(t) /u(t) ,0], t_range = [0,TAU], color = BLUE)
        self.play(Create(g))
        self.wait()
        oc1 = Circle( radius =  11**2 / 243  ).shift((11/4 - 11**2 / 243  )*RIGHT )
        oc10 = Circle( radius =  11**2 / 243  ).shift((11/4 - 11**2 / 243  )*RIGHT )
        oc1inv = Circle( radius =  10.4/77  ).shift((46.8/77 )*RIGHT )
        oc3 = Circle( radius =  11**2 / 243  ).shift((11/4 - 11**2 / 243  )*LEFT )
        oc30 = Circle( radius =  11**2 / 243  ).shift((11/4 - 11**2 / 243  )*LEFT )
        oc3inv = Circle( radius =  10.4/77  ).shift((46.8/77 )*LEFT )
        oc2 = Circle(  radius = 27/104 ).shift( ( 1/2 + 27/104 )*UP  )
        oc2inv = Circle(  radius = (1.3)*27/53 ).shift((1.3)*(79/53)*UP  )
        oc4 = Circle(  radius = 27/104 ).shift( ( 1/2 + 27/104 )*DOWN  )
        oc4inv = Circle(  radius = (1.3)*27/53 ).shift((1.3)*(79/53)*DOWN  )
        v1 = Dot([11/4,0,0], color = PINK)
        v10 = Dot([11/4,0,0], color = PINK)
        v1l = Tex(r"$v_1$", font_size=40).next_to(v1, RIGHT)
        v1b = Line([0,0,0], [0.2, 0,0], color = PINK).next_to(v1l, DOWN).shift(0.2*UP)
        v1inv = Dot( [ 5.2/11, 0, 0], color = PINK )
        v2 = Dot([0,1/2,0], color = PINK)
        v2inv = Dot([ 0, 2.6 ,0  ], color=PINK)
        v3 = Dot([-11/4,0,0], color = PINK)
        v30 = Dot([-11/4,0,0], color = PINK)
        v3l = Tex(r"$v_2$", font_size=40).next_to(v3, LEFT)
        v3b = Line([0,0,0], [0.2, 0,0], color = PINK).next_to(v3l, DOWN).shift(0.2*UP)
        v3inv = Dot( [ -5.2/11, 0, 0], color = PINK )
        v4 = Dot([0,-1/2,0], color = PINK)
        v4inv = Dot([ 0,-2.6 ,0  ], color=PINK)
        v2l = Tex(r"$v_3$", font_size=40).next_to(v2, DOWN).shift(0.3*RIGHT+ 0.3*UP)
        v2b = Line([0,0,0], [0.2, 0,0], color = PINK).next_to(v2l, DOWN).shift(0.2*UP)
        v4l = Tex(r"$v_4$", font_size=40).next_to(v4, UP).shift(0.3*LEFT + 0.2*DOWN)
        v4b = Line([0,0,0], [0.2, 0,0], color = PINK).next_to(v4l, DOWN).shift(0.2*UP)
        self.play(Create(oc1), Create(v1), Create(v1b), Write(v1l))
        self.wait()
        self.play(Create(oc3), Create(v3), Create(v3b), Write(v3l))
        self.wait()
        self.play(Create(c))
        self.wait()
        self.play(Transform(g,ginv), Transform(v1, v1inv), Transform(v3, v3inv), Transform(oc1, oc1inv), Transform(oc3, oc3inv), FadeOut(v1l), FadeOut(v1b), FadeOut(v3l), FadeOut(v3b))
        self.wait()
        self.play(Create(oc2inv), Create(v2inv), Create(v4inv), Create(oc4inv))
        self.wait()
        self.play(Transform(g, g0), Transform(v1, v10), Transform(v3, v30), Transform(oc1, oc10), Transform(oc3, oc30), Write(v1l), Create(v1b), Write(v3l), Create(v3b), Write(v4l), Create(v4b), Write(v2l), Create(v2b), Transform(v2inv, v2), Transform(v4inv, v4), Transform(oc2inv, oc2), Transform(oc4inv, oc4))
        self.wait()
        self.play(FadeOut(c))
        self.wait()
        ch = Tex(r"Check:" ,r" A curve can support its osculating circle only at a vertex", font_size = 36).shift(2.5*UP)
        ch2 = Tex(r"and inversion sends osculating circles to osculating circles.", font_size = 36).next_to(ch,DOWN)
        ul = Line([0,0,0], [1,0,0]).next_to(ch[0], DOWN).shift(0.1*UP)
        self.play(Write(ch), Write(ch2), Create(ul))
        self.wait()
        

        