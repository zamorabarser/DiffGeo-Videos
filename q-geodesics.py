from asyncio import threads
from this import d
from tkinter import E
from manim import *
from numpy import sqrt
import math



class ti(Scene):
    def construct(self):
        t1 = Text("Geodesics", font_size=60)
        self.play(Write(t1))
        self.wait()        


class dpq(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 45*DEGREES)
        sl = Tex(r"Let ", r"$p$,",  r" $q$", r" $ \in$ ", r"$\Sigma$", font_size = 34).shift(2.6*UP + 4*LEFT)
        qli = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(sl[1], DOWN).shift(0.15*UP+0.03*RIGHT)
        pli = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(sl[2], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(sl[4], DOWN).shift(0.07*UP)
        dpq = Tex(r"$d(p,q) : = \inf \{  $ length$(\gamma) \vert \gamma : [a,b] \to \Sigma , \gamma (a) = p, \gamma (b) = q \} $", font_size = 34).next_to(sl, DOWN).shift(0.5*DOWN + 4*RIGHT)
        box = SurroundingRectangle(dpq, buff = .1, color = BLUE)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 5, y_length = 5, z_length = 4)
        def f(u,v):
            return (u**2 - v**2) / 8 - 2
        s = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v) ),
            u_range = [-4,4],
            v_range = [-4,4],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        p = Sphere(radius = 0.07)
        q = p.copy()
        p.set_color(ORANGE).shift(ax.c2p(0, 2, f(0, 2) ))
        q.set_color(ORANGE).shift(ax.c2p(0, -2, f(0,-2)  ))
        g1 = ParametricFunction(lambda t : ax.c2p( (t+2)*(2-t)/2, t,f((t+2)*(2-t)/2, t) ) , t_range = [-2,2], color = BLUE)
        g2 = ParametricFunction(lambda t : ax.c2p( 1.7*np.sin(PI*t/2) , t ,f(1.7*np.sin(PI*t/2), t) ) , t_range = [-2,2], color = BLUE)
        g3 = ParametricFunction(lambda t : ax.c2p( 0, t,f(0, t) ) , t_range = [-2,2], color = BLUE)
        self.add_fixed_in_frame_mobjects(sl, qli, pli, dpq, sli, box)
        self.remove(sl, qli, pli, dpq, sli, box)
        self.play(Create(s), Create(p), Create(q), Write(sl), Create(qli), Create(pli), Create(sli))
        self.wait()
        self.play(Create(dpq), Create(box))
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/8)
        self.wait()
        self.play(Create(g1))
        self.wait()
        self.play(FadeOut(g1))
        self.wait()
        self.play(Create(g2))
        self.wait()
        self.play(FadeOut(g2))
        self.wait()
        self.play(Create(g3))
        self.wait()
        self.stop_ambient_camera_rotation()
        self.wait()
        #self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(sl), FadeOut(n2), FadeOut(sli), FadeOut(pli))



class minim(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 13*PI/8)
        sl = Tex(r"Let ", r"$p$,",  r" $q$", r" $ \in$ ", r"$\Sigma$", font_size = 34).shift(2.6*UP + 4*LEFT)
        qli = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(sl[1], DOWN).shift(0.15*UP+0.03*RIGHT)
        pli = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(sl[2], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(sl[4], DOWN).shift(0.07*UP)
        dpq = Tex(r"$d(p,q) : = \inf \{  $ length$(\gamma) \vert \gamma : [a,b] \to \Sigma , \gamma (a) = p, \gamma (b) = q \} $", font_size = 34).next_to(sl, DOWN).shift(0.5*DOWN + 4*RIGHT)
        box = SurroundingRectangle(dpq, buff = .1, color = BLUE)
        mindef = Tex(r"$\gamma : [a,b] \to \Sigma$ is called ", r"minimizing", r" if it attains the infimum.", font_size = 34).next_to(dpq, DOWN).shift(UP)
        minimli = Line([0,0,0], [1.7,0,0], color = BLUE).next_to(mindef[1], DOWN).shift(0.15*UP)
        minal = Tex(r"For all $\alpha : [a,b] \to \Sigma$ with $\alpha (a) = p , \alpha (b) = q$, one has", font_size = 34).next_to(mindef, DOWN)#.shift(0.1*UP)
        minal2 = Tex(r"length$(\alpha) \geq $ length$(\gamma)$.", font_size = 34).next_to(minal, DOWN)#.shift(0.1*UP)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 5, y_length = 5, z_length = 4)
        def f(u,v):
            return (u**2 - v**2) / 8 - 2
        s = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v) ),
            u_range = [-4,4],
            v_range = [-4,4],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        p = Sphere(radius = 0.07)
        q = p.copy()
        p.set_color(ORANGE).shift(ax.c2p(0, 2, f(0, 2) ))
        q.set_color(ORANGE).shift(ax.c2p(0, -2, f(0,-2)  ))
        g3 = ParametricFunction(lambda t : ax.c2p( 0, t,f(0, t) ) , t_range = [-2,2], color = BLUE)
        self.add_fixed_in_frame_mobjects(sl, qli, pli, dpq, sli, box, mindef, minimli, minal, minal2)
        self.remove(mindef, minimli, minal, minal2)
        self.add( s, p, q, g3 )
        self.wait()
        self.play(FadeOut(sl), FadeOut(qli), FadeOut(pli), FadeOut(sli), box.animate.shift(UP), dpq.animate.shift(UP), s.animate.shift([0,0,-0.5]), p.animate.shift([0,0,-0.5]), q.animate.shift([0,0,-0.5]), g3.animate.shift([0,0,-0.5]))
        self.wait()
        self.play(Write(mindef), Create(minimli))
        self.wait()
        self.play(Write(minal), Create(minal2))
        self.wait()
        #self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(sl), FadeOut(n2), FadeOut(sli), FadeOut(pli))
        



class ugex(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= - PI/4)
        sl = Tex(r"$\Sigma$", r" $ = \mathbb{R}^2 \backslash \{ 0 \} $, ", r"$p$", r" $ = (1,0)$, ", r"$q$", r" $ = (-1,0)$", font_size = 34).shift(2.1*UP)
        qli = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(sl[4], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(sl[2], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(sl[0], DOWN).shift(0.07*UP)
        uge = Tex(r"$\Rightarrow$ There is no minimizing curve between $p$ and $q$.", font_size = 34).next_to(sl, DOWN).shift(0.5*DOWN)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 10, y_length = 10, z_length = 8)
        h = -1
        s = Surface(
            lambda u, v: ax.c2p( u , v , h ),
            u_range = [-3,3],
            v_range = [-3,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        p = Sphere(radius = 0.07)
        q = p.copy()
        p.set_color(ORANGE).shift(ax.c2p(2, 0, h ))
        q.set_color(ORANGE).shift(ax.c2p(-2, 0, h  ))
        z = Dot( ax.c2p(0,0,h) ,color = BLACK, radius = 0.1)
        g1 = ParametricFunction(lambda t : ax.c2p( t, np.cos(PI*t/4), h ) , t_range = [-2,2], color = BLUE)
        g2 = ParametricFunction(lambda t : ax.c2p( t, np.cos(PI*t/4)/2, h ) , t_range = [-2,2], color = BLUE)
        g3 = ParametricFunction(lambda t : ax.c2p( t, np.cos(PI*t/4)/4, h ) , t_range = [-2,2], color = BLUE)
        g4 = ParametricFunction(lambda t : ax.c2p( t, 0, h ) , t_range = [-2,2], color = RED)
        self.add_fixed_in_frame_mobjects(sl, qli, pli, uge, sli)
        self.remove(sl, qli, pli, sli, uge)
        self.play(Create(sl), Create(pli), Create(qli), Create(sli), Create(uge), Create(s), Create(p), Create(q), Create(z))
        self.wait()
        self.play(Create(g1))
        self.wait()
        self.play(FadeOut(g1))
        self.wait()
        self.play(Create(g2))
        self.wait()
        self.play(FadeOut(g2))
        self.wait()
        self.play(Create(g3))
        self.wait()
        self.play(FadeOut(g3))
        self.wait()
        self.play(Create(g4))
        self.wait()
        self.play(FadeOut(g4))
        self.wait()
        #self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(sl), FadeOut(n2), FadeOut(sli), FadeOut(pli))
        



class spex(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= - 3*PI/16)
        sl = Tex(r"$\Sigma$", r" $ = \mathbb{S}^2 $, ", r"$p$", r" $ = (1,0)$, ", r"$q$", r" $ = (-1,0)$", font_size = 34).shift(2.4*UP)
        qli = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(sl[4], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(sl[2], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(sl[0], DOWN).shift(0.07*UP)
        uge = Tex(r"$\Rightarrow$ There are infinitely many minimizing curves between $p$ and $q$.", font_size = 34).next_to(sl, DOWN).shift(0.3*DOWN)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 15, y_length = 15, z_length = 12)
        h = - 0.8
        s = Surface(
            lambda u, v: ax.c2p( np.cos(u)*np.sin(v) , np.sin(u)*np.sin(v) , np.cos(v)  + h),
            u_range = [-PI, PI],
            v_range = [0, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        p = Sphere(radius = 0.07)
        q = p.copy()
        p.set_color(ORANGE).shift(ax.c2p(0, 0, 1 + h ))
        q.set_color(ORANGE).shift(ax.c2p(0, 0, -1 + h ))
        g1 = ParametricFunction(lambda t : ax.c2p(  np.sin(t) , 0 , h + np.cos(t)  ) , t_range = [0,PI], color = BLUE)
        g2 = ParametricFunction(lambda t : ax.c2p(  np.sin(t) /sqrt(2) , np.sin(t) /sqrt(2)  , h + np.cos(t)  ) , t_range = [0,PI], color = BLUE)
        g3 = ParametricFunction(lambda t : ax.c2p(  - np.sin(t) , 0 , h + np.cos(t)  ) , t_range = [0,PI], color = BLUE)
        g4 = ParametricFunction(lambda t : ax.c2p(  - np.sin(t) /sqrt(2) , np.sin(t) /sqrt(2)  , h + np.cos(t)  ) , t_range = [0,PI], color = BLUE)
        g5 = ParametricFunction(lambda t : ax.c2p(  0 ,  np.sin(t)  , h + np.cos(t)  ) , t_range = [0,PI], color = BLUE)
        g6 = ParametricFunction(lambda t : ax.c2p(  np.sin(t) /sqrt(2) , - np.sin(t) /sqrt(2)  , h + np.cos(t)  ) , t_range = [0,PI], color = BLUE)
        g7 = ParametricFunction(lambda t : ax.c2p(  0 , - np.sin(t) , h + np.cos(t)  ) , t_range = [0,PI], color = BLUE)
        g8 = ParametricFunction(lambda t : ax.c2p(  - np.sin(t) /sqrt(2) , - np.sin(t) /sqrt(2)  , h + np.cos(t)  ) , t_range = [0,PI], color = BLUE)
        self.add_fixed_in_frame_mobjects(sl, qli, pli, uge, sli)
        self.remove(sl, qli, pli, sli, uge)
        self.play(Create(sl), Create(pli), Create(qli), Create(sli), Create(uge), Create(s), Create(p), Create(q))
        self.wait()
        self.play(Create(g1), run_time = 0.5)
        self.play(Create(g6), run_time = 0.5)
        self.play(Create(g7), run_time = 0.5)
        self.play(Create(g8), run_time = 0.5)
        self.play(Create(g3), run_time = 0.5)
        self.play(Create(g4), run_time = 0.5)
        self.play(Create(g5), run_time = 0.5)
        self.play(Create(g2), run_time = 0.5)
        self.wait()
        #self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(sl), FadeOut(n2), FadeOut(sli), FadeOut(pli))
        



class hr1(Scene):
    def construct(self):        
        pr = Tex(r"Proposition", font_size=34, color = BLUE).to_edge(UL).shift(1.4*DOWN + 0.5*RIGHT)
        l = pr.get_left()[0]
        pr1 = Tex(r"If a surface $\Sigma$ is a closed subset of $\mathbb{R}^3$, then $(\Sigma , d) $ is a complete metric space.", font_size=34).next_to(pr,DOWN)
        hi = Tex(r"Hint: ", r"A Cauchy sequence on $(\Sigma , d)$  is also a Cauchy sequence on $(\mathbb{R}^3, \vert \cdot \vert )$.", font_size=34).next_to(pr1, DOWN)
        hi[0].set_color(BLUE)
        th = Tex(r"Theorem (Hopf-Rinow I)", font_size=34, color = BLUE).next_to(hi, DOWN).shift(0.7*DOWN)
        th1 = Tex(r"If  $(\Sigma , d)$ is a complete metric space, then for all $p,q \in \Sigma$,", font_size=34).next_to(th,DOWN)
        th2 = Tex(r"there is a minimizing curve between $p$ and $q$.", font_size=34).next_to(th1,DOWN)
        pr1.shift((l - pr1.get_left()[0])*RIGHT)
        hi.shift((l - hi.get_left()[0])*RIGHT)
        th.shift((l - th.get_left()[0])*RIGHT)
        th1.shift((l - th1.get_left()[0])*RIGHT)    
        th2.shift((l - th2.get_left()[0])*RIGHT)
        self.play(Write(pr), Write(pr1))
        self.wait()
        self.play(Write(hi))
        self.wait()
        self.play(Write(th), Write(th1), Write(th2))
        self.wait()
        



class ener(Scene):
    def construct(self):        
        ende = Tex(r"Definition", font_size=34, color = BLUE).to_edge(UL).shift(0.5*DOWN + 0.2*RIGHT)
        l = ende.get_left()[0]
        ende1 = Tex(r"For a curve $\gamma : [a,b] \to \Sigma$, we define its " , r"energy", r" as", font_size=34).next_to(ende,DOWN)
        ende1[1].set_color(BLUE)
        ende2 = Tex(r"$E (\gamma ) : =  \frac{1}{2}  \int_a^b \vert \gamma ^{\prime} (t) \vert ^2 dt$.", font_size=34).next_to(ende1, DOWN)
        enex = Tex(r"Exercise", font_size=34, color = BLUE).next_to(ende2, DOWN)#.shift(0.3*DOWN)
        enex1 = Tex(r"A smooth curve $\gamma : [a,b] \to \Sigma$ minimizes the energy among all amooth curves", font_size=34).next_to(enex,DOWN)
        enex2 = Tex(r"$\alpha : [a,b] \to \Sigma$ with the same endpoints if and only if it has constant speed", font_size=34).next_to(enex1,DOWN)
        enex3 = Tex(r"and it minimizes the length among smooth curves with the same endpoints.", font_size=34).next_to(enex2,DOWN)
        hi = Tex(r"Hint: ",r"For each continuous function $f : [a,b] \to \mathbb{R}$ one has", font_size=34).next_to(enex3, DOWN)
        hi[0].set_color(BLUE)
        hi1 = Tex(r"$\int_a^b f^2 \geq \frac{ \left( \int_a^b f \right) ^2 }{b-a}$,", font_size=34).next_to(hi,DOWN)
        hi2 = Tex(r"and equality holds if and only if $f$ is constant.", font_size=34).next_to(hi1,DOWN)
        ende1.shift((l - ende1.get_left()[0])*RIGHT)
        ende2.shift(ende2.get_center()[0]*LEFT)
        hi1.shift((l - hi1.get_left()[0])*RIGHT)
        enex.shift((l - enex.get_left()[0])*RIGHT)
        enex1.shift((l - enex1.get_left()[0])*RIGHT)    
        enex2.shift((l - enex2.get_left()[0])*RIGHT)  
        enex3.shift((l - enex3.get_left()[0])*RIGHT)
        hi.shift((l - hi.get_left()[0])*RIGHT)
        hi2.shift((l - hi2.get_left()[0])*RIGHT)    
        hi1.shift(hi1.get_center()[0]*LEFT)
        self.play(Write(ende), Write(ende1), Write(ende2))
        self.wait()
        self.play(Write(enex), Write(enex1), Write(enex2), Write(enex3))
        self.wait()
        self.play(Write(hi), Write(hi1), Write(hi2))
        self.wait()
        


class varde(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta= - PI/3)
        thel = Tex(r"$\theta : (-\varepsilon, \varepsilon ) \times [a,b] \to \Sigma $ is a ", r"variation", r" of $\gamma : [a,b] \to \Sigma$", font_size = 34).shift(2.4*UP)
        thel[1].set_color(BLUE)
        thel2 = Tex(r"if $\theta (0,t) = \gamma (t) $ for all $t \in [a,b]$.", font_size = 34).next_to(thel, DOWN)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,2,10], x_length = 5, y_length = 5, z_length = 2, tips = False)
        axd = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 5, y_length = 5, z_length = 0.01, tips = False)
        ax.shift([3, sqrt(3), -1])
        axd.shift([-3, -sqrt(3), -1.7])
        def f(u,v):
            return (np.sin(u) + np.cos(v))/1.2
        s = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v)),
            u_range = [-4,4],
            v_range = [-4,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        dom = Surface(
            lambda u, v: axd.c2p( v , u , 0 ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7, 
            resolution=(3,6)
        )
        thes = Surface(
            lambda u, v: ax.c2p( v , u -1  , f(v,u -1 ) ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        thar = CurvedArrow( [-1.4 ,-0.5 ,0],[0.2 ,-0.2,0], radius = -3 , tip_length = 0.2)
        tharl = Tex(r"$\theta$", font_size = 34).next_to(thar, UP) 
        g0 = ParametricFunction(lambda t : axd.c2p( t ,  0 , 0  ) , t_range = [-3,3], color = BLUE)
        g1 = ParametricFunction(lambda t : ax.c2p(  t, -1 , f(t,1)  ) , t_range = [-3,3], color = BLUE)
        self.add_fixed_in_frame_mobjects(thel, thel2, thar, tharl)
        self.remove(thel, thel2, thar, tharl)
        self.wait()
        self.play(Create(thel), Create(thel2), Create(ax), Create(axd), Create(s), Create(dom), Create(g1), Create(thar), Create(tharl), Create(g0))
        self.wait()
        self.play(Transform(dom, thes))
        self.play(FadeOut(s))
        self.wait()
        #self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(sl), FadeOut(n2), FadeOut(sli), FadeOut(pli))
        


class varde(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta= - PI/3)
        thel = Tex(r"$\theta : (-\varepsilon, \varepsilon ) \times [a,b] \to \Sigma $ is a ", r"variation", r" of $\gamma : [a,b] \to \Sigma$ if", font_size = 34).shift(2.4*UP)
        thel[1].set_color(BLUE)
        thel2 = Tex(r"$\theta (0,t) = \gamma (t) $ for all $t \in [a,b]$.", font_size = 34).next_to(thel, DOWN)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,2,10], x_length = 5, y_length = 5, z_length = 2, tips = False)
        axd = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 5, y_length = 5, z_length = 0.01, tips = False)
        ax.shift([3, sqrt(3), -1])
        axd.shift([-3, -sqrt(3), -1.7])
        def f(u,v):
            return (np.sin(u) + np.cos(v))/1.2
        s = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v)),
            u_range = [-4,4],
            v_range = [-4,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        str = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v)),
            u_range = [-4,4],
            v_range = [-4,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.2,
            stroke_width = 0.1
        )
        dom = Surface(
            lambda u, v: axd.c2p( v , u , 0 ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7, 
            resolution=(3,6)
        )
        dom1 = dom.copy()
        thes = Surface(
            lambda u, v: ax.c2p( v , u -1  , f(v,u -1 ) ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7,
            resolution=(3,24)
        )
        thar = CurvedArrow( [-1.4 ,-0.5 ,0],[0.2 ,-0.2,0], radius = -3 , tip_length = 0.2)
        tharl = Tex(r"$\theta$", font_size = 34).next_to(thar, UP).shift(0.2*DOWN)
        g0 = ParametricFunction(lambda t : axd.c2p( t ,  0 , 0  ) , t_range = [-3,3], color = BLUE)
        g1 = ParametricFunction(lambda t : ax.c2p(  t, -1 , f(t,1)  ) , t_range = [-3,3], color = BLUE)
        self.add_fixed_in_frame_mobjects(thel, thel2, thar, tharl)
        self.remove(thel, thel2, thar, tharl)
        self.wait()
        self.play(Create(thel), Create(thel2), Create(ax), Create(axd), Create(s), Create(dom), Create(g1), Create(thar), Create(tharl), Create(g0))
        self.wait()
        self.play(Transform(dom1, thes), Transform(s,str))
        self.wait()
        #self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(sl), FadeOut(n2), FadeOut(sli), FadeOut(pli))
        


class vdef(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta= - PI/3)
        thel = Tex(r"$\theta : (-\varepsilon, \varepsilon ) \times [a,b] \to \Sigma $ is a ", r"variation", r" of $\gamma : [a,b] \to \Sigma$ if", font_size = 34).shift(2.4*UP)
        thel[1].set_color(BLUE)
        thel2 = Tex(r"$\theta (0,t) = \gamma (t) $ for all $t \in [a,b]$.", font_size = 34).next_to(thel, DOWN)
        vla = Tex(r"$V ( t)$", r" $   : = \theta _ s (0,t) $ is called the ", r"direction", r" of $\theta$.", font_size = 34).next_to(thel2, DOWN)
        vla[2].set_color(BLUE)
        vli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(vla[0], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,2,10], x_length = 5, y_length = 5, z_length = 2, tips = False)
        axd = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 5, y_length = 5, z_length = 0.01, tips = False)
        ax.shift([3, sqrt(3), -1])
        axd.shift([-3, -sqrt(3), -1.7])
        def f(u,v):
            return (np.sin(u) + np.cos(v))/1.2
        dom = Surface(
            lambda u, v: axd.c2p( v , u , 0 ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7, 
            resolution=(3,6)
        )
        thes = Surface(
            lambda u, v: ax.c2p( v , u -1  , f(v,u -1 ) ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7,
            resolution=(3,24)
        )
        str = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v)),
            u_range = [-4,4],
            v_range = [-4,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.2,
            stroke_width = 0.1
        )
        v = Arrow(ax.c2p(-3,-0.88 ,f(-3,-1) + np.sin(1)/10 ), ax.c2p( -3, -2.92 , f(-3,-1) - 1.6*np.sin(1)  ) , color = YELLOW, max_tip_length_to_length_ratio=0.2)
        vd = Arrow(axd.c2p( -3 , 0.12 , 0 ), axd.c2p( -3, -2.8 , 0  ) , color = YELLOW, max_tip_length_to_length_ratio=0.2)        
        z = ValueTracker(-3)
        v.add_updater(lambda x : x.become(Arrow(ax.c2p(z.get_value(),-0.88 ,f(z.get_value(),-1) + np.sin(1)/10 ), ax.c2p( z.get_value(), -2.92 , f(z.get_value(),-1) - 1.6* np.sin(1)  ),  color = YELLOW , max_tip_length_to_length_ratio=0.2)))
        vd.add_updater(lambda x : x.become(Arrow(axd.c2p(z.get_value(), 0.12 , 0 ), axd.c2p( z.get_value(), -2.8 ,  0  ),  color = YELLOW , max_tip_length_to_length_ratio=0.2 )))
        thar = CurvedArrow( [-1.4 ,-0.5 ,0],[0.2 ,-0.2,0], radius = -3 , tip_length = 0.2)
        tharl = Tex(r"$\theta$", font_size = 34).next_to(thar, UP).shift(0.2*DOWN)
        g0 = ParametricFunction(lambda t : axd.c2p( t ,  0 , 0  ) , t_range = [-3,3], color = BLUE)
        g1 = ParametricFunction(lambda t : ax.c2p(  t, -1 , f(t,-1)  ) , t_range = [-3,3], color = BLUE)
        self.add_fixed_in_frame_mobjects(thel, thel2, thar, tharl, vla, vli)
        self.add(g0,g1, ax, axd, dom, thes, str)
        self.remove(vla, vli)
        self.wait()
        self.play(Create(vla), Create(v), Create(vli), Create(vd))
        self.wait()
        self.play( z.animate.set_value(3) , run_time = 3 )
        self.wait()



class vvf(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta= - PI/3)
        thel = Tex(r"$\theta : (-\varepsilon, \varepsilon ) \times [a,b] \to \Sigma $ is a ", r"variation", r" of $\gamma : [a,b] \to \Sigma$ if", font_size = 34).shift(2.4*UP)
        thel[1].set_color(BLUE)
        thel2 = Tex(r"$\theta (0,t) = \gamma (t) $ for all $t \in [a,b]$.", font_size = 34).next_to(thel, DOWN)
        vla = Tex(r"$V ( t)$", r" $   : = \theta _ s (0,t) $ is called the ", r"direction", r" of $\theta$.", font_size = 34).next_to(thel2, DOWN)
        vla[2].set_color(BLUE)
        vli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(vla[0], DOWN).shift(0.15*UP)
        vvf = Tex(r"$V(t) \in T_{\gamma (t)} \Sigma$ for all $t \in [a,b]$", font_size = 34).next_to(vla, RIGHT).shift(1.35*LEFT + 0.4*DOWN)
        box = SurroundingRectangle(vvf, buff = .1, color = YELLOW)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,2,10], x_length = 5, y_length = 5, z_length = 2, tips = False)
        axd = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 5, y_length = 5, z_length = 0.01, tips = False)
        ax.shift([3, sqrt(3), -1])
        axd.shift([-3, -sqrt(3), -1.7])
        def f(u,v):
            return (np.sin(u) + np.cos(v))/1.2
        dom = Surface(
            lambda u, v: axd.c2p( v , u , 0 ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7, 
            resolution=(3,6)
        )
        thes = Surface(
            lambda u, v: ax.c2p( v , u -1  , f(v,u -1 ) ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7,
            resolution=(3,24)
        )
        str = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v)),
            u_range = [-4,4],
            v_range = [-4,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.2,
            stroke_width = 0.1
        )
        v = Arrow(ax.c2p(3,-0.88 ,f(3,-1) + np.sin(1)/10 ), ax.c2p( 3, -2.92 , f(3,-1) - 1.6*np.sin(1)  ) , color = YELLOW, max_tip_length_to_length_ratio=0.2)
        vd = Arrow(axd.c2p( 3 , 0.12 , 0 ), axd.c2p( 3, -2.8 , 0  ) , color = YELLOW, max_tip_length_to_length_ratio=0.2)  
        thar = CurvedArrow( [-1.4 ,-0.5 ,0],[0.2 ,-0.2,0], radius = -3 , tip_length = 0.2)
        tharl = Tex(r"$\theta$", font_size = 34).next_to(thar, UP).shift(0.2*DOWN)
        g0 = ParametricFunction(lambda t : axd.c2p( t ,  0 , 0  ) , t_range = [-3,3], color = BLUE)
        g1 = ParametricFunction(lambda t : ax.c2p(  t, -1 , f(t,-1)  ) , t_range = [-3,3], color = BLUE)
        self.add_fixed_in_frame_mobjects(thel, thel2, thar, tharl, vla, vli, vvf, box)
        self.add(g0,g1, ax, axd, dom, thes, str, v, vd)
        self.remove(vvf, box)
        self.wait()
        self.play(Create(vvf), Create(box), thel.animate.shift(2.2*LEFT), thel2.animate.shift(2.2*LEFT), vla.animate.shift(2.2*LEFT), vli.animate.shift(2.2*LEFT))
        self.wait()


class varnear(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta= - PI/3)
        thel = Tex(r"$\theta : (-\varepsilon, \varepsilon ) \times [a,b] \to \Sigma $ is a ", r"variation", r" of $\gamma : [a,b] \to \Sigma$ if", font_size = 34).shift(2.4*UP + 2.2*LEFT)
        thel[1].set_color(BLUE)
        thel2 = Tex(r"$\theta (0,t) = \gamma (t) $ for all $t \in [a,b]$.", font_size = 34).next_to(thel, DOWN)
        vla = Tex(r"$V ( t)$", r" $   : = \theta _ s (0,t) $ is called the ", r"direction", r" of $\theta$.", font_size = 34).next_to(thel2, DOWN)
        vla[2].set_color(BLUE)
        vli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(vla[0], DOWN).shift(0.15*UP + 2.2*LEFT)
        vvf = Tex(r"$V(t) \in T_{\gamma (t)} \Sigma$ for all $t \in [a,b]$", font_size = 34).next_to(vla, RIGHT).shift(0.85*RIGHT + 0.4*DOWN)
        box = SurroundingRectangle(vvf, buff = .1, color = YELLOW)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,2,10], x_length = 5, y_length = 5, z_length = 2, tips = False)
        axd = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 5, y_length = 5, z_length = 0.01, tips = False)
        ax.shift([3, sqrt(3), -1])
        axd.shift([-3, -sqrt(3), -1.7])
        def f(u,v):
            return (np.sin(u) + np.cos(v))/1.2
        dom = Surface(
            lambda u, v: axd.c2p( v , u , 0 ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7, 
            resolution=(3,6)
        )
        thes = Surface(
            lambda u, v: ax.c2p( v , u -1  , f(v,u -1 ) ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7,
            resolution=(3,24)
        )
        s = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v)),
            u_range = [-4,4],
            v_range = [-4,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        str = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v)),
            u_range = [-4,4],
            v_range = [-4,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.2,
            stroke_width = 0.1
        )
        v = Arrow(ax.c2p(3,-0.88 ,f(3,-1) + np.sin(1)/10 ), ax.c2p( 3, -2.92 , f(3,-1) - 1.6*np.sin(1)  ) , color = YELLOW, max_tip_length_to_length_ratio=0.2)
        vd = Arrow(axd.c2p( 3 , 0.12 , 0 ), axd.c2p( 3, -2.8 , 0  ) , color = YELLOW, max_tip_length_to_length_ratio=0.2)  
        thar = CurvedArrow( [-1.4 ,-0.5 ,0],[0.2 ,-0.2,0], radius = -3 , tip_length = 0.2)
        tharl = Tex(r"$\theta$", font_size = 34).next_to(thar, UP).shift(0.2*DOWN)
        g0 = ParametricFunction(lambda t : axd.c2p( t ,  0 , 0  ) , t_range = [-3,3], color = BLUE)
        g1 = ParametricFunction(lambda t : ax.c2p(  t, -1 , f(t,-1)  ) , t_range = [-3,3], color = BLUE)
        ga0 = ParametricFunction(lambda t : axd.c2p( t ,  -0.6 , 0  ) , t_range = [-3,3], color = PINK, stroke_width = 2)
        ga1 = ParametricFunction(lambda t : ax.c2p(  t, -1.6 , f(t,-1.6)  ) , t_range = [-3,3], color = PINK, stroke_width = 2)
        gb0 = ParametricFunction(lambda t : axd.c2p( t ,  -0.3 , 0  ) , t_range = [-3,3], color = GREEN, stroke_width = 2)
        gb1 = ParametricFunction(lambda t : ax.c2p(  t, -1.3 , f(t,-1.3)  ) , t_range = [-3,3], color = GREEN, stroke_width = 2)
        gc0 = ParametricFunction(lambda t : axd.c2p( t ,  0.3 , 0  ) , t_range = [-3,3], color = ORANGE, stroke_width = 2)
        gc1 = ParametricFunction(lambda t : ax.c2p(  t, -0.7 , f(t,-0.7)  ) , t_range = [-3,3], color = ORANGE, stroke_width = 2)
        gd0 = ParametricFunction(lambda t : axd.c2p( t , 0.6 , 0  ) , t_range = [-3,3], color = RED, stroke_width = 2)
        gd1 = ParametricFunction(lambda t : ax.c2p(  t, -0.4 , f(t,-0.4)  ) , t_range = [-3,3], color = RED, stroke_width = 2)
        self.add_fixed_in_frame_mobjects(thel, thel2, thar, tharl, vla, vli, vvf, box)
        self.add(g0,g1, ax, axd, dom, thes, str, v, vd)
        self.wait()
        self.play(Create(ga0), Create(ga1), run_time = 0.7)
        self.play(Create(gb0), Create(gb1), run_time = 0.7)
        self.play(Create(gc0), Create(gc1), run_time = 0.7)
        self.play(Create(gd0), Create(gd1), run_time = 0.7)
        self.wait()
        self.play(FadeOut(ga0), FadeOut(ga1), FadeOut(gb0), FadeOut(gb1), FadeOut(gc0), FadeOut(gc1), FadeOut(gd0), FadeOut(gd1), FadeOut(thel), FadeOut(thel2), FadeOut(vla), FadeOut(vli), FadeOut(vvf), FadeOut(box), FadeOut(thes), FadeOut(v) , FadeOut(vd), Transform(str, s))
        self.wait()

class propde(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta= - PI/3)
        prop = Tex(r"A variation $\theta$ is ", r"proper", r" if $\theta (s,a) = \gamma (a), \theta (s,b) = \gamma (b) $ for all $s \in ( - \varepsilon , \varepsilon )$", font_size = 34).shift(2.4*UP)
        prop[1].set_color(BLUE)
        l = prop.get_left()[0]
        vprop = Tex(r"For a proper variation, $V(a ) = V ( b) = 0$.", font_size = 34).next_to(prop, DOWN)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,2,10], x_length = 5, y_length = 5, z_length = 2, tips = False)
        axd = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 5, y_length = 5, z_length = 0.01, tips = False)
        ax.shift([3, sqrt(3), -1])
        axd.shift([-3, -sqrt(3), -1.7])
        def f(u,v):
            return (np.sin(u) + np.cos(v))/1.2
        dom = Surface(
            lambda u, v: axd.c2p( v , u , 0 ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7, 
            resolution=(3,6)
        )
        thes = Surface(
            lambda u, v: ax.c2p(  v , -1 + u * np.cos(PI * v /6)  , f(v, -1 + u  * np.cos(PI * v /6)  ) ),
            u_range = [-0.6,0.6],
            v_range = [-3,3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7,
            resolution=(3,24)
        )
        s = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v)),
            u_range = [-4,4],
            v_range = [-4,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        str = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v)),
            u_range = [-4,4],
            v_range = [-4,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.2,
            stroke_width = 0.1
        )
        domo = dom.copy()
        #v = Arrow(ax.c2p(3,-1 ,f(3,-1) ), ax.c2p( 3, -1 , f(3,-1) ) , color = YELLOW, max_tip_length_to_length_ratio=0.2)
        #vd = Arrow(axd.c2p( 3 , 0.12 , 0 ), axd.c2p( 3, -2.8 , 0  ) , color = YELLOW, max_tip_length_to_length_ratio=0.2)  
        thar = CurvedArrow( [-1.4 ,-0.5 ,0],[0.2 ,-0.2,0], radius = -3 , tip_length = 0.2)
        tharl = Tex(r"$\theta$", font_size = 34).next_to(thar, UP).shift(0.2*DOWN)
        g0 = ParametricFunction(lambda t : axd.c2p( t ,  0 , 0  ) , t_range = [-3,3], color = BLUE)
        g1 = ParametricFunction(lambda t : ax.c2p(  t, -1 , f(t,-1)  ) , t_range = [-3,3], color = BLUE)
        ga0 = ParametricFunction(lambda t : axd.c2p( t ,  -0.6 , 0  ) , t_range = [-3,3], color = PINK, stroke_width = 2)
        ga1 = ParametricFunction(lambda t : ax.c2p(  t, -1 - 0.6* np.cos(PI*t/6) , f(t,-1 - 0.6* np.cos(PI*t/6))  ) , t_range = [-3,3], color = PINK, stroke_width = 2)
        gb0 = ParametricFunction(lambda t : axd.c2p( t ,  -0.3 , 0  ) , t_range = [-3,3], color = GREEN, stroke_width = 2)
        gb1 = ParametricFunction(lambda t : ax.c2p(  t, -1 -0.3* np.cos(PI*t/6) , f(t,-1 - 0.3* np.cos(PI*t/6))  ) , t_range = [-3,3], color = GREEN, stroke_width = 2)
        gc0 = ParametricFunction(lambda t : axd.c2p( t ,  0.3 , 0  ) , t_range = [-3,3], color = ORANGE, stroke_width = 2)
        gc1 = ParametricFunction(lambda t : ax.c2p(  t, -1 + 0.3 * np.cos(PI*t/6) , f(t,-1 + 0.3* np.cos(PI*t/6))  ) , t_range = [-3,3], color = ORANGE, stroke_width = 2)
        gd0 = ParametricFunction(lambda t : axd.c2p( t , 0.6 , 0  ) , t_range = [-3,3], color = RED, stroke_width = 2)
        gd1 = ParametricFunction(lambda t : ax.c2p(  t, -1 + 0.6 * np.cos(PI*t/6) , f(t,-1 + 0.6* np.cos(PI*t/6))  ) , t_range = [-3,3], color = RED, stroke_width = 2)
        self.add_fixed_in_frame_mobjects(thar, tharl, prop,  vprop)
        self.add(g0,g1, s, ax, axd, dom)
        self.remove(prop,  vprop)
        self.wait()
        self.play(Create(prop),  Transform(s, str), Transform(domo, thes))
        self.wait()
        self.play(Create(ga0), Create(ga1), run_time = 0.7)
        self.play(Create(gb0), Create(gb1), run_time = 0.7)
        self.play(Create(gc0), Create(gc1), run_time = 0.7)
        self.play(Create(gd0), Create(gd1), run_time = 0.7)
        self.wait()
        self.play(Create(vprop))
        self.wait()





class fvf(Scene):
    def construct(self):        
        th = Tex(r"First variation formula", font_size=34, color = BLUE).to_edge(UL).shift(0.3*DOWN + 0.2*RIGHT)
        l = th.get_left()[0]
        th1 = Tex(r"For a proper variation $\theta : (- \varepsilon, \varepsilon ) \times [a,b] \to \Sigma$ of a curve $\gamma : [a,b] \to \Sigma $,", font_size=34).next_to(th,DOWN)
        th2 = Tex(r"we define $E : (-\varepsilon , \varepsilon ) \to \mathbb{R} $ as ", font_size=34).next_to(th1, DOWN)
        th3 = Tex(r"$E ( s) : = E( t \mapsto \theta (s,t)) = \frac{1}{2} \int_a^b \vert \theta_t  (s,t)\vert ^2 dt.$", font_size=34).next_to(th2,DOWN)
        th4 = Tex(r"Then ", font_size=34).next_to(th3,DOWN).shift(0.1*UP)
        th5 = Tex(r"$E^{\prime} (0) = - \int_a^b V(t) \cdot \gamma ^{\prime \prime} (t) $d$t$.", font_size=34).next_to(th4, DOWN)
        pr = Tex(r"Proof", font_size=34, color = BLUE).next_to(th5, DOWN)
        pr1 = Tex(r"$E^{\prime} (0) $ ", r"$=  \frac{1}{2} \int_a^b \frac{\partial }{\partial s} ( \vert \theta_t \vert ^2 ) (0,t)  $d$t$", r" $= \int_a^b ( \theta_t \cdot \theta_{st} ) (0,t) $d$t$", font_size=34).next_to(pr, DOWN)
        pr2 = Tex(r"$=   \int_a^b \frac{\partial}{\partial t} (\theta_t \cdot \theta_{s}  )(0,t)  $d$ t$ $-$ $ \int_a^b  (\theta_{tt} \cdot \theta_{s}  )(0,t)  $d$ t$  ", font_size=34).next_to(pr1, DOWN)
        pr3 = Tex(r"$= $ ", r"$   (\theta_t \cdot \theta_{s}  )(0,b)$", r" $ -$ ", r"$ (\theta_t \cdot \theta_{s}  )(0,a) $", r" $-  \int_a^b  \gamma ^{\prime \prime}(t) \cdot V (t)  $d$ t$.", font_size=34).next_to(pr2, DOWN)
        cancel = Line([pr3[1].get_left()[0],pr3[1].get_left()[1] - 0.2 , 0] , [pr3[1].get_right()[0],pr3[1].get_right()[1] + 0.2 , 0], color = RED )
        th1.shift((l - th1.get_left()[0])*RIGHT)
        th2.shift((l-th2.get_left()[0])*RIGHT)
        th3.shift(th3.get_center()[0]*LEFT)
        th4.shift((l - th4.get_left()[0])*RIGHT)
        th5.shift(th5.get_center()[0]*LEFT)
        th5b = SurroundingRectangle(th5, buff = .1, color = BLUE)
        pr.shift((l - pr.get_left()[0])*RIGHT)
        pr1.shift(pr1.get_center()[0]*LEFT)
        pr2.shift((pr1[1].get_left()[0] - pr2.get_left()[0])*RIGHT)
        pr3.shift((pr2.get_left()[0] - pr3.get_left()[0])*RIGHT)
        cancel = Line([pr3[1].get_left()[0],pr3[1].get_left()[1] - 0.1 , 0] , [pr3[1].get_right()[0],pr3[1].get_right()[1] + 0.1 , 0], color = RED )
        cancel2 = Line([pr3[3].get_left()[0],pr3[3].get_left()[1] - 0.1 , 0] , [pr3[3].get_right()[0],pr3[3].get_right()[1] + 0.1 , 0], color = RED )
        self.play(Write(th), Write(th1), Write(th2), Write(th3))
        self.wait()
        self.play(Write(th4), Write(th5), Create(th5b))
        self.wait()
        self.play(Write(pr), Write(pr1[0]), Write(pr1[1]))
        self.wait()
        self.play(Write(pr1[2]))
        self.wait()
        self.play(Write(pr2))
        self.wait()
        self.play(Write(pr3))
        self.wait()
        self.play(Create(cancel), Create(cancel2))
        self.wait()





class fvfbend(ThreeDScene):
    def construct(self):  
        self.set_camera_orientation(phi=75 * DEGREES, theta= 45*DEGREES)      
        th = Tex(r"First variation formula", font_size=34, color = BLUE).to_edge(UL).shift(0.3*DOWN + 0.2*RIGHT)
        l = th.get_left()[0]
        th1 = Tex(r"For a proper variation $\theta : (- \varepsilon, \varepsilon ) \times [a,b] \to \Sigma$ of a curve $\gamma : [a,b] \to \Sigma $,", font_size=34).next_to(th,DOWN)
        th2 = Tex(r"we define $E : (-\varepsilon , \varepsilon ) \to \mathbb{R} $ as ", font_size=34).next_to(th1, DOWN)
        th3 = Tex(r"$E ( s) : = E( t \mapsto \theta (s,t)) = \frac{1}{2} \int_a^b \vert \theta_t  (s,t)\vert ^2 dt.$", font_size=34).next_to(th2,DOWN)
        th4 = Tex(r"Then ", font_size=34).next_to(th3,DOWN).shift(0.1*UP)
        th5 = Tex(r"$E^{\prime} (0) = - \int_a^b V(t) \cdot \gamma ^{\prime \prime} (t) $d$t$.", font_size=34).next_to(th4, DOWN)
        th1.shift((l - th1.get_left()[0])*RIGHT)
        th2.shift((l-th2.get_left()[0])*RIGHT)
        th3.shift(th3.get_center()[0]*LEFT)
        th4.shift((l - th4.get_left()[0])*RIGHT)
        th5.shift(th5.get_center()[0]*LEFT)
        th5b = SurroundingRectangle(th5, buff = .1, color = BLUE)
        self.add_fixed_in_frame_mobjects(th, th1, th2, th3, th4, th5, th5b)
        self.wait()
        self.play(FadeOut(th), FadeOut(th1), FadeOut(th2), FadeOut(th3), FadeOut(th4), th5.animate.shift(2*UP), th5b.animate.shift(2*UP))
        self.wait()
        h = -2
        r = 2.5
        s = Surface(
            lambda u, v: [ r* np.sin(v)*np.cos(u) , r*np.sin(v)*np.sin(u) ,r* np.cos(v) + h ],
            u_range = [-PI,PI],
            v_range = [0,PI/2],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.2, 
            stroke_width = 0.1
        )
        def f(u,v):
            return sqrt( r**2 - u**2 - v**2 ) + h
        p = Sphere(radius = 0.07)
        q = p.copy()
        p.set_color(ORANGE).shift([0, 1,  f(0,1) ])
        q.set_color(ORANGE).shift([0, -1, f(0,-1) ])
        g = ParametricFunction(lambda t : [ sqrt(1 - t**2) , t , f( sqrt(1 - t**2) , t) ] , t_range = [-1,1], color = BLUE)
        self.play(Create(s), Create(p), Create(q), Create(g))
        self.wait()
        z = ValueTracker(1)
        g.add_updater(lambda x : x.become( ParametricFunction(lambda t : [ z.get_value()* sqrt(1 - t**2) , z.get_value() *t + (1 - z.get_value())*t**3  , f( z.get_value()* sqrt(1 - t**2), z.get_value() *t + (1 - z.get_value())*t**3) ] , t_range = [-1,1], color = BLUE)  ))
        self.play(z.animate.set_value(0))




class gede(Scene):
    def construct(self):        
        ged = Tex(r"Definition", font_size=34, color = BLUE).to_edge(UL).shift(0.6*DOWN +  RIGHT)
        l = ged.get_left()[0]
        ged1 = Tex(r"A $C^2$-curve $\gamma : [a,b] \to \Sigma$ is called a ", r"geodesic", r" if ", font_size=34).next_to(ged,DOWN)
        ged1[1].set_color(BLUE)
        ged2 = Tex(r"$ \gamma ^{\prime \prime} (t) \perp T_{\gamma (t) } \Sigma $ for all $t \in [a,b]$.", font_size=34).next_to(ged1, DOWN)    
        fac = Tex(r"Facts", font_size=34, color = BLUE).next_to(ged2, DOWN).shift(0.5*DOWN)
        fac1 = Tex(r"$\bullet$ If $\gamma$ is a geodesic, then it has constant speed.", font_size=34).next_to(fac,DOWN)
        fac2 = Tex(r"$\bullet$ If $\gamma$ minimizes the energy among smooth curves with the same  ", font_size=34).next_to(fac1, DOWN)
        fac3 = Tex(r"doimain and endpoints, then $\gamma$ is a geodesic.", font_size=34).next_to(fac2, DOWN)
        fac4 = Tex(r"$\bullet$ If $\gamma$ has constant speed and minimizes the length among smooth ", font_size=34).next_to(fac3, DOWN)
        fac5 = Tex(r"curves with the same endpoints, then $\gamma$ is a geodesic.", font_size=34).next_to(fac4, DOWN)
        ged1.shift((l - ged1.get_left()[0])*RIGHT)
        ged2.shift(ged2.get_center()[0]*LEFT)
        fac.shift((l - fac.get_left()[0])*RIGHT)
        fac1.shift((l - fac1.get_left()[0])*RIGHT)
        fac2.shift((l - fac2.get_left()[0])*RIGHT)
        fac3.shift((l - fac3.get_left()[0] + 0.4)*RIGHT)
        fac4.shift((l - fac4.get_left()[0])*RIGHT)
        fac5.shift((l - fac5.get_left()[0] + 0.4)*RIGHT)
        self.play(Create(ged), Create(ged1), Create(ged2))
        self.wait()
        self.play(Create(fac), Create(fac1))
        self.wait()
        self.play(Create(fac2), Create(fac3))
        self.wait()
        self.play(Create(fac4), Create(fac5))
        self.wait()
        



class geok(ThreeDScene):
    def construct(self):  
        self.set_camera_orientation(phi=80 * DEGREES, theta= -PI/3 )      
        gm = Tex(r"If $\Sigma$ is oriented, there is ", r"$N$", r" $ : \Sigma \to \mathbb{S}^2$ with", font_size=34).to_edge(UL).shift(RIGHT +   DOWN)
        l = gm.get_left()[0]
        nli = Line([0,0,0], [0.2,0,0], color = RED).next_to(gm[1], DOWN).shift(0.15*UP)
        gm1 = Tex(r"$N (p) \perp T_p \Sigma $ for all $p \in \Sigma$.", font_size=34).next_to(gm, DOWN)
        gk = Tex(r"For $\gamma : [a,b] \to \Sigma$ parametrized by arc length, its ", r"geodesic curvature", r" is ", font_size=34).next_to(gm1, DOWN)
        gk[1].set_color(BLUE)
        gk1 = Tex(r"$ k _{\gamma} (t) : = \gamma ^{\prime \prime } (t) \cdot (N (\gamma (t) )\times \gamma ^{\prime}(t) ) $", font_size=34).next_to(gk, DOWN)
        gm1.shift(gm1.get_center()[0]*LEFT)
        gk.shift((l - gk.get_left()[0])*RIGHT)
        gk1.shift(gk1.get_center()[0]*LEFT)
        gkb = SurroundingRectangle(gk1, buff = .1, color = BLUE)
        kneg = Tex(r"$k_{\gamma} < 0$", font_size=34).next_to(gk1, DOWN).shift(4.5*LEFT + 0.5*DOWN)
        kzer = Tex(r"$k_{\gamma} = 0$", font_size=34).next_to(gk1, DOWN).shift(0.5*DOWN)
        kpos = Tex(r"$k_{\gamma} > 0$", font_size=34).next_to(gk1, DOWN).shift(4.5*RIGHT+ 0.5*DOWN)
        knb = SurroundingRectangle(kneg, buff = .1, color = PINK)
        kzb = SurroundingRectangle(kzer, buff = .1, color = ORANGE)
        kpb = SurroundingRectangle(kpos, buff = .1, color = YELLOW)
        h = -2.7
        len = 1.5
        s1 = Surface(
            lambda u, v: [ u - 2.2*sqrt(3) ,v -2.2 , h ],
            u_range = [-1.4,1.4],
            v_range = [-1.4,1.4],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.3, 
            stroke_width = 0.2,
            resolution=(14,14)
        )
        s2 = Surface(
            lambda u, v: [ u,v, h ],
            u_range = [-1.4,1.4],
            v_range = [-1.4,1.4],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.3, 
            stroke_width = 0.2,
            resolution=(14,14)
        )
        s3 = Surface(
            lambda u, v: [ u + 2.2* sqrt(3) ,v + 2.2 , h ],
            u_range = [-1.4,1.4],
            v_range = [-1.4,1.4],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.3, 
            stroke_width = 0.2,
            resolution=(14,14)
        )
        n1 = Arrow([- 2.2*sqrt(3) , -2.2 ,h - 0.2 ]  , [- 2.2*sqrt(3) , -2.2 , h + len ] , color = RED,  max_tip_length_to_length_ratio=0.2   )
        n2 = Arrow([ 0, 0 ,h - 0.2 ]  , [ 0, 0 , h + len ] , color = RED,  max_tip_length_to_length_ratio=0.2   )
        n3 = Arrow([ 2.2*sqrt(3) , 2.2 ,h-0.2 ]  , [ 2.2*sqrt(3) , 2.2 , h + len ] , color = RED,  max_tip_length_to_length_ratio=0.2   )
        g1 = ParametricFunction(lambda t : [ t**2 - 2.2*sqrt(3) ,  t - 2.2 , h ] , t_range = [-1,1], color = PINK)
        g2 = ParametricFunction(lambda t : [ 0 , t , h ] , t_range = [-1,1], color = ORANGE)
        g3 = ParametricFunction(lambda t : [ t**2 + 2.2*sqrt(3) ,  - t + 2.2 , h ] , t_range = [-1,1], color = YELLOW)
        g1t1 = Line([ 1- 2.2*sqrt(3), -1.2, h], [ 1- 2.2*sqrt(3) - 0.2, -1.2, h ], color = PINK)
        g1t2 = Line([ 1- 2.2*sqrt(3), -1.2, h], [  1- 2.2*sqrt(3) - 3/25, -1.2 - 4/25 , h ], color = PINK)
        g2t1 = Line([ 0, 1 , h], [ -1/sqrt(125), 1 - 2/sqrt(125) , h], color = ORANGE)
        g2t2 = Line([ 0, 1 , h], [ 1/sqrt(125) , 1 - 2/sqrt(125) , h], color = ORANGE)
        g3t1 = Line([ 1 + 2.2*sqrt(3) , 1.2 , h], [ 1 + 2.2*sqrt(3) - 0.2 , 1.2 , h ], color = YELLOW)
        g3t2 = Line([ 1 + 2.2*sqrt(3) , 1.2 , h ], [  1 + 2.2*sqrt(3) - 3/25 , 1.2 + 4/25 , h ], color = YELLOW)
        self.add_fixed_in_frame_mobjects(gm, gm1, gk, gk1, nli, gkb, kneg, kzer, kpos, knb, kzb, kpb)
        self.remove(gm, gm1, gk, gk1, nli, gkb, kneg, kzer, kpos, knb, kzb, kpb)
        self.play(Create(gm), Create(gm1), Create(nli))
        self.wait()
        self.play(Write(gk), Write(gk1), Write(gkb))
        self.wait()
        self.play(Create(s1), Create(s2), Create(s3), Create(n1), Create(n2), Create(n3))
        self.wait()
        self.play(Create(kneg), Create(g1), Create(g1t1), Create(g1t2), Create(knb))
        self.wait()
        self.play(Create(kzer), Create(g2), Create(g2t1), Create(g2t2), Create(kzb))
        self.wait()
        self.play(Create(kpos), Create(g3), Create(g3t1), Create(g3t2), Create(kpb))
        self.wait()






class vfd(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=73 * DEGREES, theta= - 3*PI/16)
        de = Tex(r"Definition", color = BLUE, font_size = 34).to_edge(UL).shift(0.8*DOWN + 0.7*RIGHT)
        l = de.get_left()[0]
        de1 = Tex(r"A curve ", r"$V$" , r" $: [a,b ] \to \mathbb{R}^3$ is called a ", r"vector field", r" over ", r"$\gamma$", r" $: [a,b] \to $ ", r"$\Sigma $",r" if " ,  font_size = 34).next_to(de, DOWN)
        de1[3].set_color(BLUE)
        de2 = Tex( r"$V(t) \in T_{\gamma (t)  } \Sigma $ for all $t \in [a,b]$", font_size = 34).next_to(de1, DOWN).shift(0.2*DOWN)
        de1.shift((l - de1.get_left()[0])*RIGHT)
        de2.shift(de2.get_center()[0]*LEFT)
        vli = Line([0,0,0], [0.2,0,0], color = YELLOW).next_to(de1[1], DOWN).shift(0.15*UP)
        gli = Line([0,0,0], [0.2,0,0], color = BLUE).next_to(de1[5], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.2,0,0], color = PURPLE).next_to(de1[7], DOWN).shift(0.15*UP)
        vfb = SurroundingRectangle(de2, buff = .1, color = YELLOW)
        s = Surface(
            lambda u, v: [ u , v , - (u**2 + v**2) / 10 -0.6 ],
            u_range = [-2, 2],
            v_range = [-2, 2],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.3, 
            stroke_width = 0.3
        )
        def base(t):
            return [  t, t , -t**2 / 5  -0.6 ]
        def tip(t):
            return [ 1.1 *  np.sin(PI * t / 2 ) , - 1.1 *  np.sin(PI * t / 2 ) , 0 ] 
        def par(t):
            return -1.5 + 0.3*t 
        g = ParametricFunction(lambda t : base(t) , t_range = [-1.5,1.5], color = BLUE) 
        arrows = []
        for j in range(0,11):
            arrows.append( Arrow(  [ base(par(j))[0] - 0.1* tip(par(j))[0] , base(par(j))[1] - 0.1* tip(par(j))[1] , base(par(j))[2] - 0.1* tip(par(j))[2] ] , [ base(par(j))[0] + tip(par(j))[0] , base(par(j))[1] + tip(par(j))[1] , base(par(j))[2] + tip(par(j))[2] ] , color = YELLOW ,  max_tip_length_to_length_ratio=0.2 ))
        self.add_fixed_in_frame_mobjects(de, de1, de2, vli, gli, sli, vfb)
        self.remove(de, de1, de2, vli, gli, sli, vfb)
        self.play( Create(de), Create(de1), Create(de2), Create(vli), Create(s), Create(g), Create(gli), Create(sli), Create(vfb)) 
        self.wait()
        self.play(Create(arrows[0]), Create(arrows[1]), Create(arrows[2]), Create(arrows[3]), Create(arrows[4]), Create(arrows[5]), Create(arrows[6]), Create(arrows[7]), Create(arrows[8]), Create(arrows[9]), Create(arrows[10]))
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/8)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.wait()
        #self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(sl), FadeOut(n2), FadeOut(sli), FadeOut(pli))
        




class code(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi= 80 * DEGREES, theta= - PI/8)
        de = Tex(r"Definition", color = BLUE, font_size = 34).to_edge(UL).shift(0.2*DOWN + 0.5*RIGHT)
        l = de.get_left()[0]
        de1 = Tex(r"If ", r"$V$", r" is a vector field along ", r"$\gamma $", r", then ", r"$V^{\prime}$", r" is not necessarily a vector field along ", r"$\gamma $.", font_size = 34).next_to(de, DOWN)
        de2 = Tex( r"The vector field " , r"$D_tV : = \pi _t ( V^{\prime}(t) )$", r" is called the ", r"covariant derivative of $V$", r",",  font_size = 34).next_to(de1, DOWN)
        de3 = Tex( r"where $\pi_t : \mathbb{R}^3 \to T_{\gamma(t)}\Sigma$ is the orthogonal projection" ,  font_size = 34).next_to(de2, DOWN)
        de2[3].set_color(BLUE)
        de1.shift((l - de1.get_left()[0])*RIGHT)
        de2.shift((l - de2.get_left()[0])*RIGHT)
        de3.shift((l - de3.get_left()[0])*RIGHT)
        vli = Line([0,0,0], [0.2,0,0], color = YELLOW).next_to(de1[1], DOWN).shift(0.15*UP)
        gli = Line([0,0,0], [0.2,0,0], color = BLUE).next_to(de1[3], DOWN).shift(0.2*UP)
        vpli = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(de1[5], DOWN).shift(0.15*UP)
        gli2 = Line([0,0,0], [0.2,0,0], color = BLUE).next_to(de1[7], DOWN).shift(0.2*UP)
        cdli = Line([0,0,0], [2.4,0,0], color = PINK).next_to(de2[1], DOWN).shift(0.2*UP)
        r = 2
        h = -1.3
        w = 1.2
        s = Surface(
            lambda u, v: [ r* np.cos(u) , r* np.sin(u) , h  + v ],
            u_range = [-PI, PI],
            v_range = [-w, w ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.3, 
            stroke_width = 0.2, 
            resolution = ( 32, 12 )
        )
        g = ParametricFunction(lambda t : [ r* np.cos(t) , r*np.sin(t), h  ] , t_range = [-PI/2,PI/2], color = BLUE) 
        def base(t):
            return [ r* np.cos(t) , r*np.sin(t), h  ]
        def tip(t):
            return [ - 0.6* (t + PI/2) * np.sin(t) , 0.6* (t + PI/2) * np.cos(t) , 0 ]
        def tip2(t):
            return [ - 0.6*  (t + PI/2) * np.cos(t) - 1.3*  np.sin(t) , - 0.6*  (t + PI/2) * np.sin(t)  + 1.3* np.cos(t) , 0 ]
        def tip3(t):
            return [ - 1.3* np.sin(t) , 1.3* np.cos(t) , 0 ]
        def par(t):
            return -PI/2 + (PI/10)*t
        def vts1(t):
            return [base(t)[0] - 0.1 * tip(t)[0], base(t)[1] - 0.1 * tip(t)[1], base(t)[2] - 0.1 * tip(t)[2]]
        def vts2(t):
            return [base(t)[0] + tip(t)[0], base(t)[1] + tip(t)[1], base(t)[2] + tip(t)[2]]
        def vps1(t):
            return [base(t)[0] - 0.1 * tip2(t)[0], base(t)[1] - 0.1 * tip2(t)[1], base(t)[2] - 0.1 * tip2(t)[2]]
        def vps2(t):
            return [base(t)[0] + tip2(t)[0], base(t)[1] + tip2(t)[1], base(t)[2] + tip2(t)[2]]
        def dtvs1(t):
            return [base(t)[0] - 0.1 * tip3(t)[0], base(t)[1] - 0.1 * tip3(t)[1], base(t)[2] - 0.1 * tip3(t)[2]]
        def dtvs2(t):
            return [base(t)[0] + tip3(t)[0], base(t)[1] + tip3(t)[1], base(t)[2] + tip3(t)[2]]
        vt = []
        vp = []
        dtv = []
        for j in range(0,11):
            vt.append( Arrow(  vts1(par(j)) , vts2(par(j)) , color = YELLOW ,  max_tip_length_to_length_ratio=0.2 ))
        for j in range(0,11):
            vp.append( Arrow(  vps1(par(j)) , vps2(par(j)) , color = ORANGE ,  max_tip_length_to_length_ratio=0.2 ))
        for j in range(0,11):
            dtv.append( Arrow(  dtvs1(par(j)) , dtvs2(par(j)) , color = PINK ,  max_tip_length_to_length_ratio=0.2 ))
        vpp = []
        for j in range(0,11):
            vpp.append( vp[j].copy() )
        self.add_fixed_in_frame_mobjects(de, de1, de2, de3, vli, gli, vpli, gli2, cdli)
        self.remove(de, de1, de2, de3, vli, gli, vpli, gli2, cdli)
        self.play(Create(de), Create(de1), Create(s), Create(g), Create(vli), Create(gli), Create(vpli), Create(gli2))
        self.wait()
        self.play(Create(vt[0]), Create(vt[1]), Create(vt[2]), Create(vt[3]), Create(vt[4]), Create(vt[5]), Create(vt[6]), Create(vt[7]), Create(vt[8]), Create(vt[9]), Create(vt[10]))
        self.wait()
        self.play(Create(vp[0]), Create(vp[1]), Create(vp[2]), Create(vp[3]), Create(vp[4]), Create(vp[5]), Create(vp[6]), Create(vp[7]), Create(vp[8]), Create(vp[9]), Create(vp[10]))
        self.wait()
        self.play(FadeOut(vt[0]), FadeOut(vt[1]), FadeOut(vt[2]), FadeOut(vt[3]), FadeOut(vt[4]), FadeOut(vt[5]), FadeOut(vt[6]), FadeOut(vt[7]), FadeOut(vt[8]), FadeOut(vt[9]), FadeOut(vt[10]))
        self.wait()
        self.play(Create(de2), Create(de3), Create(cdli))
        self.wait()
        self.play(Transform(vpp[0], dtv[0]), Transform(vpp[1], dtv[1]), Transform(vpp[2], dtv[2]), Transform(vpp[3], dtv[3]), Transform(vpp[4], dtv[4]), Transform(vpp[5], dtv[5]), Transform(vpp[6], dtv[6]), Transform(vpp[7], dtv[7]), Transform(vpp[8], dtv[8]), Transform(vpp[9], dtv[9]), Transform(vpp[10], dtv[10]) )
        self.wait()
        self.play(FadeOut(vp[0]), FadeOut(vp[1]), FadeOut(vp[2]), FadeOut(vp[3]), FadeOut(vp[4]), FadeOut(vp[5]), FadeOut(vp[6]), FadeOut(vp[7]), FadeOut(vp[8]), FadeOut(vp[9]), FadeOut(vp[10]))
        self.wait()
        #self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(sl), FadeOut(n2), FadeOut(sli), FadeOut(pli))
        



class codeg(Scene):
    def construct(self):
        de = Tex(r"Definition", color = BLUE, font_size = 34).to_edge(UL).shift(0.2*DOWN + 0.5*RIGHT)
        l = de.get_left()[0]
        de1 = Tex(r"If ", r"$V$", r" is a vector field along ", r"$\gamma $", r", then ", r"$V^{\prime}$", r" is not necessarily a vector field along ", r"$\gamma $.", font_size = 34).next_to(de, DOWN)
        de2 = Tex( r"The vector field " , r"$D_tV : = \pi _t ( V^{\prime}(t) )$", r" is called the ", r"covariant derivative of $V$", r",",  font_size = 34).next_to(de1, DOWN)
        de3 = Tex( r"where $\pi_t : \mathbb{R}^3 \to T_{\gamma(t)}\Sigma$ is the orthogonal projection" ,  font_size = 34).next_to(de2, DOWN)
        de2[3].set_color(BLUE)
        de1.shift((l - de1.get_left()[0])*RIGHT)
        de2.shift((l - de2.get_left()[0])*RIGHT)
        de3.shift((l - de3.get_left()[0])*RIGHT)
        vli = Line([0,0,0], [0.2,0,0], color = YELLOW).next_to(de1[1], DOWN).shift(0.15*UP)
        gli = Line([0,0,0], [0.2,0,0], color = BLUE).next_to(de1[3], DOWN).shift(0.2*UP)
        vpli = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(de1[5], DOWN).shift(0.15*UP)
        gli2 = Line([0,0,0], [0.2,0,0], color = BLUE).next_to(de1[7], DOWN).shift(0.2*UP)
        cdli = Line([0,0,0], [2.4,0,0], color = PINK).next_to(de2[1], DOWN).shift(0.2*UP)
        pr = Tex(r"Proposition", color = BLUE, font_size = 34).next_to(de3,DOWN).shift(1.2*DOWN)
        pr1 = Tex(r"For a curve $\gamma : [a,b] \to \Sigma$ parametrized by arc length, one has", font_size = 34).next_to(pr,DOWN)
        pr2 = Tex(r"$  \vert D_t \gamma ^{\prime} (t) \vert  = \vert k_{\gamma}(t) \vert  $ for all $t \in [a,b]$.", font_size = 34).next_to(pr1,DOWN)
        pr.shift((l-pr.get_left()[0])*RIGHT)
        pr1.shift((l-pr1.get_left()[0])*RIGHT)
        pr2.shift(pr2.get_center()[0]*LEFT)
        self.add(de, de1, de2, de3, vli, gli, vpli, gli2, cdli)
        self.wait()
        self.play(Create(pr), Create(pr1), Create(pr2))
        self.wait()




class vfdd(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=73 * DEGREES, theta= - 7*PI/16)
        de = Tex(r"Definition", color = BLUE, font_size = 34).to_edge(UL).shift(0.8*DOWN + 0.7*RIGHT)
        l = de.get_left()[0]
        de1 = Tex(r"A map ", r"$V$" , r" $: $ ", r"$\Sigma$" , r" $ \to \mathbb{R}^3$ is called a ", r"vector field", r" over ", r"$\Sigma$", r" if " ,  font_size = 34).next_to(de, DOWN)
        de1[5].set_color(BLUE)
        de2 = Tex( r"$V(p) \in T_{ p } \Sigma $ for all $p \in \Sigma $", font_size = 34).next_to(de1, DOWN).shift(0.2*DOWN)
        de1.shift((l - de1.get_left()[0])*RIGHT)
        de2.shift(de2.get_center()[0]*LEFT)
        vli = Line([0,0,0], [0.2,0,0], color = YELLOW).next_to(de1[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.2,0,0], color = PURPLE).next_to(de1[3], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.2,0,0], color = PURPLE).next_to(de1[7], DOWN).shift(0.15*UP)
        vfb = SurroundingRectangle(de2, buff = .1, color = YELLOW)
        s = Surface(
            lambda u, v: [ u , v , - (u**2 + v**2) / 10 -0.6 ],
            u_range = [-2, 2],
            v_range = [-2, 2],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.3, 
            stroke_width = 0.3
        ) 
        def ini(u,v):
            return [u,v,-(u**2 + v**2)/10 - 0.6]
        def arro(u,v):
            return [1, 0 , - u / 5 ]
        def ba(u,v):
            return [ ini(u,v)[0] - 0.1 * arro(u,v)[0],  ini(u,v)[1] - 0.1 * arro(u,v)[1], ini(u,v)[2] - 0.1 * arro(u,v)[2]  ]
        def tip(u,v):
            return [ ini(u,v)[0] + arro(u,v)[0],  ini(u,v)[1] + arro(u,v)[1], ini(u,v)[2] + arro(u,v)[2]  ]
        def pa(t):
            return - 2 + 2 * t / 3
        a = []
        vf = VGroup()
        for j in range(0,6):
            a.append( [] )
        for j in range(0,7):
            for i in range(0,6):
                a[i].append( Arrow(  ba( pa(i) , pa(j) ) , tip(pa(i), pa(j)) , color = YELLOW ,  max_tip_length_to_length_ratio=0.2 ))
                vf.add(a[i][j])
        self.add_fixed_in_frame_mobjects(de, de1, de2, vli, sli2,  sli, vfb)
        self.remove(de, de1, de2, vli, sli, sli2, vfb)
        self.play( Create(de), Create(de1), Create(de2), Create(vli), Create(s), Create(sli2), Create(sli), Create(vfb)) 
        self.wait()
        #self.play( Create(a1[0]), Create(a1[1]), Create(a1[2]), Create(a1[3]), Create(a1[4]), Create(a1[5]), Create(a1[6]), Create(a1[7]), Create(a1[8]), Create(a1[9]), Create(a1[10]) )
        self.play(Create(vf))
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/8)
        self.wait(1.5)
        self.stop_ambient_camera_rotation()
        self.wait()
        




class vfdd2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=73 * DEGREES, theta= -PI/4)
        de = Tex(r"Definition", color = BLUE, font_size = 34).to_edge(UL).shift(0.8*DOWN + 0.7*RIGHT)
        l = de.get_left()[0]
        de1 = Tex(r"A map ", r"$V$" , r" $: $ ", r"$\Sigma$" , r" $ \to \mathbb{R}^3$ is called a ", r"vector field", r" over ", r"$\Sigma$", r" if " ,  font_size = 34).next_to(de, DOWN)
        de1[5].set_color(BLUE)
        de2 = Tex( r"$V(p) \in T_{ p } \Sigma $ for all $p \in \Sigma $", font_size = 34).next_to(de1, DOWN).shift(0.2*DOWN)
        de1.shift((l - de1.get_left()[0])*RIGHT)
        de2.shift(de2.get_center()[0]*LEFT)
        vli = Line([0,0,0], [0.2,0,0], color = YELLOW).next_to(de1[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.2,0,0], color = PURPLE).next_to(de1[3], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.2,0,0], color = PURPLE).next_to(de1[7], DOWN).shift(0.15*UP)
        vfb = SurroundingRectangle(de2, buff = .1, color = YELLOW)
        s = Surface(
            lambda u, v: [ u , v , - (u**2 + v**2) / 10 -0.6 ],
            u_range = [-2, 2],
            v_range = [-2, 2],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.3, 
            stroke_width = 0.3
        ) 
        def ini(u,v):
            return [u,v,-(u**2 + v**2)/10 - 0.6]
        def arro(u,v):
            return [1, 0 , - u / 5 ]
        def ba(u,v):
            return [ ini(u,v)[0] - 0.1 * arro(u,v)[0],  ini(u,v)[1] - 0.1 * arro(u,v)[1], ini(u,v)[2] - 0.1 * arro(u,v)[2]  ]
        def tip(u,v):
            return [ ini(u,v)[0] + arro(u,v)[0],  ini(u,v)[1] + arro(u,v)[1], ini(u,v)[2] + arro(u,v)[2]  ]
        def pa(t):
            return - 2 + 2 * t / 3
        a = []
        vf = VGroup()
        for j in range(0,6):
            a.append( [] )
        for j in range(0,7):
            for i in range(0,6):
                a[i].append( Arrow(  ba( pa(i) , pa(j) ) , tip(pa(i), pa(j)) , color = YELLOW ,  max_tip_length_to_length_ratio=0.2 ))
                vf.add(a[i][j])
        vf.add(s)
        shi = 3.5
        axd = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-2,6,10], x_length = 3, y_length = 3, z_length = 0.01, tips = False)
        axd.shift( [-shi /sqrt(2), -shi /sqrt(2) , -1.7])
        dom = Surface(
            lambda u, v: axd.c2p(u,v,0),
            u_range = [-2, 2],
            v_range = [-2, 2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.3, 
            stroke_width = 0.3
        ) 
        phim = Surface( 
            lambda u, v: [ u + shi/sqrt(2) , v + shi/sqrt(2) , - (u**2 + v**2) / 10 -0.6 ],
            u_range = [-1, 1],
            v_range = [-1, 1],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.3, 
            stroke_width = 0.3
        ) 
        dom2 = dom.copy()
        phiar = CurvedArrow( [-1.1 ,-0.6 ,0],[0.5 ,-0.4,0], radius = -3 , tip_length = 0.2)
        phiarl = Tex(r"$\phi$", font_size = 34).next_to(phiar, UP).shift(0.15*DOWN)
        self.add_fixed_in_frame_mobjects(de, de1, de2, vli, sli2,  sli, vfb, phiar, phiarl)
        self.remove(phiar, phiarl)
        self.add(vf)
        self.wait()
        self.play(vf.animate.shift([shi / sqrt(2) , shi / sqrt(2), 0]))
        self.wait()
        self.play(Create(axd), Create(dom), Create(phiar), Create(phiarl))
        self.wait()
        self.play(Transform(dom2, phim))
        self.wait()
        



class codep(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=73 * DEGREES, theta= -PI/4)
        coor = Tex(r"For a vector field  ", r"$V$", r" $  : $ ", r"$\Sigma$", r" $ \to \mathbb{R}^3$, and a parametrization $\phi :$ ", r"$ U$", r" $ \to $ ", r"$\Sigma$,", font_size = 34).to_edge(UL).shift(0.5*DOWN + 0.7*RIGHT)
        l = coor.get_left()[0]
        coor2 = Tex(r"we can write ", r"$V$", r" $ : $ ", r"$ U$", r" $ \to \mathbb{R}^3$ as $V(u,v)$ on $\phi (U)$. " ,  font_size = 34).next_to(coor, DOWN)
        coor2.shift((l-coor2.get_left()[0])*RIGHT)
        vli = Line([0,0,0], [0.2, 0,0], color = YELLOW).next_to(coor[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.2, 0,0], color = PURPLE).next_to(coor[3], DOWN).shift(0.15*UP)
        uli = Line([0,0,0], [0.2, 0,0], color = GREEN).next_to(coor[5], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.2, 0,0], color = PURPLE).next_to(coor[7], DOWN).shift(0.15*UP)
        coor3 = Tex(r"$D_u V (p) : = \pi_p \left( \frac{\partial}{\partial u} V (p) \right)$" ,  font_size = 34).next_to(coor2, DOWN)
        coor4 = Tex(r"$D_v V (p) : = \pi_p \left( \frac{\partial}{\partial v} V (p) \right)$" ,  font_size = 34).next_to(coor3, DOWN)
        coor3.shift(coor3.get_center()[0]*LEFT)
        coor4.shift(coor4.get_center()[0]*LEFT)
        coor5 = Tex(r"$\gamma (t)  = (u (t) , v(t))$" ,  font_size = 34).next_to(coor3, RIGHT).shift(0.3*LEFT)
        coor6 = Tex(r"$D_t V (t) : = u^{\prime} D_uV(\gamma (t) ) + v^{\prime} D_vV(\gamma (t))  $" ,  font_size = 34).next_to(coor5, DOWN).shift(0.2*DOWN)
        gli = Line([0,0,0], [2.5, 0,0], color = BLUE).next_to(coor5, DOWN).shift(0.2*UP)
        box = SurroundingRectangle(coor6, buff = .1, color = PINK)
        s = Surface(
            lambda u, v: [ u , v , - (u**2 + v**2) / 10 -0.6 ],
            u_range = [-2, 2],
            v_range = [-2, 2],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.3, 
            stroke_width = 0.3
        ) 
        def ini(u,v):
            return [u,v,-(u**2 + v**2)/10 - 0.6]
        def arro(u,v):
            return [1, 0 , - u / 5 ]
        def ba(u,v):
            return [ ini(u,v)[0] - 0.1 * arro(u,v)[0],  ini(u,v)[1] - 0.1 * arro(u,v)[1], ini(u,v)[2] - 0.1 * arro(u,v)[2]  ]
        def tip(u,v):
            return [ ini(u,v)[0] + arro(u,v)[0],  ini(u,v)[1] + arro(u,v)[1], ini(u,v)[2] + arro(u,v)[2]  ]
        def pa(t):
            return - 2 + 2 * t / 3
        a = []
        vf = VGroup()
        for j in range(0,6):
            a.append( [] )
        for j in range(0,7):
            for i in range(0,6):
                a[i].append( Arrow(  ba( pa(i) , pa(j) ) , tip(pa(i), pa(j)) , color = YELLOW ,  max_tip_length_to_length_ratio=0.2 ))
                vf.add(a[i][j])
        vf.add(s)
        shi = 3.5
        axd = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-2,6,10], x_length = 3, y_length = 3, z_length = 0.01, tips = False)
        axd.shift( [-shi /sqrt(2), -shi /sqrt(2) , -1.7])
        dom = Surface(
            lambda u, v: axd.c2p(u,v,0),
            u_range = [-2, 2],
            v_range = [-2, 2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.3, 
            stroke_width = 0.3
        ) 
        def f(u,v):
            return - (u**2 + v**2) / 10 -0.6
        phim = Surface( 
            lambda u, v: [ u + shi/sqrt(2) , v + shi/sqrt(2) , - (u**2 + v**2) / 10 -0.6 ],
            u_range = [-1, 1],
            v_range = [-1, 1],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.3, 
            stroke_width = 0.3
        ) 
        g = ParametricFunction(lambda t : axd.c2p( t , 1.4*np.sin(2*t) , 0  ) , t_range = [-1.5,1.5], color = BLUE)
        g2 = ParametricFunction(lambda t : [ t/2 + shi /sqrt(2) , 0.7*np.sin(2*t) + shi /sqrt(2) ,  - ((t/2)**2 + (0.7*np.sin(2*t))**2) / 10 -0.6  ] , t_range = [-1.5,1.5], color = BLUE)
        vf.shift([shi / sqrt(2) , shi / sqrt(2), 0])
        dom2 = dom.copy()
        phiar = CurvedArrow( [-1.1 ,-0.6 ,0],[0.5 ,-0.4,0], radius = -3 , tip_length = 0.2)
        phiarl = Tex(r"$\phi$", font_size = 34).next_to(phiar, UP).shift(0.15*DOWN)
        self.add_fixed_in_frame_mobjects(coor, coor2, phiar, phiarl, vli, uli, sli, sli2, coor3, coor4, coor5, coor6, box, gli)
        self.add(s,dom, axd, vf, phim)
        self.remove(coor, coor2, sli, uli, sli2, vli, coor3, coor4, coor5, coor6, box, gli)
        self.wait()
        self.play(Create(coor), Create(coor2), Create(sli), Create(uli), Create(sli2), Create(vli))
        self.wait()
        self.play(Create(coor3), Create(coor4))
        self.wait()
        self.play(Create(g), Create(g2), Create(coor5), coor3.animate.shift(3.5*LEFT), coor4.animate.shift(3.5*LEFT), Create(gli))
        self.wait()
        self.play(Create(coor6), Create(box))
        self.wait()




class chrisde(Scene):
    def construct(self):
        de = Tex(r"Definition", color = BLUE, font_size = 34).to_edge(UL).shift(DOWN + RIGHT)
        l = de.get_left()[0]
        de1 = Tex(r"For a parametrization $\phi : U \to \Sigma$, the ", r"Christoffel symbols",  font_size = 34).next_to(de, DOWN)
        de2 = Tex(r"$\Gamma_{11}^1 , \Gamma_{11}^2 , \Gamma_{12}^1 , \Gamma_{12}^2 , \Gamma_{21}^1 , \Gamma_{21}^2 , \Gamma_{22}^1 , \Gamma_{22}^2  : U \to \mathbb{R}  $ are defined by", font_size = 34).next_to(de1, DOWN)
        de1[1].set_color(BLUE)
        de1.shift((l - de1.get_left()[0])*RIGHT)
        de2.shift((l - de2.get_left()[0])*RIGHT)
        c1 = Tex(r"$D_u \phi_u = \Gamma_{11}^1 \phi_u + \Gamma _{11}^2 \phi_v$", font_size = 34).next_to(de2, DOWN).shift(0.5*DOWN)
        c2 = Tex(r"$D_u \phi_v = \Gamma_{12}^1 \phi_u + \Gamma _{12}^2 \phi_v$", font_size = 34).next_to(c1, DOWN)
        c3 = Tex(r"$D_v \phi_u = \Gamma_{21}^1 \phi_u + \Gamma _{21}^2 \phi_v$", font_size = 34).next_to(c2, DOWN)
        c4 = Tex(r"$D_v \phi_v = \Gamma_{22}^1 \phi_u + \Gamma _{22}^2 \phi_v$", font_size = 34).next_to(c3, DOWN)
        c1.shift(c1.get_center()[0]*LEFT)
        c2.shift(c2.get_center()[0]*LEFT)
        c3.shift(c3.get_center()[0]*LEFT)
        c4.shift(c4.get_center()[0]*LEFT)
        self.play(Write(de), Write(de1), Write(de2), Write(c1), Write(c2), Write(c3), Write(c4))
        self.wait()
        self.play(FadeOut(de), FadeOut(de1), FadeOut(de2), FadeOut(c1), FadeOut(c2), FadeOut(c3), FadeOut(c4))
        self.wait()
        c212 = MathTex("\Gamma_{21}^2")
        c = MathTex("D_v\phi_u = \Gamma_{21}^1 \phi_u +")
        c212n = MathTex("\Gamma_{21}^2").next_to(c, RIGHT)
        cc = MathTex("\phi_v").next_to(c212n, RIGHT)
        eqc = VGroup(c, c212n, cc)
        eqc.shift(eqc.get_center()[0]*LEFT)
        self.play(Create(c212))
        self.wait()
        self.play(TransformMatchingTex(c212, c212n), Write(c), Write(cc))
        self.wait()
        self.play(FadeOut(c), FadeOut(cc), FadeOut(c212n))
        self.wait()



class chriss(Scene):
    def construct(self):
        de1 = Tex(r"Since $\phi_{uv} = \phi_{vu}$, then",  font_size = 34).to_edge(UL).shift(RIGHT+2*DOWN)
        l = de1.get_left()[0]
        de2 = Tex(r"$D_u\phi_v = D_v \phi_u$ ", r" $\Rightarrow $ ", r"$\Gamma_{12}^1 = \Gamma_{21}^1 $, $\Gamma_{12}^2 = \Gamma_{21}^2$", font_size = 34).next_to(de1, DOWN).shift(0.8*DOWN)
        de2.shift( de1.get_center()[0]*LEFT)
        line = Line([0,0,0], [3,0,0], color = PINK).next_to(de2[2], DOWN).shift(0.15*UP)
        self.play(Create(de1), Write(de2[0]))
        self.wait()
        self.play(Create(de2[1]), Create(de2[2]), Create(line))
        self.wait()




class geeq(Scene):
    def construct(self):
        gee = Tex(r"A smooth curve $\gamma$ is a geodesic if and only if",  font_size = 34).to_edge(UL).shift(0.6*DOWN)
        l = gee.get_left()[0]
        gee1 = Tex(r"$ 0 = D_t \gamma^{\prime} $ ", r"$=  D_t ( u ^{\prime} \phi _u + v^{\prime } \phi_v ) $", font_size = 34).next_to(gee, DOWN).shift(0.3*DOWN)
        gee1.shift( ( 1 +  l - gee1.get_left()[0])*RIGHT)
        gee2 = Tex(r"$=  u ^{\prime \prime } \phi _ u  + u^{\prime} D_t \phi_u  + v^{\prime \prime } \phi_v  + v^{\prime } D_t \phi_v $", font_size = 34).next_to(gee1, DOWN)
        gee2.shift(( gee1[1].get_left()[0] - gee2.get_left()[0])*RIGHT)
        gee3 = Tex(r"$=  u ^{\prime \prime } \phi _ u  + u^{\prime}  (u^{\prime} D_u \phi _u + v^{\prime } D_v \phi_u  ) + v^{\prime \prime } \phi_v  + v ^{\prime} (u^{\prime} D_u \phi _v + v^{\prime } D_v \phi_v  )  $", font_size = 34).next_to(gee2, DOWN)
        gee3.shift(( gee1[1].get_left()[0] - gee3.get_left()[0])*RIGHT)
        gee4 = Tex(r"$=  $ ", r"$u ^{\prime \prime } \phi _ u $", r" $ + u^{\prime}  (u^{\prime} ( $ ", r"$\Gamma_{11}^1 \phi_u$", r" $ + $ ", r"$\Gamma_{11}^2 \phi_v$", r" $ ) + v^{\prime } ( $ ", r"$ \Gamma_{12}^1 \phi_u$", r" $ + $ ", r"$\Gamma_{12}^2 \phi_v $", r" $)  )   $", font_size = 34).next_to(gee3, DOWN)
        gee4.shift(( gee1[1].get_left()[0] - gee4.get_left()[0])*RIGHT)
        gee5 = Tex(r"$  + $ ", r"$v^{\prime \prime } \phi_v $", r" $  + v ^{\prime} (u^{\prime} ( $ ", r"$ \Gamma_{12}^1 \phi_u $", r" $ + $ ", r"$\Gamma_{12}^2 \phi_v$", r" $ ) + v^{\prime } ( $ ", r"$ \Gamma_{22}^1 \phi_u$", r" $ + $ ", r"$\Gamma_{22}^2 \phi_v$", r" $ )  )  $", font_size = 34).next_to(gee4, DOWN)
        gee5.shift(( 0.4 + gee1[1].get_left()[0] - gee5.get_left()[0])*RIGHT)
        lu1 = Line([0,0,0], [1,0,0], color = ORANGE).next_to(gee4[1], DOWN).shift(0.15*UP)
        lu2 = Line([0,0,0], [1,0,0], color = ORANGE).next_to(gee4[3], DOWN).shift(0.15*UP)
        lu3 = Line([0,0,0], [1,0,0], color = ORANGE).next_to(gee4[7], DOWN).shift(0.15*UP)
        lu4 = Line([0,0,0], [1,0,0], color = ORANGE).next_to(gee5[3], DOWN).shift(0.15*UP)
        lu5 = Line([0,0,0], [1,0,0], color = ORANGE).next_to(gee5[7], DOWN).shift(0.15*UP)
        lv1 = Line([0,0,0], [1,0,0], color = PINK).next_to(gee5[1], DOWN).shift(0.15*UP)
        lv2 = Line([0,0,0], [1,0,0], color = PINK).next_to(gee4[5], DOWN).shift(0.15*UP)
        lv3 = Line([0,0,0], [1,0,0], color = PINK).next_to(gee4[9], DOWN).shift(0.15*UP)
        lv4 = Line([0,0,0], [1,0,0], color = PINK).next_to(gee5[5], DOWN).shift(0.15*UP)
        lv5 = Line([0,0,0], [1,0,0], color = PINK).next_to(gee5[9], DOWN).shift(0.15*UP)
        gef = Tex(r"$ u^{\prime \prime } + u^{\prime} u^{\prime} \Gamma_{11}^1  + 2 u^{\prime } v^{\prime } \Gamma_{12}^1 + v^{\prime} v^{\prime} \Gamma_{22}^1 =0 $", font_size = 34).next_to(gee5, DOWN).shift(0.3*DOWN)
        gef.shift( gef.get_center()[0]*LEFT)
        gefv = Tex(r"$ v^{\prime \prime } + u^{\prime} u^{\prime} \Gamma_{11}^2  + 2 u^{\prime } v^{\prime } \Gamma_{12}^2 + v^{\prime} v^{\prime} \Gamma_{22}^2 = 0  $", font_size = 34).next_to(gef, DOWN).shift(0.2*DOWN)
        gefv.shift( gefv.get_center()[0]*LEFT)
        ubox = SurroundingRectangle(gef, buff = .1, color = ORANGE)
        vbox = SurroundingRectangle(gefv, buff = .1, color = PINK)
        self.play(Create(gee))
        self.wait()
        self.play(Create(gee1))
        self.wait()
        self.play(Create(gee2))
        self.wait()
        self.play(Create(gee3))
        self.wait()
        self.play(Create(gee4), Create(gee5))
        self.wait()
        self.play(Create(lu1), Create(lu2), Create(lu3), Create(lu4), Create(lu5), Create(lv1), Create(lv2), Create(lv3), Create(lv4), Create(lv5))
        self.wait()
        self.play(Write(gef), Write(ubox))
        self.wait()
        self.play(Write(gefv), Write(vbox))
        self.wait()
        



class geext(Scene):
    def construct(self):
        gee = Tex(r"Theorem", color = BLUE,  font_size = 34).to_edge(UL).shift(0.6*DOWN)
        l = gee.get_left()[0]
        gee1 = Tex(r"For a parametrization $\phi : U \to \Sigma$, and $p = \phi (0)$, there is an open set $A \subset \mathbb{R}^4$", font_size = 34).next_to(gee, DOWN)
        gee1.shift( ( l - gee1.get_left()[0])*RIGHT)
        gee2 = Tex(r"with $0 \in A$ and $\varepsilon > 0 $ such that for all $(u,v,a,b) \in A$, there is a geodesic ", font_size = 34).next_to(gee1, DOWN)
        gee2.shift( ( l - gee2.get_left()[0])*RIGHT)
        gee3 = Tex(r"$\gamma_{u,v,a,b} : ( - \varepsilon , \varepsilon ) \to \Sigma $ with $\gamma_{u,v,a,b}(0) = \phi (u,v)$, $\gamma_{u,v,a,b}^{\prime}(0) = a \phi_u + b \phi_v $,", font_size = 34).next_to(gee2, DOWN)
        gee3.shift( ( l - gee3.get_left()[0])*RIGHT)
        gee4 = Tex(r"and the point $\phi_{u,v,a,b}(t)$ depends smoothly on $(u,v,a,b,t) \in A \times (-\varepsilon , \varepsilon)$.", font_size = 34).next_to(gee3, DOWN)
        gee4.shift( ( l - gee4.get_left()[0])*RIGHT)
        self.play(Write(gee), Write(gee1), Write(gee2), Write(gee3), Write(gee4))
        self.wait()
        



class geex(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= -PI/4)
        gee = Tex(r"Theorem", color = BLUE,  font_size = 34).to_edge(UL).shift(0.6*DOWN)
        l = gee.get_left()[0]
        gee1 = Tex(r"For a parametrization $\phi : U \to \Sigma$, and $p = \phi (0)$, there is an open set $A \subset \mathbb{R}^4$", font_size = 34).next_to(gee, DOWN)
        gee1.shift( ( l - gee1.get_left()[0])*RIGHT)
        gee2 = Tex(r"with $0 \in A$ and $\varepsilon > 0 $ such that for all $(u,v,a,b) \in A$, there is a geodesic ", font_size = 34).next_to(gee1, DOWN)
        gee2.shift( ( l - gee2.get_left()[0])*RIGHT)
        gee3 = Tex(r"$\gamma_{u,v,a,b} : ( - \varepsilon , \varepsilon ) \to \Sigma $ with $\gamma_{u,v,a,b}(0) = \phi (u,v)$, $\gamma_{u,v,a,b}^{\prime}(0) = a \phi_u + b \phi_v $,", font_size = 34).next_to(gee2, DOWN)
        gee3.shift( ( l - gee3.get_left()[0])*RIGHT)
        gee4 = Tex(r"and the point $\phi_{u,v,a,b}(t)$ depends smoothly on $(u,v,a,b,t) \in A \times (-\varepsilon , \varepsilon)$.", font_size = 34).next_to(gee3, DOWN)
        gee4.shift( ( l - gee4.get_left()[0])*RIGHT)
        self.add_fixed_in_frame_mobjects(gee, gee1, gee2, gee3, gee4)
        self.wait()
        s = Surface(
            lambda u, v: [ u , v , - (u**2 + v**2) / 10 -0.8 ],
            u_range = [-2, 2],
            v_range = [-2, 2],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.6, 
            stroke_width = 0.3
        )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift( [0, 0, -0.8 ] )
        a = Arrow( [-0.1, 0,-0.8], [1.5,0,-0.8] , color = YELLOW , max_tip_length_to_length_ratio=0.2 )
        a1 = Arrow( [0, 0.1, -0.8], [ 0, -1.5,-0.8] , color = YELLOW , max_tip_length_to_length_ratio=0.2 )
        a2 = Arrow( [ -0.1 , -0.5 , -1/40 - 0.8 ], [1.5, -0.5, -1/40 -0.8] , color = YELLOW , max_tip_length_to_length_ratio=0.2 )
        g = ParametricFunction(lambda t : [  t, 0, - t **2 / 10 - 0.8  ] , t_range = [-1.5,1.5], color = BLUE)
        g1 = ParametricFunction(lambda t : [ 0, -t,  - t **2 / 10 - 0.8  ] , t_range = [-1.5,1.5], color = BLUE)
        w = 9.5
        g2 = ParametricFunction(lambda t : [ t, -0.5 + t**2 /w ,  - t **2 / 10 - (-0.5 +  t**2 /w )**2/ 10 - 0.8  ] , t_range = [-1.5,1.5], color = BLUE)
        self.play(Create(s), Create(p))
        self.wait()
        self.play(Create(a))
        self.wait()
        self.play(Create(g))
        self.wait()
        self.play(FadeOut(a), FadeOut(g))
        self.wait()
        self.play(Create(a1))
        self.wait()
        self.play(Create(g1))
        self.wait()
        self.play(FadeOut(a1), FadeOut(g1))
        self.wait()
        self.play(Create(a2))
        self.wait()
        self.play(Create(g2))
        self.wait()
        self.play(FadeOut(a2), FadeOut(g2))
        self.wait()
        




class geextt(Scene):
    def construct(self):
        gee = Tex(r"Theorem", color = BLUE,  font_size = 34).to_edge(UL).shift(0.6*DOWN)
        l = gee.get_left()[0]
        gee1 = Tex(r"For a parametrization $\phi : U \to \Sigma$, and $p = \phi (0)$, there is an open set $A \subset \mathbb{R}^4$", font_size = 34).next_to(gee, DOWN)
        gee1.shift( ( l - gee1.get_left()[0])*RIGHT)
        gee2 = Tex(r"with $0 \in A$ and $\varepsilon > 0 $ such that for all $(u,v,a,b) \in A$, there is a geodesic ", font_size = 34).next_to(gee1, DOWN)
        gee2.shift( ( l - gee2.get_left()[0])*RIGHT)
        gee3 = Tex(r"$\gamma_{u,v,a,b} : ( - \varepsilon , \varepsilon ) \to \Sigma $ with $\gamma_{u,v,a,b}(0) = \phi (u,v)$, $\gamma_{u,v,a,b}^{\prime}(0) = a \phi_u + b \phi_v $,", font_size = 34).next_to(gee2, DOWN)
        gee3.shift( ( l - gee3.get_left()[0])*RIGHT)
        gee4 = Tex(r"and the point $\phi_{u,v,a,b}(t)$ depends smoothly on $(u,v,a,b,t) \in A \times (-\varepsilon , \varepsilon)$.", font_size = 34).next_to(gee3, DOWN)
        gee4.shift( ( l - gee4.get_left()[0])*RIGHT)
        self.add(gee, gee1, gee2, gee3, gee4)
        self.wait()
        co = Tex(r"Corollary", color = BLUE,  font_size = 34).next_to(gee4, DOWN).shift(0.2*DOWN)
        co1 = Tex(r"For each $p \in \Sigma$, there is a ball $B \subset T_p \Sigma$ and a map $ \text{exp}_p : B \to \Sigma$", font_size = 34).next_to(co, DOWN)
        co2 = Tex(r"such that $\text{exp}_p (v) = \gamma _v (1)$, where $ \gamma _v $ is a geodesic with $\gamma_v(0) = p$, $\gamma_v^{\prime}(0) = v$.", font_size = 34).next_to(co1, DOWN)
        co3 = Tex(r"Moreover, $\text{exp}_p$ is smooth, and $d_p \text{exp} = \text{Id}_{T_p \Sigma}$.", font_size = 34).next_to(co2, DOWN)
        co.shift( ( l - co.get_left()[0])*RIGHT)
        co1.shift( ( l - co1.get_left()[0])*RIGHT)
        co2.shift( ( l - co2.get_left()[0])*RIGHT)
        co3.shift( ( l - co3.get_left()[0])*RIGHT)
        self.play(Write(co), Write(co1), Write(co2), Write(co3))
        self.wait()
        shi = gee.get_top()[1] - co.get_top()[1]
        self.play(co.animate.shift(shi*UP), co1.animate.shift(shi*UP),co2.animate.shift(shi*UP),co3.animate.shift(shi*UP), FadeOut(gee), FadeOut(gee1), FadeOut(gee2), FadeOut(gee3), FadeOut(gee4))     
        self.wait()



class exp(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= -3*PI/16)
        co = Tex(r"Corollary", color = BLUE,  font_size = 34).to_edge(UL).shift(0.6*DOWN)
        l = co.get_left()[0]
        co1 = Tex(r"For each $p \in \Sigma$, there is a ball $B \subset T_p \Sigma$ and a map $ \text{exp}_p : B \to \Sigma$", font_size = 34).next_to(co, DOWN)
        co2 = Tex(r"such that $\text{exp}_p (v) = \gamma _v (1)$, where $ \gamma _v $ is a geodesic with $\gamma_v(0) = p$, $\gamma_v^{\prime}(0) = v$.", font_size = 34).next_to(co1, DOWN)
        co3 = Tex(r"Moreover, $\text{exp}_p$ is smooth, and $d_p \text{exp} = \text{Id}_{T_p \Sigma}$.", font_size = 34).next_to(co2, DOWN)
        co.shift( ( l - co.get_left()[0])*RIGHT)
        co1.shift( ( l - co1.get_left()[0])*RIGHT)
        co2.shift( ( l - co2.get_left()[0])*RIGHT)
        co3.shift( ( l - co3.get_left()[0])*RIGHT)
        h = -0.8
        s = Surface(
            lambda u, v: [ u , v , - (u**2 + v**2) / 10  +h ],
            u_range = [-2, 2],
            v_range = [-2, 2],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.6, 
            stroke_width = 0.3
        )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift( [0, 0, h ] )
        a = Arrow( [-0.1, 0,h], [1.5,0,h] , color = YELLOW , max_tip_length_to_length_ratio=0.2 )
        a1 = Arrow( [0, 0.1, h], [ 0, -1.5,h] , color = ORANGE , max_tip_length_to_length_ratio=0.2 )
        a2 = Arrow( [0.1/sqrt(2), 0.1/sqrt(2), h], [ -1.5 / sqrt(2), -1.5 / sqrt(2),h ] , color = PINK , max_tip_length_to_length_ratio=0.2 )
        a3 = Arrow( [-0.1/sqrt(2), - 0.1/sqrt(2), h], [ 1.5 / sqrt(2), 1.5 / sqrt(2),h ] , color = GREEN_B , max_tip_length_to_length_ratio=0.2 )
        g = ParametricFunction(lambda t : [  t, 0, - t **2 / 10  +h  ] , t_range = [ 0,1.2], color = YELLOW)
        g1 = ParametricFunction(lambda t : [ 0, -t,  - t **2 / 10  +h  ] , t_range = [ 0,1.2], color = ORANGE)
        g2 = ParametricFunction(lambda t : [ -t / sqrt(2), -t /sqrt(2) ,  - t **2 / 10   +h ] , t_range = [ 0,1.2], color = PINK)
        g3 = ParametricFunction(lambda t : [ t / sqrt(2), t /sqrt(2) ,  - t **2 / 10   +h ] , t_range = [ 0,1.2], color = GREEN_B )
        self.add_fixed_in_frame_mobjects(co, co1, co2, co3)
        self.wait()
        self.play(Create(s), Create(p))
        self.wait()
        self.play(Create(a), Create(a1), Create(a2), Create(a3))
        self.wait()
        self.play(Create(g), Create(g1), Create(g2), Create(g3))
        self.wait()





class expu(Scene):
    def construct(self):
        co = Tex(r"Corollary", color = BLUE,  font_size = 34).to_edge(UL).shift(0.6*DOWN)
        l = co.get_left()[0]
        co1 = Tex(r"For each $p \in \Sigma$, there is a ball $B \subset T_p \Sigma$ and a map $ \text{exp}_p : B \to \Sigma$", font_size = 34).next_to(co, DOWN)
        co2 = Tex(r"such that $\text{exp}_p (v) = \gamma _v (1)$, where $ \gamma _v $ is a geodesic with $\gamma_v(0) = p$, $\gamma_v^{\prime}(0) = v$.", font_size = 34).next_to(co1, DOWN)
        co3 = Tex(r"Moreover, $\text{exp}_p$ is smooth, and $d_p \text{exp} = \text{Id}_{T_p \Sigma}$.", font_size = 34).next_to(co2, DOWN)
        co.shift( ( l - co.get_left()[0])*RIGHT)
        co1.shift( ( l - co1.get_left()[0])*RIGHT)
        co2.shift( ( l - co2.get_left()[0])*RIGHT)
        co3.shift( ( l - co3.get_left()[0])*RIGHT)
        self.add(co, co1, co2, co3)
        self.wait()
        ex = Tex(r"Exercise", color = BLUE,  font_size = 34).next_to(co3, DOWN).shift(0.2*DOWN)
        ex1 = Tex(r"For any $v \in T_p \Sigma$, there is an open interval $I$ with $0 \in I$ and a smooth", font_size = 34).next_to(ex,DOWN)
        ex2 = Tex(r"geodesic $\gamma_v : I \to \Sigma$ with $\gamma_v (0 ) = p$, $\gamma_v^{\prime}(0) = v$.  Moreover, if there is ", font_size = 34).next_to(ex1,DOWN)
        ex3 = Tex(r"another geodesic $\alpha : J \to \Sigma $ for some open interval $J$ with $\alpha (0 ) = p$, $\alpha^{\prime}(0) = v$, ", font_size = 34).next_to(ex2,DOWN)
        ex4 = Tex(r"then $J \subset I$, and $\alpha ( t) = \gamma_v ( t) $ for all $t \in J$.", font_size = 34).next_to(ex3,DOWN)
        ex.shift( ( l - ex.get_left()[0])*RIGHT)
        ex1.shift( ( l - ex1.get_left()[0])*RIGHT)
        ex2.shift( ( l - ex2.get_left()[0])*RIGHT)
        ex3.shift( ( l - ex3.get_left()[0])*RIGHT)
        ex4.shift( ( l - ex4.get_left()[0])*RIGHT)
        self.play(Create(ex), Create(ex1), Create(ex2), Create(ex3), Create(ex4))
        self.wait()




class expli(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= -3*PI/16)
        co = Tex(r"$\text{exp}_p : T_p \Sigma \to \Sigma$ sends lines passing through the origin",  font_size = 34).shift(2.5*UP)
        co1 = Tex(r"to geodesics passing through $p$.", font_size = 34).next_to(co, DOWN)
        h = -1.5
        s = Surface(
            lambda u, v: [ u , v , - (u**2 + v**2) / 10  +h ],
            u_range = [-2, 2],
            v_range = [-2, 2],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.6, 
            stroke_width = 0.3
        )
        tps = Surface(
            lambda u, v: [ u , v , 0 ],
            u_range = [-2, 2],
            v_range = [-2, 2],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.6, 
            stroke_width = 0.3
        )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift( [0, 0, h ] )
        q = Sphere(radius = 0.07)
        q.set_color(ORANGE).shift( [0, 0, 0 ] )
        li = ParametricFunction(lambda t : [  t, 0, 0  ] , t_range = [ -1.2,1.2], color = YELLOW)
        li1 = ParametricFunction(lambda t : [ 0, -t, 0  ] , t_range = [  -1.2,1.2], color = ORANGE)
        li2 = ParametricFunction(lambda t : [ -t / sqrt(2), -t /sqrt(2) ,  0 ] , t_range = [  -1.2,1.2], color = PINK)
        li3 = ParametricFunction(lambda t : [ t /sqrt(2) , t /sqrt(2) ,  0 ] , t_range = [  -1.2,1.2], color = GREEN_B )
        g = ParametricFunction(lambda t : [  t, 0, - t **2 / 10  +h  ] , t_range = [  -1.2,1.2], color = YELLOW)
        g1 = ParametricFunction(lambda t : [ 0, -t,  - t **2 / 10  +h  ] , t_range = [  -1.2,1.2], color = ORANGE)
        g2 = ParametricFunction(lambda t : [ -t / sqrt(2), -t /sqrt(2) ,  - t **2 / 10   +h ] , t_range = [  -1.2,1.2], color = PINK)
        g3 = ParametricFunction(lambda t : [ t / sqrt(2), t /sqrt(2) ,  - t **2 / 10   +h ] , t_range = [  -1.2,1.2], color = GREEN_B )
        self.add_fixed_in_frame_mobjects(co, co1)
        self.wait()
        self.play(Create(co), Create(co1))
        self.wait()
        self.play(Create(s), Create(p), Create(tps), Create(q))
        self.wait()
        self.play(Create(g), Create(li))
        self.wait()
        self.play(FadeOut(g), FadeOut(li))
        self.wait()
        self.play(Create(g1), Create(li1))
        self.wait()
        self.play(FadeOut(g1), FadeOut(li1))
        self.wait()
        self.play(Create(g2), Create(li2))
        self.wait()
        self.play(FadeOut(g2), FadeOut(li2))
        self.wait()
        self.play(Create(g3), Create(li3))
        self.wait()
        self.play(FadeOut(g3), FadeOut(li3))
        self.wait()





class hrii(Scene):
    def construct(self):
        gee = Tex(r"Theorem (Hopf-Rinow II)", color = BLUE,  font_size = 34).to_edge(UL).shift(2*DOWN + 1.3*RIGHT)
        l = gee.get_left()[0]
        gee1 = Tex(r"For a surface $\Sigma \subset \mathbb{R}^3$, the following are equivalent:", font_size = 34).next_to(gee, DOWN)
        gee1.shift( ( l - gee1.get_left()[0])*RIGHT)
        gee2 = Tex(r"$\bullet$ The metric space $(\Sigma, d)$ is complete.", font_size = 34).next_to(gee1, DOWN)
        gee2.shift( ( l - gee2.get_left()[0])*RIGHT)
        gee3 = Tex(r"$\bullet$ For all $q \in \Sigma$, the map $\exp_q : T_q \Sigma \to \Sigma$ is defined.", font_size = 34).next_to(gee2, DOWN)
        gee3.shift( ( l - gee3.get_left()[0])*RIGHT)
        gee4 = Tex(r"$\bullet$ There is $p \in \Sigma$ such that the map $\exp_p : T_p \Sigma \to \Sigma$ is defined.", font_size = 34).next_to(gee3, DOWN)
        gee4.shift( ( l - gee4.get_left()[0])*RIGHT)
        self.play(Create(gee), Create(gee1), Create(gee2), Create(gee3), Create(gee4))
        self.wait()

