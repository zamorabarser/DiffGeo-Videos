from this import d
from tkinter import E
from manim import *
from numpy import sqrt
import math



class ti(Scene):
    def construct(self):
        t1 = Text("Surfaces of revolution", font_size=60).shift(0.5*UP)
        t2 = Text("Lagunov's fishbowl", font_size=60).next_to(t1, DOWN)
        self.play(Write(t1), Write(t2))
        self.wait()        


class curverev(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 0 * DEGREES, theta= 0*DEGREES, gamma = 90 * DEGREES)
        gl = Tex(r"$\gamma $", r" $: [a,b] \to \mathbb{R}^2 $", font_size = 34).shift(2.5*UP)
        li = Line([0,0,0], [0.3, 0, 0], color = ORANGE).next_to(gl[0], DOWN).shift(0.15*UP)
        g = ParametricFunction(lambda t : [ 1.5 - ( 2/3 - 2*t / 7  ) * np.cos( 2 * t + 0.2 )   , t - 0.5 , 0  ] , t_range = [-PI/2,PI/2], color =ORANGE )
        l = Line([0,-2.5 ,0], [0, 1.5,0], color = BLUE)
        ll = Tex(r"$L$", font_size = 34).shift(2.3*DOWN + 0.5*LEFT)
        lli = Line([0,0,0], [0.3, 0, 0], color = BLUE).next_to(ll, DOWN).shift(0.15*UP)
        self.wait()
        self.play(Create(gl), Create(li), Create(g))
        self.wait()
        self.play(Create(l), Write(ll), Create(lli))
        self.wait()
        

class revdef(ThreeDScene):
    def construct(self):
        gl = Tex(r"$\gamma $", r" $: [a,b] \to \mathbb{R}^2 $", font_size = 34).shift(2.5*UP)
        li = Line([0,0,0], [0.3, 0, 0], color = ORANGE).next_to(gl[0], DOWN).shift(0.15*UP)
        g = ParametricFunction(lambda t : [ 1.5 - ( 2/3 - 2*t / 7  ) * np.cos( 2 * t + 0.2 )   , t - 0.5 , 0  ] , t_range = [-PI/2,PI/2], color =ORANGE )
        l = Line([0,-2.5 ,0], [0, 1.5,0], color = BLUE)
        ll = Tex(r"$L$", font_size = 34).shift(2.3*DOWN + 0.5*LEFT)
        lli = Line([0,0,0], [0.3, 0, 0], color = BLUE).next_to(ll, DOWN).shift(0.15*UP)
        self.add(gl, li, g, l, ll, lli)
        self.wait()
        self.move_camera(phi = 70 * DEGREES, theta = -60 * DEGREES, run_time =2)
        s = Surface(
            lambda u, v: [ (1.5 - ( 2/3 - 2*v / 7  ) * np.cos( 2 * v + 0.2 ))*np.cos(u) ,  v - 0.5   , (1.5 - ( 2/3 - 2*v / 7  ) * np.cos( 2 * v + 0.2 ))*np.sin(u) ],
            u_range = [0, 2*PI],
            v_range = [ -PI/2,PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        self.play(Create(s))
        self.wait()
        self.play(FadeOut(gl), FadeOut(li), FadeOut(ll), FadeOut(lli))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/5)
        self.wait(3)
        self.stop_ambient_camera_rotation()
        self.wait()


class revphi(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = - 60 * DEGREES)
        g = ParametricFunction(lambda t : [ 1.3 - ( 2/3 - 2*t / 8  ) * np.cos( 2 * t + 0.2 )   ,  0 , t    ] , t_range = [-PI/2,PI/2], color =ORANGE )
        l = Line([0, 0 , -2], [0,0, 2], color = BLUE)
        s = Surface(
            lambda u, v: [ (1.3 - ( 2/3 - 2*v / 8  ) * np.cos( 2 * v + 0.2 ))*np.cos(u) ,  (1.3 - ( 2/3 - 2*v / 8  ) * np.cos( 2 * v + 0.2 ))*np.sin(u) , v  ],
            u_range = [0, 2*PI],
            v_range = [ -PI/2,PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        group = VGroup(g,l,s)
        group.shift([2*sqrt(3),2,0])
        gl = Tex(r"$\gamma $", r" $: [a,b] \to \mathbb{R}^2 $", font_size = 34).shift(4.5*LEFT + 1.3*UP)
        li = Line([0,0,0], [0.3, 0, 0], color = ORANGE).next_to(gl[0], DOWN).shift(0.15*UP)
        left = gl.get_left()[0]
        geq = Tex(r"$\gamma (v) = (\gamma_1(v), \gamma _ 2 (v))$", font_size = 34).next_to(gl, DOWN)
        geq.shift((left - geq.get_left()[0])*RIGHT)
        leq = Tex(r"$L$", r" $= z$-axis", font_size = 34).next_to(geq, DOWN)
        leq.shift((left - leq.get_left()[0])*RIGHT)
        leqli = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(leq[0], DOWN).shift(0.15*UP)
        phidec = Tex(r"$\phi : (0, 2\pi )\times (a,b ) \to $ ", r"$\Sigma$", font_size = 34).next_to(leq, DOWN)
        phidec.shift((left - phidec.get_left()[0])*RIGHT)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(phidec[1], DOWN).shift(0.15*UP)
        phieq = Tex(r"$\phi (u,v) = (\gamma_1(v)\cos(u), \gamma_1(v)\sin(u), \gamma_2(v) )$", font_size = 34).next_to(phidec, DOWN)
        phieq.shift((left - phieq.get_left()[0])*RIGHT)
        self.add_fixed_in_frame_mobjects(gl, li, geq, leq, leqli, phidec, sli, phieq)
        self.remove(gl, li, geq, leq, leqli, phidec, sli, phieq)
        self.play(Create(g), Create(l))
        self.wait()
        self.play(Create(gl), Create(li), Create(geq), Create(leq), Create(leqli))
        self.wait()
        self.play(Create(s))
        self.wait()
        self.play(Create(phidec), Create(phieq), Create(sli))
        self.wait()



class exsp(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = - 60 * DEGREES)
        g = ParametricFunction(lambda t : [ np.sin(t)   ,  0 , np.cos(t)    ] , t_range = [0,PI], color =ORANGE )
        l = Line([0, 0 , -2], [0,0, 2], color = BLUE)
        s = Surface(
            lambda u, v: [ ( np.sin(v) )*np.cos(u) ,  (np.sin(v))*np.sin(u) ,np.cos(v)  ],
            u_range = [0, 2*PI],
            v_range = [ 0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        group = VGroup(g,l,s)
        group.shift([2*sqrt(3),2,0])
        gl = Tex(r"$\gamma $", r" $: [ 0 , \pi ] \to \mathbb{R}^2 $", font_size = 34).shift(4.5*LEFT + 1.3*UP)
        li = Line([0,0,0], [0.3, 0, 0], color = ORANGE).next_to(gl[0], DOWN).shift(0.15*UP)
        left = gl.get_left()[0]
        geq = Tex(r"$\gamma (v) = ( \sin(v), \cos(v) )$", font_size = 34).next_to(gl, DOWN)
        geq.shift((left - geq.get_left()[0])*RIGHT)
        leq = Tex(r"$L$", r" $= z$-axis", font_size = 34).next_to(geq, DOWN)
        leq.shift((left - leq.get_left()[0])*RIGHT)
        leqli = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(leq[0], DOWN).shift(0.15*UP)
        phidec = Tex(r"$\phi : (0, 2\pi )\times (0, \pi ) \to $ ", r"$\Sigma$", font_size = 34).next_to(leq, DOWN)
        phidec.shift((left - phidec.get_left()[0])*RIGHT)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(phidec[1], DOWN).shift(0.15*UP)
        phieq = Tex(r"$\phi (u,v) = (\sin (v)\cos(u), \sin (v)\sin(u), \cos (v) )$", font_size = 34).next_to(phidec, DOWN)
        phieq.shift((left - phieq.get_left()[0])*RIGHT)
        self.add_fixed_in_frame_mobjects(gl, li, geq, leq, leqli, phidec, sli, phieq)
        self.remove(gl, li, geq, leq, leqli, phidec, sli, phieq)
        self.play(Create(gl), Create(li), Create(geq), Create(leq), Create(leqli), Create(g), Create(l))
        self.wait()
        self.play(Create(s), Create(phidec), Create(phieq), Create(sli))
        self.wait()



class exto(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = - 60 * DEGREES)
        g = ParametricFunction(lambda t : [ 2 *np.sin(t) /3 + 1.5   ,  0 , 2*np.cos(t) /3    ] , t_range = [0,2*PI], color =ORANGE )
        l = Line([0, 0 , -1.5], [0,0, 1.5], color = BLUE)
        s = Surface(
            lambda u, v: [ ( 2* np.sin(v) /3 + 1.5 )*np.cos(u) ,  ( 2* np.sin(v) /3 + 1.5)*np.sin(u) , 2*np.cos(v)/3  ],
            u_range = [0, 2*PI],
            v_range = [ 0, 2*PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        group = VGroup(g,l,s)
        group.shift([1.8*sqrt(3),1.8,1])
        gl = Tex(r"$\gamma $", r" $: [ 0 , 2 \pi ] \to \mathbb{R}^2 $", font_size = 34).shift(5.2*LEFT + 0.4*UP)
        li = Line([0,0,0], [0.3, 0, 0], color = ORANGE).next_to(gl[0], DOWN).shift(0.15*UP)
        left = gl.get_left()[0]
        geq = Tex(r"$\gamma (v) = ( 2 \sin(v) /3 + 3/2 , 2\cos(v) /3 )$", font_size = 34).next_to(gl, DOWN)
        geq.shift((left - geq.get_left()[0])*RIGHT)
        leq = Tex(r"$L$", r" $= z$-axis", font_size = 34).next_to(geq, DOWN)
        leq.shift((left - leq.get_left()[0])*RIGHT)
        leqli = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(leq[0], DOWN).shift(0.15*UP)
        phidec = Tex(r"$\phi : (0, 2\pi )\times (0, 2\pi ) \to $ ", r"$\Sigma$", font_size = 34).next_to(leq, DOWN)
        phidec.shift((left - phidec.get_left()[0])*RIGHT)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(phidec[1], DOWN).shift(0.15*UP)
        phieq = Tex(r"$\phi (u,v) = ( ( 2\sin (v)/3 + 3/2 )\cos(u), (2\sin (v) /3 + 3/2 )\sin(u), 2\cos (v) /3 )$", font_size = 34).next_to(phidec, DOWN)
        phieq.shift((left - phieq.get_left()[0])*RIGHT)
        self.add_fixed_in_frame_mobjects(gl, li, geq, leq, leqli, phidec, sli, phieq)
        self.remove(gl, li, geq, leq, leqli, phidec, sli, phieq)
        self.play(Create(gl), Create(li), Create(geq), Create(leq), Create(leqli), Create(g), Create(l))
        self.wait()
        self.play(Create(s), Create(phidec), Create(phieq), Create(sli))
        self.wait()


class mp(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = - 60 * DEGREES)
        merl = Tex(r"Meridians $ = v $-curves", font_size = 34).shift(0.5*UP + 3.5*LEFT)
        merli = Line([0,0,0], [2, 0, 0], color = ORANGE).next_to(merl, DOWN).shift(0.15*UP)
        parl = Tex(r"Parallels $ = u $-curves", font_size = 34).next_to(merl, DOWN)
        parli = Line([0,0,0], [2, 0, 0], color = YELLOW).next_to(parl, DOWN).shift(0.15*UP)
        g = ParametricFunction(lambda t : [ 1.5 - ( 2/3 - 2*t / 7  ) * np.cos( 2 * t + 0.2 )  , 0 , t  ] , t_range = [-PI/2,PI/2], color =ORANGE )
        l = Line([0 , 0 ,-2 ], [0 , 0,  2 ], color = BLUE)
        s = Surface(
            lambda u, v: [ (1.5 - ( 2/3 - 2*v / 7  ) * np.cos( 2 * v + 0.2 ))*np.cos(u) ,   (1.5 - ( 2/3 - 2*v / 7  ) * np.cos( 2 * v + 0.2 ))*np.sin(u), v   ],
            u_range = [0, 2*PI],
            v_range = [ -PI/2,PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        group = VGroup(g,l,s)
        mgroup = VGroup()
        pgroup = VGroup()
        mer = []
        par = []
        for j in range(1,16):
            mer.append( ParametricFunction(lambda t : [ (1.5 - ( 2/3 - 2*t / 7  ) * np.cos( 2 * t + 0.2 )) * np.cos( j * PI /8 )  , (1.5 - ( 2/3 - 2*t / 7  ) * np.cos( 2 * t + 0.2 )) * np.sin( j * PI /8 ) , t  ] , t_range = [-PI/2,PI/2], color =ORANGE, stroke_width = 1.5 ) )
        for j in range(0, 17):
            par.append( ParametricFunction(lambda t : [ (1.5 - (2/3 - PI* ( j/8 -1 )/7 ) * np.cos( PI* (j/8 -1) + 0.2 ) )*np.cos(t) , (1.5 - (2/3 - PI* ( j/8 -1 )/7 ) * np.cos( PI* (j/8 -1) + 0.2 ))*np.sin(t) , PI*(j/8 - 1)/2  ] , t_range = [0,2*PI], color =YELLOW , stroke_width = 1.5 ) )
        group = VGroup(g,l,s)
        for x in mer:
            group.add(x)
            mgroup.add(x)
        for x in par:
            group.add(x)
            pgroup.add(x)
        group.shift([ sqrt(3),1,0])
        self.add_fixed_in_frame_mobjects(merl, merli, parl, parli)
        self.remove(merl, merli, parl, parli)
        self.play(Create(s), Create(g), Create(l) )
        self.wait()
        self.play(Create(mgroup), Create(merl), Create(merli))
        self.wait()
        self.play(Create(pgroup), Create(parl), Create(parli))
        self.wait()


class mplane(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = - 60 * DEGREES)
        g = ParametricFunction(lambda t : [ 1.5 - ( 2/3 - 2*t / 7  ) * np.cos( 2 * t + 0.2 )  , 0 , t  ] , t_range = [-PI/2,PI/2], color =ORANGE )
        l = Line([0 , 0 ,-2 ], [0 , 0,  2 ], color = BLUE)
        line2 = Line([0 , 0 ,-2 ], [0 , 0,  2 ], color = BLUE)
        s = Surface(
            lambda u, v: [ (1.5 - ( 2/3 - 2*v / 7  ) * np.cos( 2 * v + 0.2 ))*np.cos(u) ,   (1.5 - ( 2/3 - 2*v / 7  ) * np.cos( 2 * v + 0.2 ))*np.sin(u), v   ],
            u_range = [0, 2*PI],
            v_range = [ -PI/2,PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        pl = Surface(
            lambda u, v: [ u , 0  , v  ],
            u_range = [ - 3 , 3 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.3,
            resolution = 1
        )
        g2 = ParametricFunction(lambda t : [ np.sin(t)  , 0 , np.cos(t)  ] , t_range = [0, PI], color =ORANGE )
        s2 = Surface(
            lambda u, v: [ np.sin(v)*np.cos(u) ,   np.sin(v)*np.sin(u), np.cos(v) ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        pl2 = Surface(
            lambda u, v: [ u ,  0 , v  ],
            u_range = [ - 3 , 3 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.3,
            resolution = 1
        )
        group = VGroup(g,l,s,pl)
        group2 = VGroup(g2, s2, pl2, line2)
        group2.shift([1.6*sqrt(3), 1.6, 0])
        self.play(Create(s), Create(g), Create(l) )
        self.wait()
        self.play(Create(pl))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/5)
        self.wait(2)
        self.stop_ambient_camera_rotation()
        self.wait()
        self.begin_ambient_camera_rotation(rate = -PI/5)
        self.wait(2)
        self.stop_ambient_camera_rotation()
        self.wait()
        self.play(group.animate.shift([-1.6*sqrt(3), -1.6, 0]))
        self.wait()
        self.play(Create(group2))
        self.wait()


class mpgaus(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = - 60 * DEGREES)
        g = ParametricFunction(lambda t : [ 1.5 - ( 2/3 - 2*t / 7  ) * np.cos( 2 * t + 0.2 )  , 0 , t  ] , t_range = [-PI/2,PI/2], color =ORANGE )
        l = Line([0 , 0 ,-2 ], [0 , 0,  2 ], color = BLUE)
        line2 = Line([0 , 0 ,-2 ], [0 , 0,  2 ], color = BLUE)
        s = Surface(
            lambda u, v: [ (1.5 - ( 2/3 - 2*v / 7  ) * np.cos( 2 * v + 0.2 ))*np.cos(u) ,   (1.5 - ( 2/3 - 2*v / 7  ) * np.cos( 2 * v + 0.2 ))*np.sin(u), v   ],
            u_range = [0, 2*PI],
            v_range = [ -PI/2,PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        pl = Surface(
            lambda u, v: [ u , 0  , v  ],
            u_range = [ - 3 , 3 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.3,
            resolution = 1
        )
        g2 = ParametricFunction(lambda t : [ np.sin(t)  , 0 , np.cos(t)  ] , t_range = [0, PI], color =ORANGE )
        s2 = Surface(
            lambda u, v: [ np.sin(v)*np.cos(u) ,   np.sin(v)*np.sin(u), np.cos(v) ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        pl2 = Surface(
            lambda u, v: [ u ,  0 , v  ],
            u_range = [ - 3 , 3 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.3,
            resolution = 1
        )
        def zet(t):
            return -2*(np.cos(2*t + 0.2))/7 - 2*(2/3 - 2*t / 7)*np.sin(2*t + 0.2)
        def ell(t):
            return (0.85)*sqrt(1 + (zet(t))**2)
        time = ValueTracker(-PI/2)
        def gam1(t):
            return 1.5 - ( 2/3 - 2*t / 7  ) * np.cos( 2 * t + 0.2 )
        n1 = Arrow( [ gam1(-PI/2) -1.6*sqrt(3) -1 / (ell(-PI/2) * 8)  , -1.6 , -PI/2 -zet(-PI/2) / ( ell(-PI/2) * 8)  ], [gam1(-PI/2) + 1 / ell(-PI/2) -1.6*sqrt(3) , -1.6, -PI/2 + zet(-PI/2) / ell(-PI/2)] , color = PINK)
        n2 = Arrow( [ 1.6*sqrt(3)- 1 / ( ell(-PI/2) * 8),  1.6, -zet(-PI/2) / ( ell(-PI/2) * 8) ], [ 1.6*sqrt(3) +  1 / ell(-PI/2) ,  1.6 , zet(-PI/2) / ell(-PI/2)] ,  color = PINK)
        group = VGroup(g,l,s,pl)
        group2 = VGroup(g2, s2, pl2, line2)
        group.shift([-1.6*sqrt(3), -1.6,0])
        group2.shift([1.6*sqrt(3), 1.6, 0])
        self.add(group, group2)
        self.wait()
        self.play(Create(n1), Create(n2))
        self.wait()
        n1.add_updater( lambda x : x.become(  Arrow( [ gam1(time.get_value()) -1.6*sqrt(3) -1 / (ell(time.get_value()) * 8)  , -1.6 , time.get_value() - zet(time.get_value()) / ( ell(time.get_value()) * 8)  ], [gam1(time.get_value()) + 1 / ell(time.get_value()) -1.6*sqrt(3) , -1.6, time.get_value() + zet(time.get_value()) / ell(time.get_value())] , color = PINK )      ))
        n2.add_updater( lambda x : x.become(  Arrow( [ 1.6*sqrt(3)- 1 / ( ell(time.get_value()) * 8),  1.6, -zet(time.get_value()) / ( ell(time.get_value()) * 8) ], [ 1.6*sqrt(3) +  1 / ell(time.get_value()) ,  1.6 , zet(time.get_value()) / ell(time.get_value())] ,  color = PINK)   ))
        self.play(time.animate.set_value(PI/2), run_time = 2)
        self.wait()



class mpshape(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = - 60 * DEGREES)
        slab = Tex(r"$S : T_p \Sigma \to T_p \Sigma$, ", r"$W $", r" $\in T_p \Sigma, \vert W \vert = 1$, ", font_size = 34).shift(2.45*UP)
        wli = Line([0,0,0], [0.3,0,0], color = YELLOW).next_to(slab[1], DOWN).shift(0.15*UP)
        gtex = Tex(r"$\gamma$", r" $: (- \varepsilon, \varepsilon) \to \Sigma $, $\gamma ^{\prime}(0) = W$, ", r"then $S(W) : = -  (N \circ \gamma )^{\prime}(0)$", font_size = 34).next_to(slab,DOWN)
        gli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(gtex[0], DOWN).shift(0.15*UP)
        form = Tex(r"If $\gamma$ is a meridian, ", r"$ \vert S(W) \vert = \kappa _{\gamma} (0) $", font_size = 34).next_to(gtex,DOWN)
        formbox = SurroundingRectangle(form[1], color = PINK)
        g = ParametricFunction(lambda t : [ 2 - t  , 0 ,  1 - (2-t)**2  ] , t_range = [ 0 , 3/2 ], color =ORANGE )
        s = Surface(
            lambda u, v: [ (2-v)*np.cos(u) , (2-v)*np.sin(u) , (1 - (2-v)**2)  ],
            u_range = [0, 2*PI],
            v_range = [ 0, 3/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        g2 = ParametricFunction(lambda t : [ np.sin(t)  , 0 , np.cos(t) - 1 ] , t_range = [0, PI], color =ORANGE )
        s2 = Surface(
            lambda u, v: [ np.sin(v)*np.cos(u) ,   np.sin(v)*np.sin(u), np.cos(v) -1 ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        w = Arrow([ 3/2 - 1.6*sqrt(3) + 1/30 , - 1.6 , -5/4 - 1/10 ], [ 3/2 - 1.6*sqrt(3) -1/3 , - 1.6 , -5/4 + 1 ] , color = YELLOW )
        def ex(t):
            return 4 -2*t
        def ell(t):
            return (0.85)*sqrt(1 + (ex(t))**2)
        time = ValueTracker(0.2)
        self.add_fixed_in_frame_mobjects(slab, wli, gtex, gli, form, formbox)
        self.remove(slab, wli, gtex, gli, form, formbox)
        n1 = Arrow( [ 1.8 - 0.9 / (2*ell(0.2))  -1.6*sqrt(3), -1.6 , 1 - (1.8)**2 - 1 / (8*ell(0.2)) ], [1.8 + 4*(0.9)/ ell(0.2) -1.6*sqrt(3) , -1.6 , 1 - (1.8)**2  + 1 /ell(0.2) ] , color = PINK)
        n2 = Arrow( [ - 0.9 / (2*ell(0.2)) + 1.6*sqrt(3), 1.6 , - 1 / (8*ell(0.2)) - 1 ], [  4 *(0.9) / ell(0.2) + 1.6*sqrt(3) , 1.6 ,  1 /ell(0.2) - 1 ] ,  color = PINK)
        n1.add_updater( lambda x : x.become(  Arrow( [ 2 - time.get_value() - (4 - 2 * time.get_value()) / (8 *ell(time.get_value()))  -1.6*sqrt(3), -1.6 , 1 - (2 - time.get_value())**2 - 1 / (8*ell(time.get_value())) ], [2 - time.get_value() + (4 - 2 * time.get_value())/ ell(time.get_value()) -1.6*sqrt(3) , -1.6 , 1 - (2-time.get_value())**2  + 1 /ell(time.get_value()) ] , color = PINK)    ) )
        n2.add_updater( lambda x : x.become(  Arrow( [ - (4 - 2 * time.get_value()) / (8 *ell(time.get_value()))  + 1.6*sqrt(3), 1.6 , - 1 - 1 / (8*ell(time.get_value())) ], [  (4 - 2 * time.get_value())/ ell(time.get_value()) +1.6*sqrt(3) , 1.6 , -1  + 1 /ell(time.get_value()) ] , color = PINK)    ) )
        group = VGroup(g,s)
        group2 = VGroup(g2, s2)
        group.shift([-1.6*sqrt(3), -1.6,0])
        group2.shift([1.6*sqrt(3), 1.6, 0])
        self.play(Create(s), Create(s2), Write(slab[0]))
        self.wait()
        self.play(Create(w), Write(slab[1]), Write(slab[2]), Create(wli))
        self.wait()
        self.play(Create(g), Write(gtex[0]), Write(gtex[1]), Create(gli))
        self.wait()
        self.play(Create(n1), Create(n2), Write(gtex[2]))
        self.wait()
        self.play(time.animate.set_value(1.3))
        self.wait()
        self.play(Write(form), Create(formbox))
        self.wait()
        #n1.add_updater( lambda x : x.become(  Arrow( [ gam1(time.get_value()) -1.6*sqrt(3) -1 / (ell(time.get_value()) * 8)  , -1.6 , time.get_value() - zet(time.get_value()) / ( ell(time.get_value()) * 8)  ], [gam1(time.get_value()) + 1 / ell(time.get_value()) -1.6*sqrt(3) , -1.6, time.get_value() + zet(time.get_value()) / ell(time.get_value())] , color = PINK )      ))
        #n2.add_updater( lambda x : x.become(  Arrow( [ 1.6*sqrt(3)- 1 / ( ell(time.get_value()) * 8),  1.6, -zet(time.get_value()) / ( ell(time.get_value()) * 8) ], [ 1.6*sqrt(3) +  1 / ell(time.get_value()) ,  1.6 , zet(time.get_value()) / ell(time.get_value())] ,  color = PINK)   ))
        #self.play(time.animate.set_value(PI/2), run_time = 2)
        #self.wait()


class mpeigen(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = - 60 * DEGREES)
        merl = Tex(r"Eigenspaces of $S$ are the directions", font_size = 34).shift(0.5*UP + 3.5*LEFT)
        parl = Tex(r"tangent to ", r"meridians", r" and ", r"parallels", font_size = 34).next_to(merl, DOWN)
        merli = Line([0,0,0], [1.2, 0, 0], color = ORANGE).next_to(parl[1], DOWN).shift(0.15*UP)
        parli = Line([0,0,0], [1.2, 0, 0], color = YELLOW).next_to(parl[3], DOWN).shift(0.2*UP)
        s = Surface(
            lambda u, v: [ (1.5 - ( 2/3 - 2*v / 7  ) * np.cos( 2 * v + 0.2 ))*np.cos(u) ,   (1.5 - ( 2/3 - 2*v / 7  ) * np.cos( 2 * v + 0.2 ))*np.sin(u), v   ],
            u_range = [0, 2*PI],
            v_range = [ -PI/2,PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        group = VGroup(s)
        mer = []
        par = []
        for j in range(0,16):
            mer.append( ParametricFunction(lambda t : [ (1.5 - ( 2/3 - 2*t / 7  ) * np.cos( 2 * t + 0.2 )) * np.cos( j * PI /8 )  , (1.5 - ( 2/3 - 2*t / 7  ) * np.cos( 2 * t + 0.2 )) * np.sin( j * PI /8 ) , t  ] , t_range = [-PI/2,PI/2], color =ORANGE, stroke_width = 2 ) )
        for j in range(0, 17):
            par.append( ParametricFunction(lambda t : [ (1.5 - (2/3 - PI* ( j/8 -1 )/7 ) * np.cos( PI* (j/8 -1) + 0.2 ) )*np.cos(t) , (1.5 - (2/3 - PI* ( j/8 -1 )/7 ) * np.cos( PI* (j/8 -1) + 0.2 ))*np.sin(t) , PI*(j/8 - 1)/2  ] , t_range = [0,2*PI], color =YELLOW , stroke_width = 2 ) )
        for x in mer:
            group.add(x)
        for x in par:
            group.add(x)
        group.shift([1.6* sqrt(3),1.6,0])
        self.add_fixed_in_frame_mobjects(merl, merli, parl, parli)
        self.remove(merl, merli, parl, parli)
        self.play(Create(group) , Create(merl), Create(merli), Create(parl), Create(parli))
        self.wait()




class pshape(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = - 60 * DEGREES)
        slab = Tex(r"Along a ", r"parallel", r", ", r"$N$", r" is contained in a  ", r"circle", r" of radius ", r"$r$", r" $\leq 1$.", font_size = 34).shift(2.45*UP)
        pli = Line([0,0,0], [1,0,0], color = YELLOW).next_to(slab[1], DOWN).shift(0.15*UP)
        hcli = Line([0,0,0], [0.9,0,0], color = YELLOW).next_to(slab[5], DOWN).shift(0.15*UP)
        rli = Line([0,0,0], [0.3,0,0], color = RED).next_to(slab[7], DOWN).shift(0.15*UP)
        nli = Line([0,0,0], [0.3,0,0], color = PINK).next_to(slab[3], DOWN).shift(0.15*UP)
        gtex = Tex(r"The principal curvature of that direction is given by", font_size = 34).next_to(slab,DOWN)
        form = Tex(r"$\vert k_ j \vert = \frac{r}{\gamma_1} $", font_size = 34).next_to(gtex,DOWN)
        formbox = SurroundingRectangle(form, color = PINK)
        h = 1/2
        def ell(t):
            return sqrt(4 + (np.cos(t))**2)/2
        def z(t):
            return np.cos(t)/(2*ell(t))
        par = ParametricFunction(lambda t : [  (5/4) *np.cos(t) - (1.6)*sqrt(3) ,- 1.6+ (5/4)*np.sin(t) , - PI /6 - h  ] , t_range = [ 0 , 2*PI ], color = YELLOW, stroke_width = 2 )
        s = Surface(
            lambda u, v: [ (1 - np.sin(v)/2)*np.cos(u) - (1.6)*sqrt(3) , - 1.6 + (1 - np.sin(v)/2)*np.sin(u) , v - h ],
            u_range = [0, 2*PI],
            v_range = [ -PI/2 , PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        seg = Line([ 1.6*sqrt(3), 1.6, z(-PI/6) - 2* h ] ,  [1.6*sqrt(3) + 1/ell(-PI/6), 1.6, z(-PI/6) - 2* h ] , color = RED )
        g2 = ParametricFunction(lambda t : [  np.cos(t) / ell(-PI/6) + (1.6)*sqrt(3) , 1.6 + np.sin(t) / ell(-PI/6) ,  z(-PI/6) - 2* h ] , t_range = [0, 2*PI],  color = YELLOW, stroke_width = 2 )
        s2 = Surface(
            lambda u, v: [ np.sin(v)*np.cos(u) + (1.6)*sqrt(3) , 1.6 +  np.sin(v)*np.sin(u), np.cos(v) -  2*h ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        nar = CurvedArrow( 0.7* LEFT + (0.4)* DOWN, (1.3) *RIGHT + 0.4 * DOWN , radius= -3)
        narl = Tex(r"$N$", font_size = 34).next_to(nar, UP).shift(0.1*DOWN)
        time = ValueTracker(0)
        self.add_fixed_in_frame_mobjects(slab,  gtex, form, formbox, hcli, pli, rli, nar, narl, nli)
        self.remove(slab, gtex,  form, formbox, hcli, pli, rli, nar, narl, nli)
        n1 = Arrow( [5/4 -1/(8.5*ell(-PI/6)) - (1.6)*sqrt(3) , -1.6 , - PI/6 - z(-PI/6)/8.5 - h ] , [5/4 + 1/((0.8)* ell(-PI/6)) - (1.6)*sqrt(3) , -1.6 , -PI/6 + z(-PI/6)/0.8 - h ] , color = PINK)
        n2 = Arrow( [ -1/(8.5*ell(-PI/6)) + (1.6)*sqrt(3) , 1.6 , - z(-PI/6)/8.5 -  2*h ] , [ 1/((0.8)* ell(-PI/6)) + (1.6)*sqrt(3) , 1.6 ,  + z(-PI/6)/0.8 -  2*h ] , color = PINK)
        n1.add_updater( lambda x : x.become(  Arrow(  [  (5/4) *np.cos(time.get_value()) - (1.6)*sqrt(3) - np.cos(time.get_value())/(8.5*ell(-PI/6)) ,- 1.6+ (5/4)*np.sin(time.get_value()) - np.sin(time.get_value())/(8.5*ell(-PI/6)) , - PI /6 -  z(-PI/6)/8.5 - h ]  ,   [  (5/4) *np.cos(time.get_value()) - (1.6)*sqrt(3) + np.cos(time.get_value())/(0.8*ell(-PI/6))  ,- 1.6+ (5/4)*np.sin(time.get_value()) + np.sin(time.get_value())/(0.8*ell(-PI/6)) , - PI /6 + z(-PI/6)/0.8 - h ]    , color = PINK)    ) )
        n2.add_updater( lambda x : x.become(  Arrow(  [  (1.6)*sqrt(3) - np.cos(time.get_value())/(8.5*ell(-PI/6)) ,  1.6  - np.sin(time.get_value())/(8.5*ell(-PI/6)) ,  -  z(-PI/6)/8.5 - 2* h ]  ,   [  (1.6)*sqrt(3) + np.cos(time.get_value())/(0.8*ell(-PI/6)) , 1.6 + np.sin(time.get_value())/(0.8*ell(-PI/6)) ,  z(-PI/6)/0.8  - 2* h ]    , color = PINK)    ) )
        self.play(Create(s), Create(s2), Create(g2), Create(n1), Create(n2), Create(par), Create(seg), Create(slab), Create(pli), Create(hcli), Create(rli), Create(nar), Create(narl), Create(nli))
        self.wait()
        self.play(time.animate.set_value(2*PI), run_time = 3)
        self.wait()
        self.play( Create(gtex), Create(form), Create(formbox))
        self.wait()


class compute(Scene):
    def construct(self):
        ex = Tex(r"Exercise", color = BLUE, font_size = 34).to_edge(UL).shift(1.5*RIGHT + (0.8)*DOWN)
        l = ex.get_left()[0]
        ex1 = Tex(r"Use the chart", font_size = 34).next_to(ex, DOWN)
        ex2 = Tex(r"$ \phi (u,v ) = ( \gamma _1 (v) \cos (u) , \gamma_1 (v) \sin (u) , \gamma_2 (v) ), $", font_size = 34).next_to(ex1, DOWN)
        ex3 = Tex(r"to compute $k_1$ and $k_2$ in terms of $v$, $\gamma$, and its derivatives.", font_size = 34).next_to(ex2, DOWN)
        ex4 = Tex(r"In particular, ", font_size = 34).next_to(ex3, DOWN)
        ex5 = Tex(r"$    \vert k_1 \vert , \vert k_2 \vert  \leq \max \{ \kappa _{\gamma} , 1 / \gamma_1 \} $.", font_size = 34).next_to(ex4, DOWN)
        co = Tex(r"Corollary", color = BLUE, font_size = 34).next_to(ex5, DOWN).shift(0.2*DOWN)
        co1 = Tex(r"If $\gamma$ and $L$ are such that $\kappa _{\gamma }  \leq 1$ and $d(\gamma , L ) \geq 1$, then", font_size = 34).next_to(co, DOWN)
        co2 = Tex(r"the principal curvatures of $\Sigma  ( \gamma , L ) $ have absolute value $\leq 1$.", font_size = 34).next_to(co1, DOWN)
        ex1.shift((l - ex1.get_left()[0])*RIGHT)
        ex3.shift((l - ex3.get_left()[0])*RIGHT)
        ex4.shift((l - ex4.get_left()[0])*RIGHT)
        co.shift((l - co.get_left()[0])*RIGHT)
        co1.shift((l - co1.get_left()[0])*RIGHT)
        co2.shift((l - co2.get_left()[0])*RIGHT)
        ex2.shift(ex2.get_center()[0]*LEFT)
        ex5.shift(ex5.get_center()[0]*LEFT)
        group1 = VGroup(ex, ex1, ex2, ex3, ex4, ex5)
        group2 = VGroup(co, co1, co2)
        self.play(Create(group1))
        self.wait()
        self.play(Create(group2))
        self.wait()



class moon(Scene):
    def construct(self):
        thm = Tex(r"Theorem", font_size=40, color= BLUE).to_edge(UL).shift((0.8)*RIGHT)
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


class closed(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = 30 * DEGREES)
        de = Tex(r"Definition", color = BLUE, font_size = 34).to_edge(UL).shift((1.5)* RIGHT + DOWN)
        l = de.get_left()[0]
        de1 = Tex(r"A surface ", r"$\Sigma$" , r" is ", r"closed ", r" if it is compact and has no boundary.", font_size = 34).next_to(de, DOWN)
        de1[3].set_color(BLUE)
        de1.shift((l - de1.get_left()[0])*RIGHT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE_B).next_to(de1[1], DOWN).shift(0.15*UP)
        jcs = Tex(r"Theorem", color = BLUE, font_size = 34).next_to(de1, DOWN).shift(DOWN)
        jcs1 = Tex(r"If $\Sigma$ is a closed surface, $\mathbb{R}^3 \backslash \Sigma $ has two connected components,", font_size = 34).next_to(jcs, DOWN)
        jcs2 = Tex(r"one of which is bounded, and we call it the ", r"interior", r" of $\Sigma$ .", font_size = 34).next_to(jcs1, DOWN)
        jcs.shift((l - jcs.get_left()[0])*RIGHT)
        jcs1.shift((l - jcs1.get_left()[0])*RIGHT)
        jcs2.shift((l - jcs2.get_left()[0])*RIGHT)
        jcs2[1].set_color(BLUE)
        s = Surface(
            lambda u, v: [ 2* np.sin(v)*np.cos(u) , 3 * np.sin(v)*np.sin(u) ,  np.cos(v) -1.3 ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        self.add_fixed_in_frame_mobjects(de, de1, sli, jcs, jcs1, jcs2)
        self.remove(de, de1, sli, jcs, jcs1, jcs2)
        self.wait()
        self.play(Create(s), Create(de), Create(de1), Create(sli))
        self.wait()
        self.play(FadeOut(s))
        self.wait()
        self.play(Create(jcs), Create(jcs1), Create(jcs2))
        self.wait()


class moon3(Scene):
    def construct(self):
        m = Tex(r"Theorem?", color = BLUE, font_size = 34).to_edge(UL).shift((1.5)* RIGHT + (1.3)*DOWN)
        l = m.get_left()[0]
        m1 = Tex(r"If $\Sigma$ is a closed surface, and its principal curvatures $k_1, k_2$ satisfy", font_size = 34).next_to(m, DOWN)
        m2 = Tex(r"$\vert k_1 \vert , \vert k_2 \vert \leq 1,$", font_size = 34).next_to(m1, DOWN)
        m3 = Tex(r"then $\Sigma $ contains a unit sphere in its interior.", font_size = 34).next_to(m2, DOWN)
        cros = Line([-5,0,0],[5,2,0], color = RED)
        cros2 = Line([5,0,0],[-5,2,0], color = RED)
        ex = Tex(r"Example (Lagunov)", color = BLUE, font_size = 34).next_to(m3, DOWN).shift(0.3*DOWN)
        ex1 = Tex(r"There is a closed surface with principal curvatures $\leq 1$, and", font_size = 34).next_to(ex, DOWN)
        ex2 = Tex(r"the largest sphere in its interior has radius $\leq 1/6$.", font_size = 34).next_to(ex1, DOWN)
        m1.shift((l - m1.get_left()[0])*RIGHT)
        m2.shift(m2.get_center()[0]*LEFT)
        m3.shift((l - m3.get_left()[0])*RIGHT)
        ex.shift((l - ex.get_left()[0])*RIGHT)
        ex1.shift((l - ex1.get_left()[0])*RIGHT)
        ex2.shift((l - ex2.get_left()[0])*RIGHT)
        group1 = VGroup(m, m1, m2, m3)
        group2 = VGroup(cros, cros2, ex, ex1, ex2)
        self.play(Create(group1))
        self.wait()
        self.play(Create(group2))
        self.wait()


class lag0small(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 90 * DEGREES, theta = -90 * DEGREES)
        a = Dot()
        c1 = ParametricFunction(lambda t : [  2 + np.sin(t)/2 , 0 ,  0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c2 = ParametricFunction(lambda t : [  2 + np.sin(t)/2 , 0 ,  -0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c3 = ParametricFunction(lambda t : [  2 + 0.6*np.sin(t) , 0 ,  0.6  + 0.6*np.cos(t)  ] , t_range = [ 0 , 2*PI/3 ], color = ORANGE)
        c4 = ParametricFunction(lambda t : [  2 + 0.6*sqrt(3) - 0.6*np.sin(t) , 0 ,  0.6*np.cos(t)  ] , t_range = [ PI/3 , 2*PI/3 ], color = ORANGE)        
        c5 = ParametricFunction(lambda t : [  2 + 0.6*np.sin(t) , 0 ,  -0.6  + 0.6*np.cos(t)  ] , t_range = [ PI/3 , PI ], color = ORANGE)    
        l1 = ParametricFunction(lambda t : [  t , 0 , 1.2  ] , t_range = [ 0 , 2 ], color = ORANGE)
        l2 = ParametricFunction(lambda t : [  t , 0 , 1.1  ] , t_range = [ 0 , 2 ], color = ORANGE)
        l3 = ParametricFunction(lambda t : [  t , 0 , 0.1  ] , t_range = [ 0 , 2 ], color = ORANGE)
        l4 = ParametricFunction(lambda t : [  t , 0 , -0.1  ] , t_range = [ 0 , 2 ], color = ORANGE)
        l5 = ParametricFunction(lambda t : [  t , 0 , -1.1  ] , t_range = [ 0 , 2 ], color = ORANGE)
        l6 = ParametricFunction(lambda t : [  t , 0 , -1.2  ] , t_range = [ 0 , 2 ], color = ORANGE)
        g = VGroup(c1, c2, c3, c4, c5)
        self.play(Create(g))
        self.wait()
        self.move_camera(phi = 70 * DEGREES, theta = -60 * DEGREES, run_time =2)
        


class lag0(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 90 * DEGREES, theta = -90 * DEGREES)
        c1 = ParametricFunction(lambda t : [  2.5 + np.sin(t)/2 , 0 ,  0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c2 = ParametricFunction(lambda t : [  2.5 + np.sin(t)/2 , 0 ,  -0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c3 = ParametricFunction(lambda t : [  2.5 + 0.6*np.sin(t) , 0 ,  0.6  + 0.6*np.cos(t)  ] , t_range = [ 0 , 2*PI/3 ], color = ORANGE)
        c4 = ParametricFunction(lambda t : [  2.5 + 0.6*sqrt(3) - 0.6*np.sin(t) , 0 ,  0.6*np.cos(t)  ] , t_range = [ PI/3 , 2*PI/3 ], color = ORANGE)        
        c5 = ParametricFunction(lambda t : [  2.5 + 0.6*np.sin(t) , 0 ,  -0.6  + 0.6*np.cos(t)  ] , t_range = [ PI/3 , PI ], color = ORANGE)    
        sc1 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) ,  0.6  + np.cos(v)/2    ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc2 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) , - 0.6  + np.cos(v)/2     ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc3 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,   0.6  + 0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [0, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc4 = Surface( lambda u, v: [ (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.sin(u) ,  0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [PI/3, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc5 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,  -0.6  + 0.6*np.cos(v)  ],
            u_range = [0, 2*PI],
            v_range = [PI/3, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        l1 = ParametricFunction(lambda t : [  t , 0 , 1.2  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l2 = ParametricFunction(lambda t : [  t , 0 , 1.1 ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l3 = ParametricFunction(lambda t : [  t , 0 , 0.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l4 = ParametricFunction(lambda t : [  t , 0 , -0.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l5 = ParametricFunction(lambda t : [  t , 0 , -1.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l6 = ParametricFunction(lambda t : [  t , 0 , -1.2  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        sl1 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl2 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl3 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl4 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl5 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -1.1   ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl6 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        g = VGroup(c1, c2, c3, c4, c5)
        g2 = VGroup(l1, l2, l3, l4, l5, l6)
        g3 = VGroup(sc1, sc2, sc3, sc4, sc5, sl1, sl2, sl3, sl4, sl5, sl6)
        self.play(Create(g))
        self.wait()
        self.play(Create(g2))
        self.wait()
        self.move_camera(phi = 70 * DEGREES, theta = -60 * DEGREES, run_time =2)
        self.wait()
        self.play(Create(g3))
        self.wait()


class lag1(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = -60 * DEGREES)
        c1 = ParametricFunction(lambda t : [  2.5 + np.sin(t)/2 , 0 ,  0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c2 = ParametricFunction(lambda t : [  2.5 + np.sin(t)/2 , 0 ,  -0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c3 = ParametricFunction(lambda t : [  2.5 + 0.6*np.sin(t) , 0 ,  0.6  + 0.6*np.cos(t)  ] , t_range = [ 0 , 2*PI/3 ], color = ORANGE)
        c4 = ParametricFunction(lambda t : [  2.5 + 0.6*sqrt(3) - 0.6*np.sin(t) , 0 ,  0.6*np.cos(t)  ] , t_range = [ PI/3 , 2*PI/3 ], color = ORANGE)        
        c5 = ParametricFunction(lambda t : [  2.5 + 0.6*np.sin(t) , 0 ,  -0.6  + 0.6*np.cos(t)  ] , t_range = [ PI/3 , PI ], color = ORANGE)    
        sc1 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) ,  0.6  + np.cos(v)/2    ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc2 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) , - 0.6  + np.cos(v)/2     ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc3 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,   0.6  + 0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [0, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc4 = Surface( lambda u, v: [ (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.sin(u) ,  0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [PI/3, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc5 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,  -0.6  + 0.6*np.cos(v)  ],
            u_range = [0, 2*PI],
            v_range = [PI/3, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        l1 = ParametricFunction(lambda t : [  t , 0 , 1.2  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l2 = ParametricFunction(lambda t : [  t , 0 , 1.1 ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l3 = ParametricFunction(lambda t : [  t , 0 , 0.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l4 = ParametricFunction(lambda t : [  t , 0 , -0.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l5 = ParametricFunction(lambda t : [  t , 0 , -1.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l6 = ParametricFunction(lambda t : [  t , 0 , -1.2  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        sl1 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl2 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl3 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl4 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl5 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -1.1   ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl6 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        g = VGroup(c1, c2, c3, c4, c5)
        g2 = VGroup(l1, l2, l3, l4, l5, l6)
        g3 = VGroup(sc1, sc2, sc3, sc4, sc5, sl1, sl2, sl3, sl4, sl5, sl6)
        g4 = VGroup(sc3, sc4, sc5, sl1, sl6)
        self.add(g, g2, g3)
        self.wait()
        self.play(FadeOut(g4))
        self.wait()
        self.play(FadeIn(g4))
        self.wait()


class lag2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = -60 * DEGREES)
        c1 = ParametricFunction(lambda t : [  2.5 + np.sin(t)/2 , 0 ,  0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c2 = ParametricFunction(lambda t : [  2.5 + np.sin(t)/2 , 0 ,  -0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c3 = ParametricFunction(lambda t : [  2.5 + 0.6*np.sin(t) , 0 ,  0.6  + 0.6*np.cos(t)  ] , t_range = [ 0 , 2*PI/3 ], color = ORANGE)
        c4 = ParametricFunction(lambda t : [  2.5 + 0.6*sqrt(3) - 0.6*np.sin(t) , 0 ,  0.6*np.cos(t)  ] , t_range = [ PI/3 , 2*PI/3 ], color = ORANGE)        
        c5 = ParametricFunction(lambda t : [  2.5 + 0.6*np.sin(t) , 0 ,  -0.6  + 0.6*np.cos(t)  ] , t_range = [ PI/3 , PI ], color = ORANGE)    
        sc1 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) ,  0.6  + np.cos(v)/2    ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc2 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) , - 0.6  + np.cos(v)/2     ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc3 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,   0.6  + 0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [0, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc4 = Surface( lambda u, v: [ (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.sin(u) ,  0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [PI/3, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc5 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,  -0.6  + 0.6*np.cos(v)  ],
            u_range = [0, 2*PI],
            v_range = [PI/3, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        l1 = ParametricFunction(lambda t : [  t , 0 , 1.2  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l2 = ParametricFunction(lambda t : [  t , 0 , 1.1 ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l3 = ParametricFunction(lambda t : [  t , 0 , 0.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l4 = ParametricFunction(lambda t : [  t , 0 , -0.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l5 = ParametricFunction(lambda t : [  t , 0 , -1.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l6 = ParametricFunction(lambda t : [  t , 0 , -1.2  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        sl1 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl2 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl3 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl4 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl5 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -1.1   ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl6 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        g = VGroup(c1, c2, c3, c4, c5)
        g2 = VGroup(l1, l2, l3, l4, l5, l6)
        g3 = VGroup(sc1, sc2, sc3, sc4, sc5, sl1, sl2, sl3, sl4, sl5, sl6)
        g4 = VGroup(sc1, sc2, sc3, sc4, sc5)
        self.add(g, g2, g3)
        self.wait()
        self.play(FadeOut(g4), FadeOut(g))
        self.wait()
        self.play(FadeIn(g4), FadeIn(g))
        self.wait()


class lag3(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = -60 * DEGREES)
        c1 = ParametricFunction(lambda t : [  2.5 + np.sin(t)/2 , 0 ,  0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c2 = ParametricFunction(lambda t : [  2.5 + np.sin(t)/2 , 0 ,  -0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c3 = ParametricFunction(lambda t : [  2.5 + 0.6*np.sin(t) , 0 ,  0.6  + 0.6*np.cos(t)  ] , t_range = [ 0 , 2*PI/3 ], color = ORANGE)
        c4 = ParametricFunction(lambda t : [  2.5 + 0.6*sqrt(3) - 0.6*np.sin(t) , 0 ,  0.6*np.cos(t)  ] , t_range = [ PI/3 , 2*PI/3 ], color = ORANGE)        
        c5 = ParametricFunction(lambda t : [  2.5 + 0.6*np.sin(t) , 0 ,  -0.6  + 0.6*np.cos(t)  ] , t_range = [ PI/3 , PI ], color = ORANGE)    
        sc1 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) ,  0.6  + np.cos(v)/2    ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc2 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) , - 0.6  + np.cos(v)/2     ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc3 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,   0.6  + 0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [0, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc4 = Surface( lambda u, v: [ (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.sin(u) ,  0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [PI/3, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc5 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,  -0.6  + 0.6*np.cos(v)  ],
            u_range = [0, 2*PI],
            v_range = [PI/3, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        l1 = ParametricFunction(lambda t : [  t , 0 , 1.2  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l2 = ParametricFunction(lambda t : [  t , 0 , 1.1 ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l3 = ParametricFunction(lambda t : [  t , 0 , 0.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l4 = ParametricFunction(lambda t : [  t , 0 , -0.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l5 = ParametricFunction(lambda t : [  t , 0 , -1.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l6 = ParametricFunction(lambda t : [  t , 0 , -1.2  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        sl1 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl2 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl3 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl4 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl5 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -1.1   ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl6 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        g = VGroup(c1, c2, c3, c4, c5)
        g2 = VGroup(l1, l2, l3, l4, l5, l6)
        g3 = VGroup(sc1, sc2, sc3, sc4, sc5, sl1, sl2, sl3, sl4, sl5, sl6)
        g4 = VGroup(sl1, sl2, sl3, sl4, sl5, sl6)
        self.add(g, g2, g3)
        self.wait()
        self.play(FadeOut(g4), FadeOut(g2))
        self.wait()
        self.play(FadeIn(g4), FadeIn(g2))
        self.wait()


class lag4(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = -60 * DEGREES)
        c1 = ParametricFunction(lambda t : [  2.5 + np.sin(t)/2 , 0 ,  0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c2 = ParametricFunction(lambda t : [  2.5 + np.sin(t)/2 , 0 ,  -0.6  + np.cos(t)/2  ] , t_range = [ 0 , PI ], color = ORANGE)
        c3 = ParametricFunction(lambda t : [  2.5 + 0.6*np.sin(t) , 0 ,  0.6  + 0.6*np.cos(t)  ] , t_range = [ 0 , 2*PI/3 ], color = ORANGE)
        c4 = ParametricFunction(lambda t : [  2.5 + 0.6*sqrt(3) - 0.6*np.sin(t) , 0 ,  0.6*np.cos(t)  ] , t_range = [ PI/3 , 2*PI/3 ], color = ORANGE)        
        c5 = ParametricFunction(lambda t : [  2.5 + 0.6*np.sin(t) , 0 ,  -0.6  + 0.6*np.cos(t)  ] , t_range = [ PI/3 , PI ], color = ORANGE)    
        sc1 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) ,  0.6  + np.cos(v)/2    ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc2 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) , - 0.6  + np.cos(v)/2     ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        l1 = ParametricFunction(lambda t : [  t , 0 , 1.2  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l2 = ParametricFunction(lambda t : [  t , 0 , 1.1 ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l3 = ParametricFunction(lambda t : [  t , 0 , 0.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l4 = ParametricFunction(lambda t : [  t , 0 , -0.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l5 = ParametricFunction(lambda t : [  t , 0 , -1.1  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        l6 = ParametricFunction(lambda t : [  t , 0 , -1.2  ] , t_range = [ 0 , 2.5 ], color = ORANGE)
        sl2 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl3 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl4 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl5 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -1.1   ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        spint = Sphere(
            center=(2.45 + (0.6)*(sqrt(3) - 2/sqrt(3)), 0, 0),
            radius= (0.6)*(2/sqrt(3) - 1),
            resolution=(5, 5),
            u_range=[0.001, PI - 0.001],
            v_range=[0, TAU]
        )
        spint.set_color(GREEN)
        g = VGroup(c1, c2, c3, c4, c5)
        g2 = VGroup(l1, l2, l3, l4, l5, l6)
        g3 = VGroup(sl2, sl3, sl4, sl5, sc1, sc2)
        self.add(g, g2, g3)
        self.wait()
        self.play(Create(spint))
        self.wait()
        rad = Tex(r"radius", r" $\approx \frac{2}{\sqrt{3}} -1  < \frac{1}{6}$", font_size = 34).shift(2.7*UP)
        radli = Line([0,0,0], [0.7,0,0], color = GREEN).next_to(rad[0], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(rad, radli)
        self.remove(rad, radli)
        self.play(Create(rad), Create(radli))
        self.wait()




class tube0(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 90 * DEGREES, theta = -90 * DEGREES)
        c1 = ParametricFunction(lambda t : [  2.3 - np.sin(t) , 0 ,   np.cos(t)  ] , t_range = [ 0 , PI ], color = ORANGE)
        c2 = ParametricFunction(lambda t : [  2.3 - 1.3*np.sin(t) , 0 , 1.3*np.cos(t) - 0.1  ] , t_range = [ 0 , PI ], color = ORANGE)    
        sc1 = Surface( lambda u, v: [ (2.3 - np.sin(v)) *np.cos(u) , (2.3 - np.sin(v)) *np.sin(u) ,  np.cos(v)    ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc2 = Surface( lambda u, v: [ (2.3 - 1.3*np.sin(v)) *np.cos(u) , (2.3 - 1.3*np.sin(v)) *np.sin(u) , 1.3*np.cos(v) - 0.1   ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        self.play(Create(c1), Create(c2))
        self.wait()
        self.move_camera(phi = 70 * DEGREES, theta = -60 * DEGREES, run_time =2)
        self.wait()
        self.play(Create(sc1), Create(sc2))
        self.wait()


class lag5(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = -60 * DEGREES)
        sc1 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) ,  0.6  + np.cos(v)/2    ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc2 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) , - 0.6  + np.cos(v)/2     ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc3 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,   0.6  + 0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [0, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc4 = Surface( lambda u, v: [ (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.sin(u) ,  0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [PI/3, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc5 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,  -0.6  + 0.6*np.cos(v)  ],
            u_range = [0, 2*PI],
            v_range = [PI/3, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl1 = Surface( lambda u, v: [ 1.25 + v* np.cos(u) , v*np.sin(u) ,  1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 1.075 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl2 = Surface( lambda u, v: [ 1.25 + v* np.cos(u) , v*np.sin(u) ,  1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1.075 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl3 = Surface( lambda u, v: [ 1.25 + v* np.cos(u) , v*np.sin(u) ,  0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1.075 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl4 = Surface( lambda u, v: [ 1.25 + v* np.cos(u) , v*np.sin(u) ,  -0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1.075 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl5 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -1.1   ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl6 = Surface( lambda u, v: [ v* np.cos(u) , v*np.sin(u) ,  -1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 2.5 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        ##### new floor and ceilings:
        snl1 = Surface( lambda u, v: [   v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl2 = Surface( lambda u, v: [   v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl3 = Surface( lambda u, v: [   v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl4 = Surface( lambda u, v: [   v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  -0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        ##### first tube:
        tub1 = Surface( lambda u, v: [ 1.25 +  (2.3 - np.sin(v)) *np.cos(u) /2 , (2.3 - np.sin(v)) *np.sin(u)/2 ,  np.cos(v) /2  + 0.6  ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tub2 = Surface( lambda u, v: [ 1.25 + (2.3 - 1.3*np.sin(v)) *np.cos(u) /2 , (2.3 - 1.3*np.sin(v)) *np.sin(u) /2 , 1.3*np.cos(v) /2 + 0.55   ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        g = VGroup(sc1, sc2, sc3, sc4, sc5, sl1, sl2, sl3, sl4, sl5, sl6, snl1, snl2, snl3, snl4)
        lids = VGroup(sl1, sl2, sl3, sl4)
        sides = VGroup(sc1, sc2, sc3, sc4, sc5)
        self.add(g)
        self.wait()
        self.play(FadeOut(sides))
        self.wait()
        self.play(FadeOut(lids))
        self.wait()
        self.play(Create(tub1), Create(tub2))
        self.wait()



class lag6(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = -60 * DEGREES)
        sl3 = Surface( lambda u, v: [ - 1.25 + v* np.cos(u) , v*np.sin(u) ,  0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1.075 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl4 = Surface( lambda u, v: [ - 1.25 + v* np.cos(u) , v*np.sin(u) ,  -0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1.075 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl5 = Surface( lambda u, v: [ - 1.25 + v* np.cos(u) , v*np.sin(u) ,  -1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1.075 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sl6 = Surface( lambda u, v: [ - 1.25 + v* np.cos(u) , v*np.sin(u) ,  -1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 1.075 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        ##### new floor and ceilings:
        snl1 = Surface( lambda u, v: [   v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl2 = Surface( lambda u, v: [   v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl3 = Surface( lambda u, v: [ - v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl4 = Surface( lambda u, v: [ - v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  -0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl5 = Surface( lambda u, v: [ - v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  -1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl6 = Surface( lambda u, v: [ - v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  -1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        ##### first tube:
        tub1 = Surface( lambda u, v: [ 1.25 +  (2.3 - np.sin(v)) *np.cos(u) /2 , (2.3 - np.sin(v)) *np.sin(u)/2 ,  np.cos(v) /2  + 0.6  ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tub2 = Surface( lambda u, v: [ 1.25 + (2.3 - 1.3*np.sin(v)) *np.cos(u) /2 , (2.3 - 1.3*np.sin(v)) *np.sin(u) /2 , 1.3*np.cos(v) /2 + 0.55   ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        ##### second tube:
        tub3 = Surface( lambda u, v: [ - 1.25 +  (2.3 - np.sin(v)) *np.cos(u) /2 , (2.3 - np.sin(v)) *np.sin(u)/2 ,  np.cos(v) /2  - 0.6  ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tub4 = Surface( lambda u, v: [ - 1.25 + (2.3 - 1.3*np.sin(v)) *np.cos(u) /2 , (2.3 - 1.3*np.sin(v)) *np.sin(u) /2 , 1.3*np.cos(v) /2 - 0.55   ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        g = VGroup(tub1, tub2, sl3, sl4, sl5, sl6, snl1, snl2, snl3, snl4, snl5, snl6)
        lids = VGroup(sl3, sl4, sl5, sl6)
        self.add(g)
        self.wait()
        self.play(FadeOut(lids))
        self.wait()
        self.play(Create(tub3), Create(tub4))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/5)
        self.wait(5)
        self.stop_ambient_camera_rotation()
        self.wait()


class lag7(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 70 * DEGREES, theta = 120 * DEGREES)
        snl1 = Surface( lambda u, v: [   v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl2 = Surface( lambda u, v: [   v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl3 = Surface( lambda u, v: [ - v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl4 = Surface( lambda u, v: [ - v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  -0.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl5 = Surface( lambda u, v: [ - v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  -1.1    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        snl6 = Surface( lambda u, v: [ - v*1.25 +  (v * (1.075) + (1-v) * 2.5)*np.cos(u)  ,  (v * (1.075) + (1-v) * 2.5)*np.sin(u)   ,  -1.2    ],
            u_range = [0, 2*PI],
            v_range = [0, 1 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        ##### first tube:
        tub1 = Surface( lambda u, v: [ 1.25 +  (2.3 - np.sin(v)) *np.cos(u) /2 , (2.3 - np.sin(v)) *np.sin(u)/2 ,  np.cos(v) /2  + 0.6  ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tub2 = Surface( lambda u, v: [ 1.25 + (2.3 - 1.3*np.sin(v)) *np.cos(u) /2 , (2.3 - 1.3*np.sin(v)) *np.sin(u) /2 , 1.3*np.cos(v) /2 + 0.55   ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        ##### second tube:
        tub3 = Surface( lambda u, v: [ - 1.25 +  (2.3 - np.sin(v)) *np.cos(u) /2 , (2.3 - np.sin(v)) *np.sin(u)/2 ,  np.cos(v) /2  - 0.6  ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tub4 = Surface( lambda u, v: [ - 1.25 + (2.3 - 1.3*np.sin(v)) *np.cos(u) /2 , (2.3 - 1.3*np.sin(v)) *np.sin(u) /2 , 1.3*np.cos(v) /2 - 0.55   ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        g = VGroup(tub1, tub2, tub3, tub4, snl1, snl2, snl3, snl4, snl5, snl6)
        self.add(g)
        self.wait()
        sc1 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) ,  0.6  + np.cos(v)/2    ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc2 = Surface( lambda u, v: [ (2.5 + np.sin(v)/2) *np.cos(u) , (2.5 + np.sin(v)/2) *np.sin(u) , - 0.6  + np.cos(v)/2     ],
            u_range = [0, 2*PI],
            v_range = [0, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc3 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,   0.6  + 0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [0, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc4 = Surface( lambda u, v: [ (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*sqrt(3) - 0.6*np.sin(v)) *np.sin(u) ,  0.6*np.cos(v)   ],
            u_range = [0, 2*PI],
            v_range = [PI/3, 2*PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sc5 = Surface( lambda u, v: [ (2.5 + 0.6*np.sin(v)) *np.cos(u) , (2.5 + 0.6*np.sin(v)) *np.sin(u) ,  -0.6  + 0.6*np.cos(v)  ],
            u_range = [0, 2*PI],
            v_range = [PI/3, PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sides = VGroup(sc1, sc2, sc3, sc4, sc5)
        self.play(Create(sides))
        self.wait()


