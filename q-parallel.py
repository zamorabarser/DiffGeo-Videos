from asyncio import threads
from this import d
from tkinter import E
from manim import *
from numpy import sqrt
import math



class ti(Scene):
    def construct(self):
        t1 = Text("Parallel transport and the", font_size=60).shift(0.5*UP)
        t2 = Text("Gauss-Bonnet Theorem", font_size=60).shift(0.5*DOWN)
        self.play(Write(t1), Write(t2))
        self.wait()        


class testscene(Scene):
    def construct(self):
        t1 = Tex(r"$q = 2 \pi $", font_size=60).shift(0.5*UP)
        self.play(Write(t1))
        self.wait()        
        self.play(t1.animate.shift(UP).set_color(BLUE))
        self.wait()


class par0(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/4)
        gaml = Tex(r"$\gamma $", r" $: [a,b] \to $ ", r"$\Sigma$", r", ", r"$ p = \gamma (a) $", r", ", r"$q = \gamma (b)$", font_size = 34).shift(2.4*UP)
        pgaml = Tex(r"$P_{\gamma} :$ ", r"$T_p \Sigma$", r" $\to $ ", r"$T_q \Sigma$", font_size = 34).next_to(gaml, DOWN).shift(0.2*DOWN)
        gli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(gaml[0], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE).next_to(gaml[2], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [1,0,0], color = BLUE).next_to(gaml[4], DOWN).shift(0.15*UP)
        qli = Line([0,0,0], [1,0,0], color = GREEN).next_to(gaml[6], DOWN).shift(0.15*UP)
        tpli = Line([0,0,0], [0.5,0,0], color = BLUE).next_to(pgaml[1], DOWN).shift(0.15*UP)
        tqli = Line([0,0,0], [0.5,0,0], color = GREEN).next_to(pgaml[3], DOWN).shift(0.15*UP)
        r = 2
        h = - 1.5
        op = 0.2
        p = Sphere(radius = 0.05)
        p.set_color(BLUE).shift( [ r* sqrt(3)/2 , 0, r/2 + h ] )
        pp = p.copy()
        q = Sphere(radius = 0.05)
        q.set_color(GREEN).shift([ - r* sqrt(3)/4 , r* 3/4, r/2 + h])
        g = ParametricFunction(lambda t : [  r * np.cos(t) * sqrt(3)/2  , r* np.sin(t) * sqrt(3)/2 , r/2 + h  ] , t_range = [0, 2*PI/3], color = ORANGE)
        s = Surface(
            lambda u, v:  [ r* np.sin(v)*np.cos(u),r* np.sin(v)*np.sin(u) , r*  np.cos(v) + h ],
            u_range = [ 0, TAU ],
            v_range = [ 0 , PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 32, 16]
        )
        tps = Surface(
            lambda u, v: [ r * sqrt(3)/2 - v/2 , u , r/2 + h + v * sqrt(3)/2  ] ,
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = op, 
            resolution=8
        )
        tqs = Surface(
            lambda u, v: [ - r* sqrt(3)/4  -3* sqrt(3) * u/8  - 5* v /8  , r* 3/4 + u/8 - 3* sqrt(3)* v/8  , r/2 + h + v*sqrt(3)/4 - 3* u / 4  ]  ,
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = op, 
            resolution=8
        )
        tpp = tps.copy()
        t = ValueTracker(0)
        tpp.add_updater( lambda x: x.become( 
            Surface(
                lambda u, v: [r*np.cos(t.get_value())*sqrt(3)/2 + u * (np.cos(t.get_value()/2)*(-np.sin(t.get_value())) - np.sin(t.get_value()/2)*(-np.cos(t.get_value())/2) ) + v * ( np.sin(t.get_value()/2)*(-np.sin(t.get_value())) + np.cos(t.get_value()/2)*(-np.cos(t.get_value())/2)  ) , r*np.sin(t.get_value())*sqrt(3)/2 + u * (np.cos(t.get_value()/2)*(np.cos(t.get_value())) - np.sin(t.get_value()/2)*(-np.sin(t.get_value())/2)) + v * (np.sin(t.get_value()/2)*(np.cos(t.get_value())) + np.cos(t.get_value()/2)*(-np.sin(t.get_value())/2) ) , r/2 + h + u * ( - np.sin(t.get_value()/2)*(sqrt(3)/2)) + v * ( np.cos(t.get_value()/2)*(sqrt(3)/2) ) ] ,
                u_range = [ -1, 1 ],
                v_range = [ -1, 1 ],
                checkerboard_colors = [BLUE, BLUE],
                fill_opacity = op, 
                resolution=8
            )   
        ))
        pp.add_updater( lambda x: x.become( Sphere(radius = 0.05).set_color(BLUE).shift( [  r * np.cos(t.get_value()) * sqrt(3)/2  , r* np.sin(t.get_value()) * sqrt(3)/2 , r/2 + h  ] )   ))
        self.add_fixed_in_frame_mobjects(gaml, pgaml, gli, pli, qli, sli, tpli, tqli)
        self.remove(gaml, pgaml, gli, pli, qli, sli, tpli, tqli)
        self.play(Create(gaml), Create(p), Create(q), Create(s), Create(g), Create(gli), Create(pli), Create(qli))
        self.wait()
        self.play(Create(tps), Create(tqs), Create(pgaml), Create(tpli), Create(tqli))
        self.wait()
        self.add(tpp)
        self.add(pp)
        self.play(t.animate.set_value(2*PI/3), rate_func=rate_functions.linear, run_time = 5)
        self.wait()
        self.remove(tpp, pp)
        self.play(FadeOut(tps), FadeOut(tqs))
        self.wait()


class par1(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/4)
        gaml = Tex(r"$\gamma $", r" $: [a,b] \to $ ", r"$\Sigma$", r", ", r"$ p = \gamma (a) $", r", ", r"$q = \gamma (b)$", font_size = 34).shift(2.4*UP)
        pgaml = Tex(r"$P_{\gamma} :$ ", r"$T_p \Sigma$", r" $\to $ ", r"$T_q \Sigma$", font_size = 34).next_to(gaml, DOWN).shift(0.2*DOWN)
        w0l = Tex(r"$W_0$", r" $\in T_p \Sigma$", font_size = 34).shift(2.4*UP + 3 *RIGHT)
        wl = Tex(r"$W$", r" $: [a,b] \to \mathbb{R}^3 $, $W(a) = W_0,$ $D_t(W) = 0 $", font_size = 34).next_to(w0l, DOWN).shift(0.2*DOWN)
        pgw = Tex(r"$P_{\gamma} (W_0) : = W (b)$", font_size = 34).shift(4*LEFT + 0.2*DOWN)
        boxpg = SurroundingRectangle(pgw, buff = .1, color = PINK)
        gli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(gaml[0], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE).next_to(gaml[2], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [1,0,0], color = BLUE).next_to(gaml[4], DOWN).shift(0.15*UP)
        qli = Line([0,0,0], [1,0,0], color = GREEN).next_to(gaml[6], DOWN).shift(0.15*UP)
        tpli = Line([0,0,0], [0.5,0,0], color = BLUE).next_to(pgaml[1], DOWN).shift(0.15*UP)
        tqli = Line([0,0,0], [0.5,0,0], color = GREEN).next_to(pgaml[3], DOWN).shift(0.15*UP)
        wli = Line([0,0,0], [0.4,0,0], color = PINK).next_to(w0l[0], DOWN).shift(0.15*UP)
        r = 2
        h = - 1.5
        op = 0.2
        l = 1.2
        p = Sphere(radius = 0.05)
        p.set_color(BLUE).shift( [ r* sqrt(3)/2 , 0, r/2 + h ] )
        q = Sphere(radius = 0.05)
        q.set_color(GREEN).shift([ - r* sqrt(3)/4 , r* 3/4, r/2 + h])
        w0 = Arrow([r*sqrt(3)/2 - 0.2 * sqrt(3) / 4 , -0.2 /2, r/2 + h + 0.2 * 3 /4 ], [r*sqrt(3)/2 + l * sqrt(3) / 4 , l * 1/2, r/2 + h - l * 3/4 ] , color = PINK)
        wco = []
        wtco = []
        w = []
        wvf = VGroup()
        for i in range(1,6):
            wco.append([])
            wtco.append([])
            wco[i-1].append( r * np.cos(2*PI *i / 15) * sqrt(3)/2 )
            wco[i-1].append( r* np.sin(2*PI *i / 15) * sqrt(3)/2 )
            wco[i-1].append( r/2 + h  )
            wtco[i-1].append(   np.cos(PI*i/15 + PI /3 ) * (- np.sin( 2*PI*i/15 ))  - np.sin(PI*i/15 + PI /3 ) * (-np.cos(2*PI* i /15)/2) ) 
            wtco[i-1].append(   np.cos(PI*i/15 + PI /3 ) * ( np.cos( 2*PI*i/15 ))  - np.sin(PI*i/15 + PI /3 ) * (-np.sin(2*PI* i /15)/2) )  
            wtco[i-1].append(   - np.sin(PI*i/15 + PI /3 ) * ( sqrt(3)/2) ) 
            w.append(Arrow( [wco[i-1][0] - 0.2 * wtco[i-1][0] , wco[i-1][1] - 0.2 * wtco[i-1][1] , wco[i-1][2] - 0.2 * wtco[i-1][2]], [wco[i-1][0] + l * wtco[i-1][0] , wco[i-1][1] +l * wtco[i-1][1] , wco[i-1][2] +l * wtco[i-1][2]] , color = PINK))
            wvf.add(w[i-1])
        g = ParametricFunction(lambda t : [  r * np.cos(t) * sqrt(3)/2  , r* np.sin(t) * sqrt(3)/2 , r/2 + h  ] , t_range = [0, 2*PI/3], color = ORANGE)
        s = Surface(
            lambda u, v:  [ r* np.sin(v)*np.cos(u),r* np.sin(v)*np.sin(u) , r*  np.cos(v) + h ],
            u_range = [ 0, TAU ],
            v_range = [ 0 , PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 32, 16]
        )
        self.add_fixed_in_frame_mobjects(gaml, pgaml, gli, pli, qli, sli, tpli, tqli, w0l, wl, wli, pgw, boxpg)
        self.add(s, p, q, g)
        self.remove(w0l, wl, wli, pgw, boxpg)
        b1 = VGroup(gaml, pgaml, gli, pli, qli, sli, tpli, tqli)
        self.play(b1.animate.shift(3*LEFT))
        self.play(Write(w0l), Create(wli), Create(w0))
        self.wait()
        self.play(Create(wl))
        self.wait()
        self.play(Create(wvf))
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/6)
        self.wait()
        self.stop_ambient_camera_rotation()
        self.play(Create(boxpg), Write(pgw))
        self.wait()


class par1new(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/4 + PI /6 )
        gaml = Tex(r"$\gamma $", r" $: [a,b] \to $ ", r"$\Sigma$", r", ", r"$ p = \gamma (a) $", r", ", r"$q = \gamma (b)$", font_size = 34).shift(2.4*UP)
        pgaml = Tex(r"$P_{\gamma} :$ ", r"$T_p \Sigma$", r" $\to $ ", r"$T_q \Sigma$", font_size = 34).next_to(gaml, DOWN).shift(0.2*DOWN)
        w0l = Tex(r"$W_0$", r" $\in T_p \Sigma$", font_size = 34).shift(2.4*UP + 3 *RIGHT)
        wl = Tex(r"$W$", r" $: [a,b] \to \mathbb{R}^3 $, $W(a) = W_0,$ $D_t(W) = 0 $", font_size = 34).next_to(w0l, DOWN).shift(0.2*DOWN)
        pgw = Tex(r"$P_{\gamma} (W_0) : = W (b)$", font_size = 34).shift(4*LEFT + 0.2*DOWN)
        boxpg = SurroundingRectangle(pgw, buff = .1, color = PINK)
        gli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(gaml[0], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE).next_to(gaml[2], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [1,0,0], color = BLUE).next_to(gaml[4], DOWN).shift(0.15*UP)
        qli = Line([0,0,0], [1,0,0], color = GREEN).next_to(gaml[6], DOWN).shift(0.15*UP)
        tpli = Line([0,0,0], [0.5,0,0], color = BLUE).next_to(pgaml[1], DOWN).shift(0.15*UP)
        tqli = Line([0,0,0], [0.5,0,0], color = GREEN).next_to(pgaml[3], DOWN).shift(0.15*UP)
        wli = Line([0,0,0], [0.4,0,0], color = PINK).next_to(w0l[0], DOWN).shift(0.15*UP)
        r = 2
        h = - 1.5
        op = 0.2
        l = 1.2
        p = Sphere(radius = 0.05)
        p.set_color(BLUE).shift( [ r* sqrt(3)/2 , 0, r/2 + h ] )
        q = Sphere(radius = 0.05)
        q.set_color(GREEN).shift([ - r* sqrt(3)/4 , r* 3/4, r/2 + h])
        w0 = Arrow([r*sqrt(3)/2 - 0.2 * sqrt(3) / 4 , -0.2 /2, r/2 + h + 0.2 * 3 /4 ], [r*sqrt(3)/2 + l * sqrt(3) / 4 , l * 1/2, r/2 + h - l * 3/4 ] , color = PINK)
        wco = []
        wtco = []
        w = []
        wvf = VGroup()
        for i in range(1,6):
            wco.append([])
            wtco.append([])
            wco[i-1].append( r * np.cos(2*PI *i / 15) * sqrt(3)/2 )
            wco[i-1].append( r* np.sin(2*PI *i / 15) * sqrt(3)/2 )
            wco[i-1].append( r/2 + h  )
            wtco[i-1].append(   np.cos(PI*i/15 + PI /3 ) * (- np.sin( 2*PI*i/15 ))  - np.sin(PI*i/15 + PI /3 ) * (-np.cos(2*PI* i /15)/2) ) 
            wtco[i-1].append(   np.cos(PI*i/15 + PI /3 ) * ( np.cos( 2*PI*i/15 ))  - np.sin(PI*i/15 + PI /3 ) * (-np.sin(2*PI* i /15)/2) )  
            wtco[i-1].append(   - np.sin(PI*i/15 + PI /3 ) * ( sqrt(3)/2) ) 
            w.append(Arrow( [wco[i-1][0] - 0.2 * wtco[i-1][0] , wco[i-1][1] - 0.2 * wtco[i-1][1] , wco[i-1][2] - 0.2 * wtco[i-1][2]], [wco[i-1][0] + l * wtco[i-1][0] , wco[i-1][1] +l * wtco[i-1][1] , wco[i-1][2] +l * wtco[i-1][2]] , color = PINK))
            wvf.add(w[i-1])
        g = ParametricFunction(lambda t : [  r * np.cos(t) * sqrt(3)/2  , r* np.sin(t) * sqrt(3)/2 , r/2 + h  ] , t_range = [0, 2*PI/3], color = ORANGE)
        s = Surface(
            lambda u, v:  [ r* np.sin(v)*np.cos(u),r* np.sin(v)*np.sin(u) , r*  np.cos(v) + h ],
            u_range = [ 0, TAU ],
            v_range = [ 0 , PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 32, 16]
        )
        self.add_fixed_in_frame_mobjects(gaml, pgaml, gli, pli, qli, sli, tpli, tqli, w0l, wl, wli, pgw, boxpg)
        self.add(s, p, q, g)
        self.remove(w0l, wl, wli, pgw, boxpg)
        b1 = VGroup(gaml, pgaml, gli, pli, qli, sli, tpli, tqli)
        self.play(b1.animate.shift(3*LEFT))
        self.play(Write(w0l), Create(wli), Create(w0))
        self.wait()
        self.play(Create(wl))
        self.wait()
        self.play(Create(wvf))
        self.wait(2)
        self.play(Create(boxpg), Write(pgw))
        self.wait()




class par2(Scene):
    def construct(self):
        gaml = Tex(r"$\gamma $", r" $: [a,b] \to $ ", r"$\Sigma$", r", ", r"$ p = \gamma (a) $", r", ", r"$q = \gamma (b)$", font_size = 34).shift(2.4*UP+ 3*LEFT)
        pgaml = Tex(r"$P_{\gamma} :$ ", r"$T_p \Sigma$", r" $\to $ ", r"$T_q \Sigma$", font_size = 34).next_to(gaml, DOWN).shift(0.2*DOWN)
        w0l = Tex(r"$W_0$", r" $\in T_p \Sigma$", font_size = 34).shift(2.4*UP + 3 *RIGHT)
        wl = Tex(r"$W$", r" $: [a,b] \to \mathbb{R}^3 $, $W(a) = W_0,$ $D_t(W) = 0 $", font_size = 34).next_to(w0l, DOWN).shift(0.2*DOWN)
        p = Tex(r"Proposition", font_size = 34, color = BLUE).shift(5.4*LEFT+ UP)
        l = p.get_left()[0]
        p1 = Tex(r"The vector field $W$ exists and it is unique.", font_size = 34).next_to(p, DOWN).shift(0.07*UP)
        pr = Tex(r"Proof", font_size = 34, color = BLUE).next_to(p1, DOWN)
        pr1 = Tex(r"If $s : U \to \Sigma$ covers $\gamma$, then $\gamma ( t) = s (u (t) , v (t) )$ and $W (t) = \alpha (t)s_u + \beta (t) s_v$.", font_size = 34).next_to(pr, DOWN).shift(0.07*UP)
        pr2 = Tex(r"$D_t W =  \alpha ^{\prime} s_u + \alpha D_t s_u + \beta ^{\prime} s_v + \beta D_t s_v $", font_size = 34).next_to(pr1, DOWN).shift(0.07*UP)
        pr3 = Tex(r"$D_tW = 0$ if and only if", font_size = 34).next_to(pr2, DOWN).shift(0.07*UP)
        pr4 = Tex(r"$ \alpha ^{\prime} + \alpha ( \Gamma_{11}^1 u^{\prime} + \Gamma _{12} ^1 v^{\prime} )  + \beta ( \Gamma_{12}^1 u^{\prime} + \Gamma _{22}^1 v^{\prime} ) = 0  $", font_size = 34).next_to(pr3, DOWN)
        pr5 = Tex(r"$  \beta ^{\prime} + \alpha ( \Gamma_{11}^2 u^{\prime} + \Gamma _{12} ^2 v^{\prime} )  + \beta ( \Gamma_{12}^2 u^{\prime} + \Gamma _{22}^2 v^{\prime} ) = 0  $", font_size = 34).next_to(pr4, DOWN).shift(0.07*UP)
        pge = VGroup(pr4, pr5)
        p1.shift((l - p1.get_left()[0])*RIGHT)
        pr.shift((l - pr.get_left()[0])*RIGHT)
        pr1.shift((l - pr1.get_left()[0])*RIGHT)
        pr3.shift((l - pr3.get_left()[0])*RIGHT)
        pr2.shift(pr2.get_center()[0]*LEFT)
        pr4.shift(pr4.get_center()[0]*LEFT)
        pr5.shift(pr5.get_center()[0]*LEFT)
        boxpg = SurroundingRectangle(pge, buff = .1, color = PINK)
        gli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(gaml[0], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE).next_to(gaml[2], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [1,0,0], color = BLUE).next_to(gaml[4], DOWN).shift(0.15*UP)
        qli = Line([0,0,0], [1,0,0], color = GREEN).next_to(gaml[6], DOWN).shift(0.15*UP)
        tpli = Line([0,0,0], [0.5,0,0], color = BLUE).next_to(pgaml[1], DOWN).shift(0.15*UP)
        tqli = Line([0,0,0], [0.5,0,0], color = GREEN).next_to(pgaml[3], DOWN).shift(0.15*UP)
        wli = Line([0,0,0], [0.4,0,0], color = PINK).next_to(w0l[0], DOWN).shift(0.15*UP)
        self.add(gaml, pgaml, w0l, wl, gli, pli, qli, sli, tpli, tqli, wli)
        self.wait()
        self.play(Write(p), Write(p1))
        self.wait()
        self.play(Write(pr), Write(pr1))
        self.wait()
        self.play(Write(pr2))
        self.wait()
        self.play(Write(pr3), Write(pr4), Write(pr5), Create(boxpg))
        self.wait()
        

class par3(Scene):
    def construct(self):
        gaml = Tex(r"$\gamma $", r" $: [a,b] \to $ ", r"$\Sigma$", r", ", r"$ p = \gamma (a) $", r", ", r"$q = \gamma (b)$", font_size = 34).shift(2.4*UP+ 3*LEFT)
        pgaml = Tex(r"$P_{\gamma} :$ ", r"$T_p \Sigma$", r" $\to $ ", r"$T_q \Sigma$", font_size = 34).next_to(gaml, DOWN).shift(0.2*DOWN)
        w0l = Tex(r"$W_0$", r" $\in T_p \Sigma$", font_size = 34).shift(2.4*UP + 3 *RIGHT)
        wl = Tex(r"$W$", r" $: [a,b] \to \mathbb{R}^3 $, $W(a) = W_0,$ $D_t(W) = 0 $", font_size = 34).next_to(w0l, DOWN).shift(0.2*DOWN)
        p = Tex(r"Proposition", font_size = 34, color = BLUE).shift(5.4*LEFT+ UP)
        l = p.get_left()[0]
        p1 = Tex(r"The vector field $W$ exists and it is unique.", font_size = 34).next_to(p, DOWN).shift(0.07*UP)
        pr = Tex(r"Definition", font_size = 34, color = BLUE).next_to(p1, DOWN).shift(0.17*DOWN)
        pr1 = Tex(r"The ", r"parallel transport along $\gamma$ ",r" is defined as $P_{\gamma} (W_0) = W(b)$." , font_size = 34).next_to(pr, DOWN)
        pr2 = Tex(r"Definition", font_size = 34, color = BLUE).next_to(pr1, DOWN).shift(0.1*DOWN)
        pr3 = Tex(r"For $\gamma :[a,b] \to \Sigma$ piece-wise smooth, a vector field $W : [a,b] \to \mathbb{R}^3$ along $\gamma$ " , font_size = 34).next_to(pr2, DOWN)
        pr4 = Tex(r"is called ", r"parallel", r" if $D_tW = 0 $ along the smooth pieces of $\gamma$" , font_size = 34).next_to(pr3, DOWN)
        p1.shift((l - p1.get_left()[0])*RIGHT)
        pr.shift((l - pr.get_left()[0])*RIGHT)
        pr1.shift((l - pr1.get_left()[0])*RIGHT)
        pr2.shift((l - pr2.get_left()[0])*RIGHT)
        pr3.shift((l - pr3.get_left()[0])*RIGHT)
        pr4.shift((l - pr4.get_left()[0])*RIGHT)
        gli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(gaml[0], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE).next_to(gaml[2], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [1,0,0], color = BLUE).next_to(gaml[4], DOWN).shift(0.15*UP)
        qli = Line([0,0,0], [1,0,0], color = GREEN).next_to(gaml[6], DOWN).shift(0.15*UP)
        tpli = Line([0,0,0], [0.5,0,0], color = BLUE).next_to(pgaml[1], DOWN).shift(0.15*UP)
        tqli = Line([0,0,0], [0.5,0,0], color = GREEN).next_to(pgaml[3], DOWN).shift(0.15*UP)
        wli = Line([0,0,0], [0.4,0,0], color = PINK).next_to(w0l[0], DOWN).shift(0.15*UP)
        ptli = Line([0,0,0], [3.6,0,0], color = PINK).next_to(pr1[1], DOWN).shift(0.15*UP)
        pfli = Line([0,0,0], [1,0,0], color = PINK).next_to(pr4[1], DOWN).shift(0.15*UP)
        self.add(gaml, pgaml, w0l, wl, gli, pli, qli, sli, tpli, tqli, wli, p, p1)
        self.wait()
        self.play(Write(pr), Write(pr1), p1.animate.shift(0.07*DOWN), Create(ptli))
        self.wait()
        self.play(Write(pr3), Write(pr2), Write(pr4), Create(pfli))
        self.wait()
        

class pais(Scene):
    def construct(self):
        p = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift(0.8*RIGHT+ 0.5*DOWN)
        l = p.get_left()[0]
        p0 = Tex(r"For a piece-wise smooth curve $\gamma : [a,b ] \to \Sigma$, $p = \gamma (a)$, $q = \gamma (b)$,", font_size = 34).next_to(p, DOWN)
        p1 = Tex(r"the parallel transport $P_{\gamma} : T_p \Sigma \to T_q \Sigma$ is an isometry.", font_size = 34).next_to(p0, DOWN)
        pr = Tex(r"Proof", font_size = 34, color = BLUE).next_to(p1, DOWN).shift(0.1*DOWN)
        pr1 = Tex(r"Take $W_0 \in T_p \Sigma$ and $W : [a,b] \to \mathbb{R}^3$ a parallel vector field along $\gamma$", font_size = 34).next_to(pr, DOWN)
        pr2 = Tex(r"with $W ( a) = W_0$. ", r"Then $\frac{d}{dt} \vert W \vert ^2  = 2$ $ W \cdot W^{\prime}$ ", r"$= 0$.", font_size = 34).next_to(pr1, DOWN)
        pr3 = Tex(r"Consequently, $\vert P_{\gamma} (W_0) \vert = \vert W (b) \vert  = \vert W_0 \vert$.", font_size = 34).next_to(pr2, DOWN)
        c = Tex(r"Corollary", font_size = 34, color = BLUE).next_to(pr3, DOWN).shift(0.1*DOWN)
        c1 = Tex(r"If $p = q$, then $P_{\gamma} : T_p \Sigma \to T_p \Sigma$ is a rotation.", font_size = 34).next_to(c, DOWN)
        p0.shift((l - p0.get_left()[0])*RIGHT)
        p1.shift((l - p1.get_left()[0])*RIGHT)
        pr.shift((l - pr.get_left()[0])*RIGHT)
        pr1.shift((l - pr1.get_left()[0])*RIGHT)
        pr2.shift((l - pr2.get_left()[0])*RIGHT)
        pr3.shift((l - pr3.get_left()[0])*RIGHT)
        c.shift((l - c.get_left()[0])*RIGHT)
        c1.shift((l - c1.get_left()[0])*RIGHT)
        self.play(Write(p), Write(p1), Write(p0))
        self.wait()
        self.play(Write(pr), Write(pr1), Write(pr2[0]))
        self.wait()
        self.play(Write(pr2[1]))        
        self.wait()
        self.play(Write(pr2[2]))        
        self.wait()
        self.play(Write(pr3))        
        self.wait()
        self.play(Write(c), Write(c1))        
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
        gk1 = Tex(r"$ k _{g} (t) : = \gamma ^{\prime \prime } (t) \cdot (N (\gamma (t) )\times \gamma ^{\prime}(t) ) $", font_size=34).next_to(gk, DOWN)
        gm1.shift(gm1.get_center()[0]*LEFT)
        gk.shift((l - gk.get_left()[0])*RIGHT)
        gk1.shift(gk1.get_center()[0]*LEFT)
        gkb = SurroundingRectangle(gk1, buff = .1, color = BLUE)
        kneg = Tex(r"$k_{g} < 0$", font_size=34).next_to(gk1, DOWN).shift(4.5*LEFT + 0.5*DOWN)
        kzer = Tex(r"$k_{g} = 0$", font_size=34).next_to(gk1, DOWN).shift(0.5*DOWN)
        kpos = Tex(r"$k_{g} > 0$", font_size=34).next_to(gk1, DOWN).shift(4.5*RIGHT+ 0.5*DOWN)
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
        all = VGroup(gm, gm1, nli, gk, gk1, gkb, s1, s2, s3, n1, n2, n3, kneg, g1, g1t1, g1t2, knb, kzer, g2, g2t1, g2t2, kzb, kpos, g3, g3t1, g3t2, kpb )
        self.play(Create(all))
        self.wait()


class tgk(Scene):
    def construct(self):
        ti = Tex(r"Total signed curvature", font_size=34, color=BLUE).to_edge(UL).shift(0.7*DOWN + 0.5*RIGHT)
        dtc = Tex(r"For  $\gamma : [a,b] \to \Sigma$ smooth unit-speed, the total geodesic curvature is", font_size=34).next_to(ti,DOWN)
        l = ti.get_left()[0]
        dtc.shift((dtc.get_left()[0]-l)*LEFT)
        tcf = Tex(r"$ \Psi (\gamma ) : = \int_a^b k_g (t) dt $.", font_size =34).next_to(dtc, DOWN)
        tcf.shift(tcf.get_center()[0]*LEFT)
        pw = Tex(r"If  $\gamma$ is only piecewise smooth, there are $a = t_0 < \ldots < t_n = b$", font_size=34).next_to(tcf,DOWN)
        pw.shift((pw.get_left()[0]-l)*LEFT)
        pw2 = Tex(r"such that $\gamma \vert _{[t_{i-1}, t_i]}$ is smooth unit-speed for each $i$.", font_size=34).next_to(pw,DOWN)
        pw2.shift((pw2.get_left()[0]-l)*LEFT)
        pw3 = Tex(r"$\Psi (\gamma )  : =  \sum_{i=1}^n  \int_{t_{i-1}}^{t_i} k_g (t)   dt + \sum _{i=1}^{n-1} \angle ( \gamma^{\prime} (t_i^-), \gamma^{\prime}(t_i^+)  )$.", font_size=34).next_to(pw2,DOWN)
        pw3.shift((pw2.get_center()[0])*LEFT)
        pw33 = Tex(r"where $\angle ( \cdot , \cdot ) \in (- \pi , \pi) $ is a signed angle.", font_size=34).next_to(pw3,DOWN)
        pw33.shift((pw33.get_left()[0]-l)*LEFT)
        pw4 = Tex(r"And if $\gamma$ is closed, ", font_size=34).next_to(pw33,DOWN)
        pw4.shift((pw4.get_left()[0]-l)*LEFT)
        pw5 = Tex(r"$\Psi (\gamma )  : =  \sum_{i=1}^n  \int_{t_{i-1}}^{t_i} k _g (t)  dt + \sum _{i=1}^{n-1} \angle ( \gamma^{\prime} (t_i^-), \gamma^{\prime}(t_i^+)  ) + \angle (  \gamma^{\prime} (b^{-}), \gamma ^{\prime}(a^{+}) )$.", font_size=34).next_to(pw4,DOWN)
        pw5.shift((pw5.get_center()[0])*LEFT)
        self.play(Write(ti), Write(dtc), Write(tcf))
        self.wait()
        self.play(Write(pw), Write(pw2), Write(pw3), Write(pw33))
        self.wait()
        self.play(Write(pw4), Write(pw5))
        self.wait()


class pttk(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/4)
        pro = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift(0.5*RIGHT+ 0.7*DOWN)
        l = pro.get_left()[0]
        pro1 = Tex(r"For a closed piece-wise smooth curve $\gamma : [a,b ] \to \Sigma$, $p = \gamma (a) = \gamma (b)$,", font_size = 34).next_to(pro, DOWN)
        pro2 = Tex(r"the parallel transport $P_{\gamma} : T_p \Sigma \to T_p \Sigma$ is a clockwise rotation with angle $\Psi (\gamma)$.", font_size = 34).next_to(pro1, DOWN)
        pro1.shift((l - pro1.get_left()[0])*RIGHT)
        pro2.shift((l - pro2.get_left()[0])*RIGHT)
        r = 2
        len = 1
        h = - 2
        op = 0.2
        p0 = [1/sqrt(2), 0, 1/sqrt(2)]
        p1 = [0, 1/sqrt(2), 1/sqrt(2)]
        p2 = [-1/sqrt(2), 0, 1/sqrt(2)]
        p3 = [0, -1/sqrt(2), 1/sqrt(2)]
        q1 = [-1/sqrt(6), sqrt(2/3), 1/sqrt(6)]
        q2 = [-sqrt(2/3), -1/sqrt(6), 1/sqrt(6)]
        q3 = [1/sqrt(6), -sqrt(2/3), 1/sqrt(6)]
        q4 = [sqrt(2/3), 1/sqrt(6), 1/sqrt(6)]
        def gameq(t):
            dummy = []
            if t < PI/3:
                dummy.append(np.cos(t)*p0[0] + np.sin(t)*q1[0])
                dummy.append(np.cos(t)*p0[1] + np.sin(t)*q1[1])
                dummy.append(np.cos(t)*p0[2] + np.sin(t)*q1[2])
                return dummy
            if t < 2*PI/3:
                dummy.append(np.cos(t- PI/3)*p1[0] + np.sin(t- PI/3)*q2[0])
                dummy.append(np.cos(t- PI/3)*p1[1] + np.sin(t- PI/3)*q2[1])
                dummy.append(np.cos(t- PI/3)*p1[2] + np.sin(t- PI/3)*q2[2])
                return dummy
            if t < PI:
                dummy.append(np.cos(t- 2*PI/3)*p2[0] + np.sin(t- 2*PI/3)*q3[0])
                dummy.append(np.cos(t- 2*PI/3)*p2[1] + np.sin(t- 2*PI/3)*q3[1])
                dummy.append(np.cos(t- 2*PI/3)*p2[2] + np.sin(t- 2*PI/3)*q3[2])
                return dummy
            dummy.append(np.cos(t- PI)*p3[0] + np.sin(t- PI)*q4[0])
            dummy.append(np.cos(t- PI)*p3[1] + np.sin(t- PI)*q4[1])
            dummy.append(np.cos(t- PI)*p3[2] + np.sin(t- PI)*q4[2])
            return dummy
        def gp(t):
            dummy = []
            if t < PI/3:
                dummy.append(-np.sin(t)*p0[0] + np.cos(t)*q1[0])
                dummy.append(-np.sin(t)*p0[1] + np.cos(t)*q1[1])
                dummy.append(-np.sin(t)*p0[2] + np.cos(t)*q1[2])
                return dummy
            if t < 2*PI/3:
                dummy.append(-np.sin(t- PI/3)*p1[0] + np.cos(t- PI/3)*q2[0])
                dummy.append(-np.sin(t- PI/3)*p1[1] + np.cos(t- PI/3)*q2[1])
                dummy.append(-np.sin(t- PI/3)*p1[2] + np.cos(t- PI/3)*q2[2])
                return dummy
            if t < PI:
                dummy.append(-np.sin(t- 2*PI/3)*p2[0] + np.cos(t- 2*PI/3)*q3[0])
                dummy.append(-np.sin(t- 2*PI/3)*p2[1] + np.cos(t- 2*PI/3)*q3[1])
                dummy.append(-np.sin(t- 2*PI/3)*p2[2] + np.cos(t- 2*PI/3)*q3[2])
                return dummy
            dummy.append(-np.sin(t- PI)*p3[0] + np.cos(t- PI)*q4[0])
            dummy.append(-np.sin(t- PI)*p3[1] + np.cos(t- PI)*q4[1])
            dummy.append(-np.sin(t- PI)*p3[2] + np.cos(t- PI)*q4[2])
            return dummy
        def nu(t):
            if t < PI/3:
                return [-1/sqrt(3), -1/sqrt(3), 1/sqrt(3)] 
            if t < 2*PI/3:
                return [ 1/sqrt(3), -1/sqrt(3), 1/sqrt(3)]
            if t < PI:
                return [1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]
            return [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]
        def wt(t):
            dummy = []
            if t < PI/3:
                dummy.append(gp(t)[0] - nu(t)[0])
                dummy.append(gp(t)[1] - nu(t)[1])
                dummy.append(gp(t)[2] - nu(t)[2])
                return dummy
            if t < 2*PI/3:
                dummy.append( (1 - sqrt(8))*gp(t)[0]/3 - (sqrt(8) + 1)*nu(t)[0]/3)
                dummy.append( (1 - sqrt(8))*gp(t)[1]/3 - (sqrt(8) + 1)*nu(t)[1]/3)
                dummy.append( (1 - sqrt(8))*gp(t)[2]/3 - (sqrt(8) + 1)*nu(t)[2]/3)
                return dummy
            if t < PI:
                dummy.append( - (sqrt(32) + 7)*gp(t)[0]/9 + (7 - sqrt(32))*nu(t)[0]/9)
                dummy.append( - (sqrt(32) + 7)*gp(t)[1]/9 + (7 - sqrt(32))*nu(t)[1]/9)
                dummy.append( - (sqrt(32) + 7)*gp(t)[2]/9 + (7 - sqrt(32))*nu(t)[2]/9)
                return dummy
            dummy.append(  (10*sqrt(2) - 23)*gp(t)[0]/27 + (10*sqrt(2) + 23)*nu(t)[0]/27)
            dummy.append(  (10*sqrt(2) - 23)*gp(t)[1]/27 + (10*sqrt(2) + 23)*nu(t)[1]/27)
            dummy.append(  (10*sqrt(2) - 23)*gp(t)[2]/27 + (10*sqrt(2) + 23)*nu(t)[2]/27)
            return dummy
        w = Arrow( [ r* gameq(0)[0] - 0.1 * wt(0)[0] , r*gameq(0)[1] - 0.1 * wt(0)[1], r*gameq(0)[2] - 0.1 * wt(0)[2] + h ], [ r* gameq(0)[0] + len * wt(0)[0] , r*gameq(0)[1] + len * wt(0)[1], r*gameq(0)[2] + len * wt(0)[2] + h ] , color = PINK)
        w2 = w.copy()
        p = Sphere(radius = 0.05)
        p.set_color(BLUE).shift( [ r /sqrt(2) , 0, r/sqrt(2) + h ] )
        g = ParametricFunction(lambda t : [ r * gameq(t)[0] , r * gameq(t)[1] , r * gameq(t)[2] +h ], t_range = [0, 4*PI/3], color = ORANGE)
        s = Surface(
            lambda u, v:  [ r* np.sin(v)*np.cos(u),r* np.sin(v)*np.sin(u) , r*  np.cos(v) + h ],
            u_range = [ 0, TAU ],
            v_range = [ 0 , PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 32, 16]
        )
        t = ValueTracker(0)
        w.add_updater( lambda x: x.become( 
            Arrow(
                [ r* gameq(t.get_value())[0] - 0.1 * wt(t.get_value())[0] , r*gameq(t.get_value())[1] - 0.1 * wt(t.get_value())[1], r*gameq(t.get_value())[2] - 0.1 * wt(t.get_value())[2] + h ],
                [ r* gameq(t.get_value())[0] + len * wt(t.get_value())[0] , r*gameq(t.get_value())[1] + len * wt(t.get_value())[1], r*gameq(t.get_value())[2] + len * wt(t.get_value())[2] + h ],
                color = PINK
            )   
        ))
        self.add_fixed_in_frame_mobjects(pro, pro1, pro2)
        self.remove(pro, pro1, pro2)
        self.play(Create(pro), Create(pro1), Create(pro2), Create(s), Create(p), Create(g))
        self.wait()
        self.play(Create(w))
        self.wait()
        self.add(w2)
        self.play(t.animate.set_value(PI/3), rate_func=rate_functions.linear, run_time = 3)
        #leng = sqrt(2)*len
        #tan = Arrow([r*p1[0], r*p1[1], r*p1[2] + h], [r*p1[0]+leng*q2[0], r*p1[1]+leng*q2[1],r*p1[2]+leng*q2[2]+ h], color = GREEN )
        #nor = Arrow([r*p1[0], r*p1[1], r*p1[2]+ h ], [r*p1[0]+leng*nu(PI/3)[0],r*p1[1]+leng*nu(PI/3)[1], r*p1[2]+leng*nu(PI/3)[2]+h ], color = GREEN)
        #newstart = Arrow([r*p1[0], r*p1[1] , r*p1[2]+h] , [ r*p1[0] - len*(sqrt(35)+1)*q2[0] - len*(sqrt(35)-1)*nu(PI/3)[0]  ,r*p1[1]- len*(sqrt(35)+1)*q2[1] - len*(sqrt(35)-1)*nu(PI/3)[1], r*p1[2]- len*(sqrt(35)+1)*q2[2] - len*(sqrt(35)-1)*nu(PI/3)[2] +h    ] , color = RED )
        #self.play(Create(newstart))
        #self.wait()
        self.play(t.animate.set_value(2*PI/3), rate_func=rate_functions.linear, run_time = 2)
        self.play(t.animate.set_value(PI), rate_func=rate_functions.linear, run_time = 2)
        self.play(t.animate.set_value(4*PI/3), rate_func=rate_functions.linear, run_time = 2)
        self.wait()



class pttk1(Scene):
    def construct(self):
        pro = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift(0.5*RIGHT+ 0.7*DOWN)
        l = pro.get_left()[0]
        pro1 = Tex(r"For a closed piece-wise smooth curve $\gamma : [a,b ] \to \Sigma$, $p = \gamma (a) = \gamma (b)$,", font_size = 34).next_to(pro, DOWN)
        pro2 = Tex(r"the parallel transport $P_{\gamma} : T_p \Sigma \to T_p \Sigma$ is a clockwise rotation with angle $\Psi (\gamma)$.", font_size = 34).next_to(pro1, DOWN)
        prf = Tex(r"Proof", font_size = 34, color = BLUE).next_to(pro2, DOWN)
        prf1 = Tex(r"Assume ", r"$\gamma$", r" is a concatenation of geodesics. Let ", r"$X$", r" $ : [a,b] \to \mathbb{R}^3$", font_size = 34).next_to(prf, DOWN)
        prf2 = Tex(r"be a parallel field along $\gamma$ with $X(a) = \gamma ^{\prime }(a). $", font_size = 34).next_to(prf1, DOWN)
        prf3 = Tex(r"$\angle ( $ ", r"$ \gamma ^{\prime} (t) $", r" $,$ ", r"$ X(t)$", r" $ )$ is constant along smooth pieces,", font_size = 34).next_to(prf2, DOWN)
        prf4 = Tex(r"and changes by $-\angle ( \gamma^{\prime} (t_0^{-}), \gamma ^{\prime} (t_0^{+}) )$ at the corners.", font_size = 34).next_to(prf3, DOWN)
        prf5 = Tex(r"$\Rightarrow $ $\angle ( X(a) , X(b) ) = \angle (\gamma^{\prime}(b) , X(b) ) = -$ ", r"$\sum$", r" $ \angle ( \gamma^{\prime} (t_0^{-}), \gamma ^{\prime} (t_0^{+}) ) = - \Psi (\gamma)$.", font_size = 34).next_to(prf4, DOWN)
        pro1.shift((l - pro1.get_left()[0])*RIGHT)
        pro2.shift((l - pro2.get_left()[0])*RIGHT)
        prf.shift((l - prf.get_left()[0])*RIGHT)
        prf1.shift((l - prf1.get_left()[0])*RIGHT)
        prf2.shift((l - prf2.get_left()[0])*RIGHT)
        prf3.shift((l - prf3.get_left()[0])*RIGHT)
        prf4.shift((l - prf4.get_left()[0])*RIGHT)
        prf5.shift((l - prf5.get_left()[0])*RIGHT)
        prf6 = Tex(r"$t_0$", font_size = 18).next_to(prf5[1], DOWN).shift(0.2*UP)
        gli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(prf1[1], DOWN).shift(0.15*UP)
        xli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(prf1[3], DOWN).shift(0.15*UP)
        xli2 = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(prf3[3], DOWN).shift(0.15*UP)
        gpli = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(prf3[1], DOWN).shift(0.15*UP)
        g0 = Line([0,-2,0], [3,-2,0], color = ORANGE)
        g1 = Line([3,-2,0], [4.5,-0.5, 0], color = ORANGE)
        x = Arrow([-0.1*sqrt(2),-2+ 0.1*sqrt(2) ,0], [1, -3 ,0], color = YELLOW)
        gp = Arrow([-0.2,-2,0], [ sqrt(2), -2 ,0], color = BLUE)
        self.add(pro, pro1, pro2)
        self.play(Create(prf), Create(prf1), Create(prf2), Create(xli), Create(gli))
        self.wait()
        self.play(Create(prf3), Create(prf4), Create(x), Create(gp), Create(g0), Create(g1), Create(xli2), Create(gpli))
        self.wait()
        self.play(x.animate.shift(3*RIGHT), gp.animate.shift(3*RIGHT), run_time = 2)
        self.play(Rotate(gp, angle = PI/4,  about_point = [3,-2,0]), run_time = 2)
        self.play(x.animate.shift(1.5*RIGHT + 1.5*UP), gp.animate.shift(1.5*RIGHT+1.5*UP), run_time = 2)
        self.wait()
        self.play(FadeOut(gp), FadeOut(x), FadeOut(g0), FadeOut(g1), Write(prf5), Write(prf6))
        self.wait()
        self.play(FadeOut(prf1), FadeOut(prf2), FadeOut(prf3), FadeOut(prf4), FadeOut(prf5), FadeOut(prf6), FadeOut(xli), FadeOut(xli2), FadeOut(gpli), FadeOut(gli))
        self.wait()








class pttk2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/4)
        pro = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift(0.5*RIGHT+ 0.7*DOWN)
        l = pro.get_left()[0]
        pro1 = Tex(r"For a closed piece-wise smooth curve $\gamma : [a,b ] \to \Sigma$, $p = \gamma (a) = \gamma (b)$,", font_size = 34).next_to(pro, DOWN)
        pro2 = Tex(r"the parallel transport $P_{\gamma} : T_p \Sigma \to T_p \Sigma$ is a clockwise rotation with angle $\Psi (\gamma)$.", font_size = 34).next_to(pro1, DOWN)
        prf = Tex(r"Proof", font_size = 34, color = BLUE).next_to(pro2, DOWN)
        prf1 = Tex(r"Assume ", r"$\gamma$", r" is smooth and let ", r"$X$", r" $ : [a,b] \to \mathbb{R}^3$", font_size = 34).next_to(prf, DOWN)
        prf2 = Tex(r"be a parallel field along $\gamma$ with $X(a) = \gamma ^{\prime }(a). $", font_size = 34).next_to(prf1, DOWN)
        prf3 = Tex(r"Also let ", r"$Y $", r" $= N \times X$.", r" Then", font_size = 34).next_to(prf2, DOWN)
        prf4 = Tex(r"$\gamma ^{\prime} (t)$", r" $ = \cos \theta (t) X (t) + \sin \theta (t) Y(t)$", font_size = 34).next_to(prf3, DOWN)
        prf5 = Tex(r"$\Rightarrow $ $\angle ( X(a) , X(b) ) = \angle (\gamma^{\prime}(b) , X(b) ) = -$ ", r"$\sum$", r" $ \angle ( \gamma^{\prime} (t_0^{-}), \gamma ^{\prime} (t_0^{+}) ) = - \Psi (\gamma)$.", font_size = 34).next_to(prf4, DOWN)
        pro1.shift((l - pro1.get_left()[0])*RIGHT)
        pro2.shift((l - pro2.get_left()[0])*RIGHT)
        prf.shift((l - prf.get_left()[0])*RIGHT)
        prf1.shift((l - prf1.get_left()[0])*RIGHT)
        prf2.shift((l - prf2.get_left()[0])*RIGHT)
        prf3.shift((l - prf3.get_left()[0])*RIGHT)
        prf4.shift((l - prf4.get_left()[0])*RIGHT)
        prf5.shift((l - prf5.get_left()[0])*RIGHT)
        gli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(prf1[1], DOWN).shift(0.15*UP)
        xli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(prf1[3], DOWN).shift(0.15*UP)
        yli = Line([0,0,0], [0.4,0,0], color = GREEN).next_to(prf3[1], DOWN).shift(0.15*UP)
        gpli = Line([0,0,0], [0.7,0,0], color = BLUE).next_to(prf4[0], DOWN).shift(0.15*UP)
        ht =2.3
        h = - 3.2
        r = 2.5
        len = 0.5
        op = 0.2
        def gam(t):
            return [r* np.cos(t)/2 - ht , r*np.sin(t)/2 + ht , r*np.sqrt(3)/2 + h]
        def gp(t):
            return [ - r* np.sin(t) , r*np.cos(t) ,  0 ]
        def nu(t):
            return [-sqrt(3)*r*np.cos(t)/2, -sqrt(3)*r*np.sin(t)/2 , r/2 ]
        def xt(t):
            return [ np.cos(sqrt(3)*t /2) * gp(t)[0] - np.sin(sqrt(3)*t/2)* nu(t)[0] , np.cos(sqrt(3)*t /2) * gp(t)[1] - np.sin(sqrt(3)*t/2)* nu(t)[1],  np.cos(sqrt(3)*t /2) * gp(t)[2] - np.sin(sqrt(3)*t/2)* nu(t)[2]]
        def yt(t):
            return [ np.sin(sqrt(3)*t /2) * gp(t)[0] + np.cos(sqrt(3)*t/2)* nu(t)[0] , np.sin(sqrt(3)*t /2) * gp(t)[1] + np.cos(sqrt(3)*t/2)* nu(t)[1],  np.sin(sqrt(3)*t /2) * gp(t)[2] + np.cos(sqrt(3)*t/2)* nu(t)[2]]
        x = []
        y = []
        gpp = []
        xvf = VGroup()
        yvf = VGroup()
        gpvf = VGroup()
        u = 8 
        for i in range(0,u + 1):
            x.append(  Arrow( [ gam(2*PI*i/u)[0] - 0.1* xt(2*PI*i/u)[0], gam(2*PI*i/u)[1] - 0.1* xt(2*PI*i/u)[1], gam(2*PI*i/u)[2] - 0.1* xt(2*PI*i/u)[2] ], [ gam(2*PI*i/u)[0] + len* xt(2*PI*i/u)[0], gam(2*PI*i/u)[1] + len* xt(2*PI*i/u)[1], gam(2*PI*i/u)[2] + len* xt(2*PI*i/u)[2] ] , color = YELLOW)  )
            y.append(  Arrow( [ gam(2*PI*i/u)[0] - 0.1* yt(2*PI*i/u)[0], gam(2*PI*i/u)[1] - 0.1* yt(2*PI*i/u)[1], gam(2*PI*i/u)[2] - 0.1* yt(2*PI*i/u)[2] ], [ gam(2*PI*i/u)[0] + len* yt(2*PI*i/u)[0], gam(2*PI*i/u)[1] + len* yt(2*PI*i/u)[1], gam(2*PI*i/u)[2] + len* yt(2*PI*i/u)[2] ] , color = GREEN) )
            gpp.append(  Arrow( [ gam(2*PI*i/u)[0] - 0.1* gp(2*PI*i/u)[0], gam(2*PI*i/u)[1] - 0.1* gp(2*PI*i/u)[1], gam(2*PI*i/u)[2] - 0.1* gp(2*PI*i/u)[2] ], [ gam(2*PI*i/u)[0] + len* gp(2*PI*i/u)[0], gam(2*PI*i/u)[1] + len* gp(2*PI*i/u)[1], gam(2*PI*i/u)[2] + len* gp(2*PI*i/u)[2] ] , color = BLUE)  )
            xvf.add(x[i])
            yvf.add(y[i])
            gpvf.add(gpp[i])
        p = Sphere(radius = 0.05)
        p.set_color(BLUE).shift( [  r / 2 - ht, ht , r*sqrt(3)/2 + h ] )
        g = ParametricFunction(lambda t : gam(t) , t_range = [0, 2*PI], color = ORANGE)
        s = Surface(
            lambda u, v:  [ r* np.sin(v)*np.cos(u)  - ht ,r* np.sin(v)*np.sin(u) + ht , r*  np.cos(v) + h ],
            u_range = [ 0, TAU ],
            v_range = [ 0 , PI/3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 32, 16]
        )
        self.add_fixed_in_frame_mobjects(pro, pro1, pro2, prf, prf1, prf2, prf3, prf4, prf5, gli, xli, yli, gpli)
        self.remove(prf1, prf1, prf2, prf3, prf4, prf5, gli, xli, yli, gpli)
        self.play(Create(prf1), Create(prf2), Create(xli), Create(gli), Create(s), Create(p), Create(g))
        self.wait()
        self.play(Create(xvf), run_time = 2)
        self.wait()
        self.play(Create(yvf), Create(prf3[0]), Create(prf3[1]), Create(prf3[2]), Create(yli))
        self.wait()
        self.play(Create(prf3[3]), Create(prf4), Create(gpli), Create(gpvf))
        self.wait()




class pttk3(Scene):
    def construct(self):
        h0 = 0
        pro = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift(0.5*RIGHT+ 0.7*DOWN)
        l = pro.get_left()[0]
        pro1 = Tex(r"For a closed piece-wise smooth curve $\gamma : [a,b ] \to \Sigma$, $p = \gamma (a) = \gamma (b)$,", font_size = 34).next_to(pro, DOWN)
        pro2 = Tex(r"the parallel transport $P_{\gamma} : T_p \Sigma \to T_p \Sigma$ is a clockwise rotation with angle $\Psi (\gamma)$.", font_size = 34).next_to(pro1, DOWN)
        prf = Tex(r"Proof", font_size = 34, color = BLUE).next_to(pro2, DOWN)
        prf1 = Tex(r"Assume ", r"$\gamma$", r" is smooth and let ", r"$X$", r" $ : [a,b] \to \mathbb{R}^3$", font_size = 34).next_to(prf, DOWN)
        prf2 = Tex(r"be a parallel field along $\gamma$ with $X(a) = \gamma ^{\prime }(a). $", font_size = 34).next_to(prf1, DOWN)
        prf3 = Tex(r"Also let ", r"$Y $", r" $= N \times X$.", r" Then", font_size = 34).next_to(prf2, DOWN)
        prf4 = Tex(r"$\gamma ^{\prime} (t)$", r" $ = \cos \theta (t) X (t) + \sin \theta (t) Y(t)$", font_size = 34).next_to(prf3, DOWN)
        prf5 = Tex(r"Claim: $k_g(t) = \theta ^{\prime}(t)$.", font_size = 34).next_to(prf4, DOWN)
        prf6 = Tex(r"$N \times \gamma ^{\prime} = \cos \theta (N \times X) + \sin \theta (N \times Y)$", font_size = 34).next_to(prf2, DOWN).shift(h0*UP)
        prf7 = Tex(r"$= \cos \theta  Y - \sin \theta X $", font_size = 34).next_to(prf6, DOWN)
        prf8 = Tex(r"$\gamma^{\prime \prime} = \cos \theta X ^{\prime} + \sin \theta Y ^{\prime} + ( \cos \theta Y - \sin \theta X  ) \theta ^{\prime}$", font_size = 34).next_to(prf7, DOWN)
        prf9 = Tex(r"$\Rightarrow $ $k_g = \gamma^{\prime \prime} \cdot (N \times \gamma ^{\prime } ) = (\cos ^2 \theta + \sin ^2 \theta ) \theta ^{\prime} = \theta^{\prime}$", font_size = 34).next_to(prf8, DOWN)
        prf10 = Tex(r"$\Rightarrow $ $\Psi (\gamma) = \int_a^b k_g  $ d$ t = \int_a^b \theta^{\prime } $ d$ t  = \theta (b)$.", font_size = 34).next_to(prf2, DOWN).shift(h0*UP)
        prf11 = Tex(r"$P_{\gamma}( $ ", r"$\gamma^{\prime}(a) $", r" $) = P_{\gamma}(X(a)) = $ ", r"$X(b)$", font_size = 34).next_to(prf10, DOWN)
        h = prf1.get_center()[1] + h0 - prf4.get_center()[1]
        pro1.shift((l - pro1.get_left()[0])*RIGHT)
        pro2.shift((l - pro2.get_left()[0])*RIGHT)
        prf.shift((l - prf.get_left()[0])*RIGHT)
        prf1.shift((l - prf1.get_left()[0])*RIGHT)
        prf2.shift((l - prf2.get_left()[0])*RIGHT)
        prf3.shift((l - prf3.get_left()[0])*RIGHT)
        prf4.shift((l - prf4.get_left()[0])*RIGHT)
        prf5.shift((l - prf5.get_left()[0])*RIGHT)
        prf6.shift((l - prf6.get_left()[0])*RIGHT)
        prf7.shift((l - prf7.get_left()[0] + 1.15 )*RIGHT)
        prf8.shift((l - prf8.get_left()[0])*RIGHT)
        prf9.shift((l - prf9.get_left()[0])*RIGHT)
        prf10.shift((l - prf10.get_left()[0])*RIGHT)
        prf11.shift((l - prf11.get_left()[0])*RIGHT)
        gli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(prf1[1], DOWN).shift(0.15*UP)
        xli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(prf1[3], DOWN).shift(0.15*UP)
        xli2 = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(prf11[3], DOWN).shift(0.15*UP)
        yli = Line([0,0,0], [0.4,0,0], color = GREEN).next_to(prf3[1], DOWN).shift(0.15*UP)
        gpli = Line([0,0,0], [0.7,0,0], color = BLUE).next_to(prf4[0], DOWN).shift(0.15*UP)
        gpli2 = Line([0,0,0], [0.7,0,0], color = BLUE).next_to(prf11[1], DOWN).shift(0.15*UP)
        x = Arrow([2.5 - 0.2,-2,0],[5.2,-2,0], color = YELLOW)
        y = Arrow([2.5,-2 - 0.2,0], [2.5,0.7,0], color = GREEN)
        gp = Arrow([2.5 -0.2 / sqrt(2),-2-0.2 / sqrt(2),0], [2.5+2.7/sqrt(2),-2+2.7/sqrt(2),0], color = BLUE)
        th = Angle(x, gp, radius=1, other_angle=False)
        thl = Tex(r"$\theta$", font_size = 34).next_to(th, RIGHT).shift(0.1*UP)
        self.add(pro, pro1, pro2, prf, prf1, prf2, prf3, prf4, gli, xli, yli, gpli)
        self.wait()
        self.play(Write(prf5), Create(x), Create(y), Create(gp), Create(th), Create(thl))
        self.wait()
        all = VGroup(pro, pro1, pro2, prf, gp, x, y, th, thl)
        self.play(all.animate.shift(h0*UP), FadeOut(prf1), FadeOut(prf2), FadeOut(prf3), prf4.animate.shift(h*UP), prf5.animate.shift(h*UP),  gpli.animate.shift(h*UP), FadeOut(xli), FadeOut(yli), FadeOut(gli))
        self.play(Write(prf6), Write(prf7))
        self.wait()
        self.play(Write(prf8))
        self.wait()
        self.play(Write(prf9))
        self.wait()
        self.play(FadeOut(prf6), FadeOut(prf7), FadeOut(prf8), FadeOut(prf9))
        self.wait()
        self.play(Write(prf10))
        self.wait()
        self.play(Write(prf11), Create(xli2), Create(gpli2))
        self.wait()




class gbf(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=78 * DEGREES, theta=  PI/4)
        t = Tex(r"Theorem (Gauss--Bonnet Formula)", font_size = 34, color = BLUE).to_edge(UL).shift(0.7*RIGHT+ 0.5*DOWN)
        l = t.get_left()[0]
        t1 = Tex(r"Let ", r"$\Sigma$", r" be a surface, ", r"$R $", r" $\subset \Sigma $ homeomorphic to a disc, such that ",r"$\partial R$", font_size = 34).next_to(t, DOWN)
        t2 = Tex(r"is a piece-wise smooth curve travelling with $R$ on its left. Then", font_size = 34).next_to(t1, DOWN)
        t3 = Tex(r"$ \int_{R} K +  \Psi (\partial R ) = 2 \pi $.", font_size = 34).next_to(t2, DOWN)
        hu = Tex(r"Theorem (Hopf Umlaufsatz)", font_size = 34, color = BLUE).next_to(t3, DOWN)
        hu1 = Tex(r"Let $\gamma : [a,b] \to \mathbb{R}^2$ be a simple closed piece-wise smooth curve.", font_size = 34).next_to(hu, DOWN)
        hu2 = Tex(r"If it travels with the region it surrounds on its left, then", font_size = 34).next_to(hu1, DOWN)
        hu3 = Tex(r"$\Psi ( \gamma ) = 2 \pi .$", font_size = 34).next_to(hu2, DOWN)
        lem = Tex(r"Lemma", font_size = 34, color = BLUE).next_to(t3, DOWN)
        l1 = Tex(r"Let ", r"$\Sigma$", r" be a surface, ", r"$R $", r" $\subset \Sigma $ homeomorphic to a disc, such that ",r"$\partial R$", font_size = 34).next_to(lem, DOWN)
        l2 = Tex(r"is a piece-wise smooth curve travelling with $R$ on its left. Then", font_size = 34).next_to(l1, DOWN)
        l3 = Tex(r"$ \int_{R} K +  \Psi (\partial R ) = 2 \pi n $ with $n \in \mathbb{Z}$.", font_size = 34).next_to(l2, DOWN)
        s1 = Tex(r"Step 1:", font_size = 30).next_to(l3, DOWN).shift(0.1*DOWN)
        s2 = Tex(r"Step 2:", font_size = 30).next_to(l3, DOWN).shift(0.1*DOWN)
        s3 = Tex(r"Step 3:", font_size = 30).next_to(l3, DOWN).shift(0.1*DOWN)
        s4 = Tex(r"Step 4:", font_size = 30).next_to(l3, DOWN).shift(0.1*DOWN)
        s11 = Tex(r"Lemma holds for", font_size = 30).next_to(s1, DOWN).shift(0.1*UP)
        s111 = Tex(r"small enough regions", font_size = 30).next_to(s11, DOWN).shift(0.1*UP)
        s22 = Tex(r"Lemma holds for", font_size = 30).next_to(s1, DOWN).shift(0.1*UP)
        s222 = Tex(r"all regions", font_size = 30).next_to(s11, DOWN).shift(0.1*UP)
        s33 = Tex(r"Theorem holds for", font_size = 30).next_to(s1, DOWN).shift(0.1*UP)
        s333 = Tex(r"small enough regions", font_size = 30).next_to(s11, DOWN).shift(0.1*UP)
        s44 = Tex(r"Theorem holds for", font_size = 30).next_to(s1, DOWN).shift(0.1*UP)
        s444 = Tex(r"all regions", font_size = 30).next_to(s11, DOWN).shift(0.1*UP)
        t1.shift((l - t1.get_left()[0])*RIGHT)
        t2.shift((l - t2.get_left()[0])*RIGHT)
        t3.shift(t3.get_center()[0]*LEFT)
        hu.shift((l - hu.get_left()[0])*RIGHT)
        hu1.shift((l - hu1.get_left()[0])*RIGHT)
        hu2.shift((l - hu2.get_left()[0])*RIGHT)
        hu3.shift(hu3.get_center()[0]*LEFT)
        lem.shift((l - lem.get_left()[0])*RIGHT)
        l1.shift((l - l1.get_left()[0])*RIGHT)
        l2.shift((l - l2.get_left()[0])*RIGHT)
        l3.shift(l3.get_center()[0]*LEFT)
        q = 1.6
        s1.shift(s1.get_center()[0]*LEFT + 3*q*LEFT)
        s2.shift(s2.get_center()[0]*LEFT + q*LEFT)
        s3.shift(s3.get_center()[0]*LEFT + q*RIGHT)
        s4.shift(s4.get_center()[0]*LEFT + 3*q*RIGHT)
        s11.shift(s11.get_center()[0]*LEFT + 3*q*LEFT)
        s22.shift(s22.get_center()[0]*LEFT + q*LEFT)
        s33.shift(s33.get_center()[0]*LEFT + q*RIGHT)
        s44.shift(s44.get_center()[0]*LEFT + 3*q*RIGHT)
        s111.shift(s111.get_center()[0]*LEFT + 3*q*LEFT)
        s222.shift(s222.get_center()[0]*LEFT + q*LEFT)
        s333.shift(s333.get_center()[0]*LEFT + q*RIGHT)
        s444.shift(s444.get_center()[0]*LEFT + 3*q*RIGHT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(t1[1], DOWN).shift(0.15*UP)
        rli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(t1[3], DOWN).shift(0.15*UP)
        prli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(t1[5], DOWN).shift(0.15*UP)
        sli1 = Line([0,0,0], [0.8,0,0], color = TEAL).next_to(s1, DOWN).shift(0.18*UP)
        sli2 = Line([0,0,0], [0.8,0,0], color = TEAL).next_to(s2, DOWN).shift(0.18*UP)
        sli3 = Line([0,0,0], [0.8,0,0], color = TEAL).next_to(s3, DOWN).shift(0.18*UP)
        sli4 = Line([0,0,0], [0.8,0,0], color = TEAL).next_to(s4, DOWN).shift(0.18*UP)
        h = - 1.2
        ran = PI
        op = 0.2
        def f(x,y):
            return  ( np.sin(2*x) - np.cos(2*y))/6 + h
        g1 = ParametricFunction(lambda t : [ ran /2 , t , f(ran /2 , t)  ] , t_range = [ - ran /2, ran/2 ], color = ORANGE)
        g2 = ParametricFunction(lambda t : [ -ran /2 , -t , f(-ran /2 , -t)  ] , t_range = [ - ran /2, ran/2 ], color = ORANGE)
        g3 = ParametricFunction(lambda t : [ -t, ran /2 , f( -t, ran /2)  ] , t_range = [ - ran /2, ran/2 ], color = ORANGE)
        g4 = ParametricFunction(lambda t : [ t, -ran /2 , f( t, -ran /2)  ] , t_range = [ - ran /2, ran/2 ], color = ORANGE)
        g = VGroup(g1, g2, g3, g4)
        s = Surface(
            lambda u, v:  [ u, v, f(u,v) ],
            u_range = [ - ran , ran ],
            v_range = [ - ran , ran ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op*0.75 #, 
            #resolution = [ 16, 16]
        )
        r = Surface(
            lambda u, v:  [ u, v, f(u,v) ],
            u_range = [ - ran /2 , ran /2 ],
            v_range = [ - ran /2 , ran /2 ],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        self.add_fixed_in_frame_mobjects( t, t1, t2, t3, hu, hu1, hu2, hu3, lem, l1, l2, l3, s1, s11, s111, s2, s22, s222, s3, s33, s333, s4, s44, s444, sli, rli, prli, sli1, sli2, sli3, sli4)
        self.remove(t, t1, t2, t3, hu, hu1, hu2, hu3, lem, l1, l2, l3, s1, s11, s111, s2, s22, s222, s3, s33, s333, s4, s44, s444, sli, rli, prli, sli1, sli2, sli3, sli4)
        self.play(Write(t), Write(t1), Write(t2), Write(t3), Create(sli), Create(rli), Create(prli), Create(s), Create(r), Create(g))
        self.wait()
        self.play(FadeOut(s), FadeOut(r), FadeOut(g))
        self.wait()
        self.play(Write(hu), Write(hu1), Write(hu2), Write(hu3))
        self.wait()
        hopf = VGroup(hu, hu1, hu2, hu3)
        self.play(FadeOut(hopf))
        self.wait()
        lemma = VGroup(lem, l1, l2, l3)
        self.play(Write(lemma))
        self.wait()
        self.play(Create(s1), Create(s2), Create(s3), Create(s4), Create(sli1), Create(sli2), Create(sli3), Create(sli4))
        self.wait()
        self.play(Create(s11), Create(s111))
        self.wait()
        self.play(Create(s22), Create(s222))
        self.wait()
        self.play(Create(s33), Create(s333))
        self.wait()
        self.play(Create(s44), Create(s444))
        self.wait()



class tr(Scene):
    def construct(self):
        pro = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift(0.7*RIGHT+ 0.3*DOWN)
        l = pro.get_left()[0]
        pro1 = Tex(r"Let $\gamma  : [a,b ] \to \Sigma$ be a unit-speed closed piece-wise smooth curve,", font_size = 34).next_to(pro, DOWN).shift(0.05*UP)
        pro2 = Tex(r"$X : [a,b] \to \mathbb{R}^3$ a unit vector field along $\gamma$ with $X(a) = X(b) $. ", font_size = 34).next_to(pro1, DOWN).shift(0.05*UP)
        pro3 = Tex(r"Let $Y = N \times X$. Then ", font_size = 34).next_to(pro2, DOWN).shift(0.05*UP)
        pro4 = Tex(r"$\Psi (\gamma ) = \int_a^b ( X^{\prime} \cdot Y ) \text{d}t $ \,\, mod $2 \pi$.", font_size = 34).next_to(pro3, DOWN).shift(0.05*UP)
        prf = Tex(r"Proof", font_size = 34, color = BLUE).next_to(pro4, DOWN)
        prf1 = Tex(r"Let $\overline{X}, \overline{Y}: [a,b] \to \mathbb{R}^3$ be parallel vector fields along $\gamma$", font_size = 34).next_to(prf, DOWN).shift(0.05*UP)
        prf2 = Tex(r"with $\overline{X}(a) = X(a), $ $\overline{Y} (a) = Y(a)$.", r" Then", font_size = 34).next_to(prf1, DOWN).shift(0.05*UP)
        prf3 = Tex(r"$X (t)  =  \cos (\theta(t)) \overline{X} (t) + \sin (\theta(t)) \overline{Y}(t) $.", font_size = 34).next_to(prf2, DOWN).shift(0.05*UP)
        prf4 = Tex(r"$X^{\prime} \cdot Y  = \theta^{\prime }$", r" $\Rightarrow$ $\int_a^b ( X^{\prime} \cdot Y ) \text{d}t = \theta (b)$.", font_size = 34).next_to(prf3, DOWN).shift(0.05*UP)
        prf5 = Tex(r" $\Rightarrow$ $P_{\gamma} (X(a))$ is a clockwise rotation of $X(a)$ by angle $\int_a^b (X^{\prime }\cdot Y ) \text{d} t $.", font_size = 34).next_to(prf4, DOWN).shift(0.05*UP)
        pro1.shift((l - pro1.get_left()[0])*RIGHT)
        pro2.shift((l - pro2.get_left()[0])*RIGHT)
        pro3.shift((l - pro3.get_left()[0])*RIGHT)
        pro4.shift(pro4.get_center()[0]*LEFT)
        prf.shift((l - prf.get_left()[0])*RIGHT)
        prf1.shift((l - prf1.get_left()[0])*RIGHT)
        prf2.shift((l - prf2.get_left()[0])*RIGHT)
        prf3.shift(prf3.get_center()[0]*LEFT)
        prf4.shift((l - prf4.get_left()[0])*RIGHT)
        prf5.shift((l - prf5.get_left()[0])*RIGHT)
        prop = VGroup(pro, pro1, pro2, pro3, pro4)
        self.play(Create(prop))
        self.wait()
        self.play(Write(prf), Write(prf1), Write(prf2[0]))
        self.wait()
        self.play(Write(prf2[1]), Write(prf3))
        self.wait()
        self.play(Write(prf4))
        self.wait()
        self.play(Write(prf5))
        self.wait()



class step1(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=78 * DEGREES, theta=  PI/3)
        t = Tex(r"Step 1", font_size = 34).to_edge(UL).shift( RIGHT+ 0.6* DOWN)
        sli1 = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(t, DOWN).shift(0.15*UP)
        l = t.get_left()[0]
        a = 1.5
        t1 = Tex(r"Assume ", r"$R$", r" is coved by a semi-geodesic chart ", r"$s: U \to \Sigma$.", font_size = 34).next_to(t, DOWN)
        t1.shift((l - t1.get_left()[0])*RIGHT)
        rli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(t1[1], DOWN).shift(0.15*UP)
        sar = CurvedArrow( [-1.2 ,-0.4 ,0],[0.8 ,-0.2,0], radius = -3 , tip_length = 0.2)
        sarl = Tex(r"$s$", font_size = 34).next_to(sar, UP)
        ht = 2
        r = 1
        a = 0.3
        m = 1.2
        h = -1
        hu = -2
        op = 0.2
        g1 = ParametricFunction(lambda t :  [ (r + a*np.cos(t) ) * np.sin(9*PI/16)  -ht*sqrt(3) , t  + ht - PI/2 , (r + a*np.cos(t) ) * np.cos(9*PI/16) + h ] , t_range = [ 0, PI ], color = ORANGE)
        g2 = ParametricFunction(lambda t :  [ (r + a*np.cos(PI) ) * np.sin(9*PI/16 - t)  -ht*sqrt(3) , PI  + ht - PI/2 , (r + a*np.cos(PI) ) * np.cos(9*PI/16 - t) + h ] , t_range = [ 0, 3*PI/8 ], color = ORANGE)
        g3 = ParametricFunction(lambda t :  [ (r + a*np.cos(PI - t) ) * np.sin(3*PI/16)  -ht*sqrt(3) , PI - t  + ht - PI/2 , (r + a*np.cos(PI - t) ) * np.cos(3*PI/16) + h ] , t_range = [ 0, PI ], color = ORANGE)
        g4 = ParametricFunction(lambda t :  [ (r + a*np.cos(0) ) * np.sin(3*PI/16 + t)  -ht*sqrt(3) ,  ht - PI/2 , (r + a*np.cos(0) ) * np.cos(3*PI/16 + t) + h ] , t_range = [ 0 , 3*PI /8 ], color = ORANGE)
        ht = 1.8
        g = VGroup(g1, g2, g3, g4)
        gm1 = ParametricFunction(lambda t : [ t/m + ht *sqrt(3) , -PI/(2*m) - ht , h    ] , t_range = [  -3*PI/16 , 3*PI / 16 ] , color = ORANGE )
        gm2 = ParametricFunction(lambda t : [ (3*PI/16)/m + ht *sqrt(3) , t/m - ht , h    ] , t_range = [  -PI/2 , PI / 2 ] , color = ORANGE )
        gm3 = ParametricFunction(lambda t : [ (3*PI/16 - t)/m + ht *sqrt(3) , PI/(2*m) - ht , h    ] , t_range = [  0 , 3*PI / 8 ] , color = ORANGE )
        gm4 = ParametricFunction(lambda t : [ (-3*PI/16)/m + ht *sqrt(3) , (PI/2 - t)/m - ht , h    ] , t_range = [  0 , PI ] , color = ORANGE )
        gm = VGroup(gm1, gm2, gm3, gm4)
        u = Surface(
            lambda u, v:  [ v/m + ht *sqrt(3) , u/m - ht , h    ],
            u_range = [ - PI  , PI  ],
            v_range = [  -3*PI/8 , 3*PI / 8 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op*0.75 #, 
            #resolution = [ 16, 16]
        )
        sr = Surface(
            lambda u, v:  [ v/m + ht*sqrt(3), u/m - ht , h    ],
            u_range = [ - PI /2  , PI /2  ],
            v_range = [  -3*PI/16 , 3*PI / 16 ],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        ht = 2
        s = Surface(
            lambda u, v:  [ (r + a*np.cos(u) ) * np.sin(v) - ht*sqrt(3) , u + ht - PI/2 , (r + a*np.cos(u) ) * np.cos(v) + h  ],
            u_range = [ - PI / 2 , 3*PI /2 ],
            v_range = [  0 , 3*PI / 4 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op*0.75 #, 
            #resolution = [ 16, 16]
        )
        r = Surface(
            lambda u, v:  [ (r + a*np.cos(u) ) * np.sin(v)  -ht*sqrt(3) , u  + ht - PI/2 , (r + a*np.cos(u) ) * np.cos(v) + h ],
            u_range = [  0 , PI  ],
            v_range = [  3*PI/16 , 9*PI / 16  ],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        fig = VGroup(s, r, u, sr)
        self.add_fixed_in_frame_mobjects( t, t1, rli, sar, sarl, sli1)
        self.remove(t, t1, rli,  sar, sarl, sli1)
        self.play(Write(t), Write(t1), Create(rli), Create(fig), Create(g), Create(gm), Create(sar), Create(sarl), Create(sli1))
        self.wait()
        self.play(FadeOut(fig), FadeOut(g), FadeOut(gm), FadeOut(sar), FadeOut(sarl))
        self.wait()



class step11(Scene):
    def construct(self):
        t = Tex(r"Step 1", font_size = 34 ).to_edge(UL).shift( RIGHT+ 0.6* DOWN)
        sli1 = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(t, DOWN).shift(0.15*UP)
        l = t.get_left()[0]
        a = 1.5
        t1 = Tex(r"Assume ", r"$R$", r" is coved by a semi-geodesic chart ", r"$s: U \to \Sigma$.", font_size = 34).next_to(t, DOWN)
        t2 = Tex(r"$\vert s_u \vert = 1$, $\vert s_v  \vert = b $, $s_u \cdot s_v = 0$, then Jac$(s) = b$.", font_size = 34).next_to(t1, DOWN)
        t3 = Tex(r"$ bK +  b_{uu} = 0$, ", r" $X_u \cdot Y = 0$, ", r" $X_v \cdot Y = b_u $.", font_size = 34).next_to(t2, DOWN)
        t4 = Tex(r"$\int_R K  = \int_{s^{-1}(R)} bK$ ", r"$ = $", r" $ - \int_{s^{-1}(R)} b_{uu} $ ", font_size = 34).next_to(t3, DOWN).shift(0.3*DOWN)
        t5 = Tex(r"$ =  - \int_{s^{-1}(\partial R)} $ ", r"$ b_u $", r" d$v $",  font_size = 34).next_to(t4, DOWN)
        t6 = Tex(r"$ = - \int_{s^{-1} (\partial R )}$ ", r"$ ( (X_u \cdot Y)$", r" $  u^{\prime}  + $ ", r"$(X_v \cdot Y )$", r" $ v^{\prime} ) \text{d}t $ ",  font_size = 34).next_to(t5, DOWN)
        t7 = Tex(r"$ = - \int_{\partial R} (X^{\prime }\cdot Y )$ ", font_size = 34).next_to(t6, DOWN)
        t8 = Tex(r"$ = - \Psi (\partial R ) + 2 \pi n $.", font_size = 34).next_to(t7, DOWN)
        t1.shift((l - t1.get_left()[0])*RIGHT)
        t2.shift((l - t2.get_left()[0])*RIGHT)
        t3.shift(t3.get_center()[0]*LEFT)
        t4.shift((l - t4.get_left()[0])*RIGHT)
        t5.shift((l - t5.get_left()[0] + a )*RIGHT)
        t6.shift((l - t6.get_left()[0] + a )*RIGHT)
        t7.shift((l - t7.get_left()[0] + a )*RIGHT)
        t8.shift((l - t8.get_left()[0] + a )*RIGHT)
        rli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(t1[1], DOWN).shift(0.15*UP)
        jli = Line([0,0,0], [1.7,0,0], color = PURPLE).next_to(t3[0], DOWN).shift(0.15*UP)
        jli2 = Line([0,0,0], [0.5,0,0], color = PURPLE).next_to(t4[1], DOWN).shift(0.15*UP)
        xuli = Line([0,0,0], [1.7,0,0], color = BLUE).next_to(t3[1], DOWN).shift(0.15*UP)
        xuli2 = Line([0,0,0], [1,0,0], color = BLUE).next_to(t6[1], DOWN).shift(0.15*UP)
        xvli = Line([0,0,0], [1.7,0,0], color = GREEN).next_to(t3[2], DOWN).shift(0.15*UP)
        xvli2 = Line([0,0,0], [0.45,0,0], color = GREEN).next_to(t5[1], DOWN).shift(0.15*UP)
        xvli3 = Line([0,0,0], [1,0,0], color = GREEN).next_to(t6[3], DOWN).shift(0.15*UP)
        self.add(t, t1, rli, sli1)
        self.wait()
        self.play(Write(t2))
        self.wait()
        self.play(Write(t3), Create(jli), Create(xuli), Create(xvli))
        self.wait()
        self.play(Write(t4[0]))
        self.wait()
        self.play(Write(t4[1]), Write(t4[2]), Create(jli2))
        self.wait()
        self.play(Write(t5))
        self.wait()
        self.play(Write(t6), Create(xuli2), Create(xvli2), Create(xvli3))
        self.wait()
        self.play(Write(t7))
        self.wait()
        self.play(Write(t8))
        self.wait()



class step12(Scene):
    def construct(self):
        shi = 0.15
        t1 = Tex(r"GB$(R) : = \int_R K + \Psi (\partial R) -2 \pi$", font_size = 34).shift(2.7*UP)
        t2 = Tex(r"Step 1", r" : If $R$ is covered by a semi-geodesic chart, then GB$(R) = 0$ mod $2\pi$.", font_size = 34).next_to(t1, DOWN).shift(0.8*LEFT)
        l = t2.get_left()[0]
        p = Tex(r"Proposition (Additivity)", font_size = 34, color = BLUE).next_to(t2, DOWN)
        p1 = Tex(r"If by cutting ", r"$R$", r" along a curve ", r"$\alpha$", r", we obtain two regions ", r"$R_1$", r" , ", r"$R_2$", r", then ", font_size = 34).next_to(p, DOWN).shift(0.1*UP)
        p2 = Tex(r"GB$(R) = $GB$(R_1) + $GB$(R_2)$.",  font_size = 34).next_to(p1, DOWN)
        pr = Tex(r"Proof", font_size = 34, color = BLUE).next_to(p2, DOWN).shift(shi*UP)
        pr1 = Tex(r"$\partial R =$ ", r"$ \gamma _1$",  r" $+$ ", r"$ \gamma _2$", r", $\, \partial R_1 = $ ", r"$\gamma_1$", r" $+$ ", r"$\alpha$", r", $\, \partial R_2 = $ ", r"$\gamma_2$", r" $ - $ ", r"$\alpha$.", font_size = 34).next_to(pr, DOWN)
        pr2 = Tex(r"$\Psi (\partial R_1) = $ ", r"$\Psi (\gamma_1)$", r" $ +$ ", r"$ \Psi (\alpha)$", r" $ + (\pi - \theta _1) + (\pi - \theta_2)$", font_size = 34).next_to(pr1, DOWN)
        pr3 = Tex(r"$\Psi (\partial R_2) = $ ", r"$\Psi (\gamma_2)$", r" $ -$ ", r"$ \Psi (\alpha)$", r" $ + (\pi - \theta _3) + (\pi - \theta_4)$", font_size = 34).next_to(pr2, DOWN)
        pr4 = Tex(r"$\Psi (\partial R  ) = $ ", r"$\Psi (\gamma_1)$", r" $ +$ ", r"$ \Psi (\gamma_2)$", r" $ + (\pi - \theta _1 - \theta_3) + (\pi - \theta_2 - \theta_4)$", font_size = 34).next_to(pr3, DOWN)
        pr5 = Tex(r"$\Rightarrow$ ", r"$\Psi (\partial R) = \Psi (\partial R_1) + \Psi (\partial R_2 ) - 2 \pi $", font_size = 34).next_to(pr4, DOWN)
        pr6 = Tex(r"GB$(R) = \int_R K + $ ", r"$\Psi (\partial R)$", r" $ - 2 \pi $", font_size = 34).next_to(pr1, DOWN)
        pr7 = Tex(r"$= \int_{R_1} K + \int_{R_2} K  + $ ", r"$  \Psi (\partial R_1 ) + \Psi (\partial R_2 ) - 2 \pi $", r" $- 2 \pi  $", font_size = 34).next_to(pr6, DOWN)
        pr8 = Tex(r"$= $ GB$(R_1) + $GB$(R_2)$", font_size = 34).next_to(pr7, DOWN)
        p.shift((l - p.get_left()[0] )*RIGHT)
        p1.shift((l - p1.get_left()[0] )*RIGHT)
        p2.shift(p2.get_center()[0]*LEFT)
        pr.shift((l - pr.get_left()[0] )*RIGHT)
        pr1.shift((l - pr1.get_left()[0] )*RIGHT)
        pr2.shift((l - pr2.get_left()[0] )*RIGHT)
        pr3.shift((l - pr3.get_left()[0] )*RIGHT)
        pr4.shift((l - pr4.get_left()[0] )*RIGHT)
        pr6.shift((l - pr6.get_left()[0] )*RIGHT)
        pr7.shift((l - pr7.get_left()[0] + 0.5 )*RIGHT)
        pr8.shift((l - pr8.get_left()[0] + 0.5 )*RIGHT)
        pr5.shift(pr5.get_center()[0]*LEFT)
        s1li = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(t2[0], DOWN).shift(0.15*UP)
        ali = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(p1[3], DOWN).shift(0.1*UP)
        r1li = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(p1[5], DOWN).shift(0.15*UP)
        r2li = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(p1[7], DOWN).shift(0.15*UP)
        g1li1 = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(pr1[1], DOWN).shift(0.15*UP)
        g1li2 = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(pr1[5], DOWN).shift(0.15*UP)
        g2li1 = Line([0,0,0], [0.4,0,0], color = MAROON_E).next_to(pr1[3], DOWN).shift(0.15*UP)
        g2li2 = Line([0,0,0], [0.4,0,0], color = MAROON_E).next_to(pr1[9], DOWN).shift(0.15*UP)
        ali1 = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(pr1[7], DOWN).shift(0.1*UP)
        ali2 = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(pr1[11], DOWN).shift(0.1*UP)
        g1li3 = Line([0,0,0], [0.7,0,0], color = ORANGE).next_to(pr2[1], DOWN).shift(0.15*UP)
        g1li4 = Line([0,0,0], [0.7,0,0], color = ORANGE).next_to(pr4[1], DOWN).shift(0.15*UP)
        g2li3 = Line([0,0,0], [0.7,0,0], color = MAROON_E).next_to(pr3[1], DOWN).shift(0.15*UP)
        g2li4 = Line([0,0,0], [0.7,0,0], color = MAROON_E).next_to(pr4[3], DOWN).shift(0.15*UP)
        ali3 = Line([0,0,0], [0.6,0,0], color = BLUE).next_to(pr2[3], DOWN).shift(0.15*UP)
        ali4 = Line([0,0,0], [0.6,0,0], color = BLUE).next_to(pr3[3], DOWN).shift(0.15*UP)
        fil1 = Line([0,0,0], [1.2,0,0], color = PINK).next_to(pr1, DOWN).shift(0.15*UP + 2.3*LEFT)
        fil2 = Line([0,0,0], [1,0,0], color = PINK).next_to(pr6[1], DOWN).shift(0.15*UP)
        fil3 = Line([0,0,0], [3.5,0,0], color = PINK).next_to(pr7[1], DOWN).shift(0.15*UP)
        #ax = Axes( x_range=[-4, 4],  y_range=[-4, 4]  )
        a = sqrt(3)/2
        #def c(x):
        #   return (1 - x^2)**(1/2)
        #r11 = ax.get_area(gra , [-a, 1], color=YELLOW, opacity=0.15)
        r = 1
        ax = Axes( x_range=[-5, 5],  y_range=[0, 5] , x_length = 10*r, y_length = 5*r)
        def c(x):
           return (1-x**2)**(1/2)
        gr = ax.plot(c, x_range=[-1, 1])
        r11 = ax.get_area(gr, [-a, 1], color= PURPLE_B, opacity=0.2, stroke_width = 0)
        r21 = ax.get_area(gr, [-1, a], color= YELLOW, opacity=0.2, stroke_width = 0)
        r11.shift( ( 4 - r11.get_left()[0] )*RIGHT  + (r11.get_bottom()[1] + 1.5)*DOWN )
        r21.shift( ( 4 - r21.get_right()[0] )*RIGHT  + (r21.get_bottom()[1] + 1.5)*DOWN )
        r12 = r21.copy()
        r22 = r11.copy()
        r12.rotate( PI , about_point = [ r12.get_right()[0], r12.get_bottom()[1], 0 ] ).set_color(PURPLE_B)
        r22.rotate( PI , about_point = [ r11.get_left()[0], r11.get_bottom()[1], 0 ] ).set_color(YELLOW)
        alp = Arrow([ 4, - 2 , 0 ], [ 4, -1, 0  ], color = BLUE, buff = 1.5, max_tip_length_to_length_ratio=0.15)
        r2l = Tex(r"$R_2$", font_size = 34).next_to(alp, RIGHT).shift(0.4*RIGHT)
        r1l = Tex(r"$R_1$", font_size = 34).next_to(alp, LEFT).shift(0.4*LEFT)
        g1 =  ParametricFunction(lambda t :  [ 4 + np.cos(t) - sqrt(3)/2 , -1.5 + np.sin(t) , 0 ] , t_range = [ PI/6, 11*PI/6 ], color = ORANGE)
        g2 =  ParametricFunction(lambda t :  [ 4 - np.cos(t) + sqrt(3)/2, -1.5 - np.sin(t) , 0 ] , t_range = [ PI/6, 11*PI/6 ], color = MAROON_E)
        z = 0.1
        g1t1 = Line([3 - sqrt(3)/2 , -1.5, 0], [3 - sqrt(3)/2 - z ,-1.5 + z ,0], color = ORANGE)
        g1t2 = Line([3 - sqrt(3)/2 , -1.5, 0], [3 - sqrt(3)/2 + z ,-1.5 + z ,0], color = ORANGE)
        g2t1 = Line([5 + sqrt(3)/2 , -1.5, 0], [5 + sqrt(3)/2 - z ,-1.5 - z,0], color = MAROON_E)
        g2t2 = Line([5 + sqrt(3)/2 , -1.5, 0], [5 + sqrt(3)/2 + z ,-1.5 - z,0], color = MAROON_E)
        a1 = Line([ 4, - 2 , 0 ], [ 4, - 1, 0  ])
        a2 = Line([ 4, - 1 , 0 ], [ 4, - 2, 0  ])
        g1n = Line([4, - 2 , 0 ], [ 3, - 2 - sqrt(3), 0] )
        g2p = Line([4, - 2 , 0 ], [ 5, - 2 - sqrt(3), 0] )
        g1p = Line([4, - 1 , 0 ], [ 3, - 1 + sqrt(3), 0] )
        g2n = Line([4, - 1 , 0 ], [ 5, - 1 + sqrt(3), 0] )
        thr = 0.15
        th1 = Angle( g1p, a2 , radius = thr, color = ORANGE)
        th2 = Angle( a1 , g1n, radius = thr, color = ORANGE )
        th3 = Angle( a2 , g2n, radius = thr, color = PURPLE_B )
        th4 = Angle( g2p, a1 , radius = thr, color = PURPLE_B )
        th1l = Tex(r"$\theta_1$", font_size = 25).next_to(th1, LEFT).shift(0.1*RIGHT)
        th2l = Tex(r"$\theta_2$", font_size = 25).next_to(th2, LEFT).shift(0.1*RIGHT)
        th3l = Tex(r"$\theta_3$", font_size = 25).next_to(th3, RIGHT).shift(0.1*LEFT)
        th4l = Tex(r"$\theta_4$", font_size = 25).next_to(th4, RIGHT).shift(0.1*LEFT)
        tips = VGroup(g1t1, g1t2, g2t1, g2t2)
        theta = VGroup(th1, th2, th3, th4, th1l, th2l, th3l, th4l)
        self.play(Write(t1))
        self.wait()
        self.play(Write(t2), Create(s1li))
        self.wait()
        self.play(Write(p), Write(p1), Write(p2), Create(ali), Create(r11), Create(r12), Create(r21), Create(r22), Create(r1li), Create(r2li), Create(alp), Create(g1), Create(g2), Create(r1l), Create(r2l), Create(tips))
        self.wait()
        self.play(Write(pr), Write(pr1), Create(ali1), Create(ali2), Create(g1li1), Create(g1li2), Create(g2li1), Create(g2li2))
        self.wait()
        self.play(Create(theta), r1l.animate.shift(0.1*LEFT), r2l.animate.shift(0.1*RIGHT))
        self.add(alp)
        self.wait()
        self.play(Write(pr2), Write(pr3), Write(pr4))
        self.wait()
        self.play(Create(g1li3), Create(g1li4), Create(g2li3), Create(g2li4), Create(ali3), Create(ali4))
        self.wait()
        self.play(Write(pr5))
        self.wait()
        draft = VGroup(pr1, pr2, pr3, pr4, ali1, ali2, ali3, ali4, g1li1, g1li2, g1li3, g1li4, g2li1, g2li2, g2li3, g2li4, pr5[0])
        self.play(FadeOut(draft), pr5[1].animate.shift(( l - pr5[1].get_left()[0] )*RIGHT + (pr1.get_center()[1] - pr5[1].get_center()[1])*UP), Create(fil1))
        self.wait()
        self.play(Write(pr6))
        self.wait()
        self.play(Write(pr7), Create(fil2), Create(fil3))
        self.wait()
        self.play(Write(pr8))
        self.wait()




class step2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=78 * DEGREES, theta=  PI/4)
        t1 = Tex(r"GB$(R) : = \int_R K + \Psi (\partial R) -2 \pi$", font_size = 34).shift(2.7*UP)
        t2 = Tex(r"Step 1", r" : If $R$ is covered by a semi-geodesic chart, then GB$(R) = 0$ mod $2\pi$.", font_size = 34).next_to(t1, DOWN).shift(0.8*LEFT)
        l = t2.get_left()[0]
        s = Tex(r"Step 2", r" : GB$(R) = 0$ mod $2 \pi$ for all regions as in the theorem.", font_size = 34).next_to(t2, DOWN)
        s1 = Tex(r"$R = \bigcup_{i = 1}^n R_i$ with $R_i$ inside a semi-geodesic chart.", font_size = 34).next_to(s, DOWN)
        s2 = Tex(r"GB$(R) = \sum_{i=1}^n $GB$(R_i)$ ", r"$= 0 $ mod $2 \pi$.",  font_size = 34).next_to(s1, DOWN)
        s.shift((l - s.get_left()[0] )*RIGHT)
        s1.shift((l - s1.get_left()[0] )*RIGHT)
        s2.shift((l - s2.get_left()[0] )*RIGHT)
        s1li = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(t2[0], DOWN).shift(0.15*UP)
        s2li = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(s[0], DOWN).shift(0.15*UP)
        h = - 1.2
        ran = PI/2
        op = 0.2
        def f(x,y):
            return  ( np.sin(2*x) - np.cos(2*y))/6 + h
        g1 = ParametricFunction(lambda t : [ ran  , t , f(ran  , t)  ] , t_range = [ - ran , ran ], color = ORANGE)
        g2 = ParametricFunction(lambda t : [ -ran , -t , f(-ran  , -t)  ] , t_range = [ - ran , ran ], color = ORANGE)
        g3 = ParametricFunction(lambda t : [ -t, ran  , f( -t, ran )  ] , t_range = [ - ran , ran ], color = ORANGE)
        g4 = ParametricFunction(lambda t : [ t, -ran  , f( t, -ran )  ] , t_range = [ - ran , ran ], color = ORANGE)
        g = VGroup(g1, g2, g3, g4)
        si = Surface(
            lambda u, v:  [ u, v, f(u,v) ],
            u_range = [ - 2*ran , 2*ran ],
            v_range = [ - 2*ran , 2*ran ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op*0.75 #, 
            #resolution = [ 16, 16]
        )
        r = Surface(
            lambda u, v:  [ u, v, f(u,v) ],
            u_range = [ - ran  , ran  ],
            v_range = [ - ran  , ran  ],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        m = 4 #number of pieces (square root)
        ri = []
        gami = []
        risur = VGroup()
        gamic = VGroup()
        block = []
        for i in range(0,m):
            ri.append([])
            gami.append([])
            block.append([])
            for j in range(0,m):
                block[i].append(VGroup())
                ri[i].append( Surface(
                    lambda u, v:  [ u, v, f(u,v) ],
                    u_range = [ - ran  + 2* i * ran /m  , - ran + 2* (i+1) * ran /m  ],
                    v_range = [ - ran  + 2* j * ran /m  , - ran + 2* (j+1) * ran /m  ],
                    checkerboard_colors = [YELLOW, YELLOW],
                    fill_opacity = op, 
                    resolution = [ 4, 4]
                    )
                )
                risur.add(ri[i][j])
                block[i][j].add(ri[i][j])
                gami[i].append([])
                gami[i][j].append( ParametricFunction(
                    lambda t : [ - ran  + 2* i * ran /m  , t , f( - ran  + 2* i * ran /m , t)  ] , 
                    t_range = [ - ran  + 2* j * ran /m  , - ran + 2* (j+1) * ran /m ], 
                    color = ORANGE)
                )
                gami[i][j].append( ParametricFunction(
                    lambda t : [ - ran  + 2* (i+1) * ran /m  , t , f( - ran  + 2* (i+1) * ran /m , t)  ] , 
                    t_range = [ - ran  + 2* j * ran /m  , - ran + 2* (j+1) * ran /m ], 
                    color = ORANGE)
                )
                gami[i][j].append( ParametricFunction(
                    lambda t : [ t ,  - ran  + 2* j * ran /m  , f( t ,  - ran  + 2* j * ran /m )  ] , 
                    t_range = [ - ran  + 2* i * ran /m  , - ran + 2* (i+1) * ran /m ], 
                    color = ORANGE)
                )
                gami[i][j].append( ParametricFunction(
                    lambda t : [ t ,  - ran  + 2* (j+1) * ran /m  , f( t ,  - ran  + 2* (j+1) * ran /m )  ] , 
                    t_range = [ - ran  + 2* i * ran /m  , - ran + 2* (i+1) * ran /m ], 
                    color = ORANGE)
                )
                for w in range(0,4):
                    block[i][j].add(gami[i][j][w])
                    gamic.add(gami[i][j][w])
        self.add_fixed_in_frame_mobjects(t1, t2, s, s1, s2, s1li, s2li)
        self.remove(s, s1, s2, s2li)
        self.wait()
        self.play(Write(s), Create(s2li), Create(si), Create(g), Create(r))
        self.wait()
        self.add(risur)
        self.remove(r)
        self.play(Write(s1), FadeOut(si), Create(gamic))
        self.wait()
        self.remove(g)
        ht =1 #horizontal translation
        self.play(
            block[0][0].animate.shift( [ ( 0 - 1.5) * ht , ( 0 - 1.5) * ht , 0 ]  ),
            block[0][1].animate.shift( [ ( 0 - 1.5) * ht , ( 1 - 1.5) * ht , 0 ]  ),
            block[0][2].animate.shift( [ ( 0 - 1.5) * ht , ( 2 - 1.5) * ht , 0 ]  ),
            block[0][3].animate.shift( [ ( 0 - 1.5) * ht , ( 3 - 1.5) * ht , 0 ]  ),
            block[1][0].animate.shift( [ ( 1 - 1.5) * ht , ( 0 - 1.5) * ht , 0 ]  ),
            block[1][1].animate.shift( [ ( 1 - 1.5) * ht , ( 1 - 1.5) * ht , 0 ]  ),
            block[1][2].animate.shift( [ ( 1 - 1.5) * ht , ( 2 - 1.5) * ht , 0 ]  ),
            block[1][3].animate.shift( [ ( 1 - 1.5) * ht , ( 3 - 1.5) * ht , 0 ]  ),
            block[2][0].animate.shift( [ ( 2 - 1.5) * ht , ( 0 - 1.5) * ht , 0 ]  ),
            block[2][1].animate.shift( [ ( 2 - 1.5) * ht , ( 1 - 1.5) * ht , 0 ]  ),
            block[2][2].animate.shift( [ ( 2 - 1.5) * ht , ( 2 - 1.5) * ht , 0 ]  ),
            block[2][3].animate.shift( [ ( 2 - 1.5) * ht , ( 3 - 1.5) * ht , 0 ]  ),
            block[3][0].animate.shift( [ ( 3 - 1.5) * ht , ( 0 - 1.5) * ht , 0 ]  ),
            block[3][1].animate.shift( [ ( 3 - 1.5) * ht , ( 1 - 1.5) * ht , 0 ]  ),
            block[3][2].animate.shift( [ ( 3 - 1.5) * ht , ( 2 - 1.5) * ht , 0 ]  ),
            block[3][3].animate.shift( [ ( 3 - 1.5) * ht , ( 3 - 1.5) * ht , 0 ]  )
        )
        self.wait()
        self.play(Write(s2))
        self.wait()



class step3(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=78 * DEGREES, theta=  PI/4)
        t1 = Tex(r"GB$(R) : = \int_R K + \Psi (\partial R) -2 \pi$", font_size = 34).shift(2.7*UP)
        t2 = Tex(r"Step 1", r" : If $R$ is covered by a semi-geodesic chart, then GB$(R) = 0$ mod $2\pi$.", font_size = 34).next_to(t1, DOWN).shift(0.8*LEFT)
        l = t2.get_left()[0]
        s = Tex(r"Step 2", r" : GB$(R) = 0$ mod $2 \pi$ for all regions as in the theorem.", font_size = 34).next_to(t2, DOWN)
        s1 = Tex(r"$R = \bigcup_{i = 1}^n R_i$ with $R_i$ inside a semi-geodesic chart.", font_size = 34).next_to(s, DOWN)
        s2 = Tex(r"GB$(R) = \sum_{i=1}^n $GB$(R_i)$ ", r"$= 0 $ mod $2 \pi$.",  font_size = 34).next_to(s1, DOWN)
        w = Tex(r"Step 3", r" : If $R$ is contained in the graph of a smooth function, then GB$(R) = 0 .$",   font_size = 34).next_to(s, DOWN).shift(UP)
        w1 = Tex(r"By Step 2, GB$(R)$ doesn't change during the process. Then GB$(R) = 0$.",   font_size = 34).next_to(w, DOWN)
        s.shift((l - s.get_left()[0] )*RIGHT)
        s1.shift((l - s1.get_left()[0] )*RIGHT)
        s2.shift((l - s2.get_left()[0] )*RIGHT)
        w.shift((l - w.get_left()[0] )*RIGHT)
        w1.shift((l - w1.get_left()[0] )*RIGHT)
        s1li = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(t2[0], DOWN).shift(0.15*UP)
        s2li = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(s[0], DOWN).shift(0.15*UP)
        s3li = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(w[0], DOWN).shift(0.15*UP)
        h = - 0.8
        nh = - 2
        ran = PI/2
        op = 0.2
        axes = ThreeDAxes(x_range = [-4,4,1], y_range = [-4, 4,1], z_range = [-1/2,2,1], x_length = 8, y_length = 8, z_length = 2.5 , tips = False )
        axes.shift([ 0, 0 , (nh - axes.c2p(0,0,0)[2] ) ])
        def f(x,y):
            return  ( np.sin(2*x) - np.cos(2*y))/6 + h
        g1 = ParametricFunction(lambda t : [ ran , t , f(ran  , t)  ] , t_range = [ - ran , ran ], color = ORANGE)
        g2 = ParametricFunction(lambda t : [ -ran , -t , f(-ran  , -t)  ] , t_range = [ - ran , ran ], color = ORANGE)
        g3 = ParametricFunction(lambda t : [ -t, ran , f( -t, ran )  ] , t_range = [ - ran , ran ], color = ORANGE)
        g4 = ParametricFunction(lambda t : [ t, -ran , f( t, -ran )  ] , t_range = [ - ran , ran ], color = ORANGE)
        g = VGroup(g1, g2, g3, g4)
        g10 = ParametricFunction(lambda t : [ ran , t ,  nh  ] , t_range = [ - ran , ran ], color = ORANGE)
        g20 = ParametricFunction(lambda t : [ -ran , -t ,  nh  ] , t_range = [ - ran , ran ], color = ORANGE)
        g30 = ParametricFunction(lambda t : [ -t, ran , nh  ] , t_range = [ - ran , ran ], color = ORANGE)
        g40 = ParametricFunction(lambda t : [ t, -ran , nh  ] , t_range = [ - ran , ran ], color = ORANGE)
        si = Surface(
            lambda u, v:  [ u, v, f(u,v) ],
            u_range = [ - 2*ran , 2*ran ],
            v_range = [ - 2*ran , 2*ran ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op*0.75 #, 
            #resolution = [ 16, 16]
        )
        r = Surface(
            lambda u, v:  [ u, v, f(u,v) ],
            u_range = [ - ran  , ran  ],
            v_range = [ - ran  , ran  ],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        si0 = Surface(
            lambda u, v: [ u, v, nh ],
            u_range = [ - 2*ran , 2*ran ],
            v_range = [ - 2*ran , 2*ran ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op*0.75 #, 
            #resolution = [ 16, 16]
        )   
        r0 = Surface(
            lambda u, v: [ u, v, nh ],
            u_range = [ - ran , ran ],
            v_range = [ - ran , ran ],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )   
        self.add_fixed_in_frame_mobjects(t1, t2, s, s1, s2, s1li, s2li, w, s3li, w1)
        self.remove(w, s3li, w1)
        self.wait()
        bye = VGroup(t1, t2, s1, s2, s1li)
        self.play(FadeOut(bye), s.animate.shift(UP), s2li.animate.shift(UP))
        self.wait()
        self.play(Write(w), Create(si), Create(r), Create(axes), Create(s3li), Create(g))
        self.wait()
        self.play(Transform(si, si0), Transform(r, r0), Transform(g1, g10), Transform(g2, g20), Transform(g3, g30), Transform(g4, g40), Write(w1))
        self.wait()


class step4(Scene):
    def construct(self):
        t1 = Tex(r"GB$(R) : = \int_R K + \Psi (\partial R) -2 \pi$", font_size = 34).shift(2.7*UP)
        t2 = Tex(r"Step 1", r" : If $R$ is covered by a semi-geodesic chart, then GB$(R) = 0$ mod $2\pi$.", font_size = 34).next_to(t1, DOWN).shift(0.8*LEFT)
        l = t2.get_left()[0]
        s = Tex(r"Step 2", r" : GB$(R) = 0$ mod $2 \pi$ for all regions as in the theorem.", font_size = 34).next_to(t2, DOWN).shift(UP)
        w = Tex(r"Step 3", r" : If $R$ is contained in the graph of a smooth function, then GB$(R) = 0 .$",   font_size = 34).next_to(s, DOWN)
        w1 = Tex(r"By Step 2, GB$(R)$ doesn't change during the process. Then GB$(R) = 0$.",   font_size = 34).next_to(w, DOWN)
        s4 = Tex(r"Step 4", r" : For all $R$ as in the theorem,  GB$(R) = 0 .$",   font_size = 34).next_to(w1, DOWN)
        s1 = Tex(r"$R = \bigcup_{i = 1}^n R_i$ with $R_i$ in a graph portion of $\Sigma$, and such that ", font_size = 34).next_to(s4, DOWN)
        s2 = Tex(r"$\bigcup _{i = 1} ^ {k+1} R_i $ is obtained by gluing $R_{k+1}$ to $\bigcup _{i=1}^k R_i$ along a curve for each $k$.",  font_size = 34).next_to(s1, DOWN)
        s3 = Tex(r"GB$(R) = \sum_{i=1}^n $GB$(R_i)$ ", r"$= 0 $.",  font_size = 34).next_to(s2, DOWN).shift(0.1*DOWN)
        s.shift((l - s.get_left()[0] )*RIGHT)
        s1.shift((l - s1.get_left()[0] )*RIGHT)
        s2.shift((l - s2.get_left()[0] )*RIGHT)
        s4.shift((l - s4.get_left()[0] )*RIGHT)
        s3.shift(s3.get_center()[0]*LEFT)
        w.shift((l - w.get_left()[0] )*RIGHT)
        w1.shift((l - w1.get_left()[0] )*RIGHT)
        s1li = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(t2[0], DOWN).shift(0.15*UP)
        s2li = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(s[0], DOWN).shift(0.15*UP)
        s3li = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(w[0], DOWN).shift(0.15*UP)
        s4li = Line([0,0,0], [0.85,0,0], color = TEAL).next_to(s4[0], DOWN).shift(0.15*UP)
        box = SurroundingRectangle(s3, buff = .1, color = TEAL)
        self.add(s, w, w1, s2li, s3li)
        self.wait()
        self.play(Write(s4), Create(s4li))
        self.wait()
        self.play(Write(s1), Write(s2))
        self.wait()
        self.play(Write(s3), Create(box))
        self.wait()





class gbf2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/6)
        t = Tex(r"Theorem (Gauss--Bonnet Formula)", font_size = 34, color = BLUE).to_edge(UL).shift(0.7*RIGHT+ 0.2*DOWN)
        l = t.get_left()[0]
        t1 = Tex(r"Let ", r"$\Sigma$", r" be a surface, ", r"$R $", r" $\subset \Sigma $ homeomorphic to a disc, such that ",r"$\partial R$", font_size = 34).next_to(t, DOWN)
        t2 = Tex(r"is a piece-wise smooth curve travelling with $R$ on its left. Then", font_size = 34).next_to(t1, DOWN)
        t3 = Tex(r"$ \int_{R} K +  \Psi (\partial R ) = 2 \pi $.", font_size = 34).next_to(t2, DOWN)
        c = Tex(r"Corollary", font_size = 34, color = BLUE).next_to(t3, DOWN).shift(0.1*UP)
        c1 = Tex(r"Let ", r"$\Delta$", r" be a geodesic triangle in a surface $\Sigma$ with internal angles $\alpha$, $\beta$, $\gamma$.", font_size = 34).next_to(c, DOWN)
        c2 = Tex(r"Then", font_size = 34).next_to(c1, DOWN)
        c3 = Tex(r"$\int_{\Delta} K = \alpha + \beta + \gamma - \pi $.", font_size = 34).next_to(c1, DOWN)
        p = Tex(r"Proof", font_size = 34, color = BLUE).next_to(c3, DOWN).shift(0.1*UP)
        p1 = Tex(r"$\Psi (\Delta ) =  ( \pi - \alpha ) +  (\pi - \beta ) + (\pi - \gamma ) $.", font_size = 34).next_to(p, DOWN)
        t1.shift((l - t1.get_left()[0])*RIGHT)
        t2.shift((l - t2.get_left()[0])*RIGHT)
        c.shift((l - c.get_left()[0])*RIGHT)
        c1.shift((l - c1.get_left()[0])*RIGHT)
        c2.shift((l - c2.get_left()[0])*RIGHT)
        p.shift((l - p.get_left()[0])*RIGHT)
        p1.shift((l - p1.get_left()[0])*RIGHT)
        t3.shift(t3.get_center()[0]*LEFT)
        c3.shift(c3.get_center()[0]*LEFT)
        dli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(c1[1], DOWN).shift(0.15*UP)
        r = 1.5
        h = -2.6
        op = 0.2
        def f(u,v):
            return [r*np.cos(u)*np.sin(v), r*np.sin(u)*np.sin(v), r*np.cos(v) + h ]
        def dp(x,y,z):#dot between X and A
            return x/sqrt(2) + z / sqrt(2)
        def yl(x,y,z):
            return sqrt(1 - dp(x,y,z)**2 )
        def w(x,y,z):#normal from X to A
            return [ ( 1/sqrt(2) - x*dp(x,y,z) ) / yl(x,y,z) , - y*dp(x,y,z)  / yl(x,y,z)  ,   ( 1/sqrt(2) - z*dp(x,y,z) ) / yl(x,y,z)   ]
        def xcomp(x,y,z,t): 
            return (1-t) + t*dp(x,y,z)
        def seno(x,y,z,t):
            return sqrt(1 - xcomp(x,y,z,t)**2)
        def geo(x,y,z,t): # geodesic from X to A sends 0 to X, 1 to A
            return [ r*(x*xcomp(x,y,z,t) + seno(x,y,z,t)*w(x,y,z)[0]),r*( y*xcomp(x,y,z,t) + seno(x,y,z,t)*w(x,y,z)[1] ), r*( z*xcomp(x,y,z,t) + seno(x,y,z,t)*w(x,y,z)[2] ) + h ]
        bcsin = sqrt( 1 - ( (3 + 2 * sqrt(3))/8 )**2 )
        btoc = [ ( -1/4 ) / bcsin , ( sqrt(3) - 6) / (16* bcsin ) , ( 6*sqrt(3) - 3 ) / ( 16 * bcsin )  ]
        def bc(t):
            return [  sqrt(2*t - t**2) *btoc[0] ,  (1-t)* sqrt(3)/2  + sqrt(2*t - t**2) *btoc[1]  ,  (1-t)/2  + sqrt(2*t - t**2) *btoc[2]     ]
        bac = ParametricFunction(lambda t :  geo(0, sqrt(3)/2 , 1/2, t)  , t_range = [0, 1], color = ORANGE)
        cac = ParametricFunction(lambda t :  geo(-1/4, sqrt(3)/4 , sqrt(3)/2, t)  , t_range = [0, 1], color = ORANGE)
        bcc = ParametricFunction(lambda t : [ r*bc(t)[0], r*bc(t)[1], r*bc(t)[2] + h ] , t_range = [0, 1 - ( 3 + 2*sqrt(3)) / 8 ], color = ORANGE)
        s = Surface(
            lambda u, v:  f(u,v),
            u_range = [ 0, 2*PI ],
            v_range = [ 0 , PI / 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op*0.75, 
            resolution = [ 32, 16]
        )
        delta = Surface(
            lambda u, v:  geo( bc(u)[0], bc(u)[1], bc(u)[2], v),
            u_range = [ 0.001, 1 - ( 3 + 2*sqrt(3)) / 8 ],
            v_range = [ 0.001 , 1 ],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = 2*op, 
            resolution = [8,8]
        )
        self.add_fixed_in_frame_mobjects( t, t1, t2, t3, c, c1, c2, c3, p, p1, dli)
        self.remove(c, c1, c2, c3, p, p1, dli)
        self.wait()
        self.play(Write(c), Write(c1), Write(c2), Write(c3), Create(dli), Create(s), Create(delta), Create(bac), Create(cac), Create(bcc))
        self.wait()
        self.play(FadeOut(s), FadeOut(delta), FadeOut(bac), FadeOut(cac), FadeOut(bcc), Write(p), Write(p1))
        self.wait()
        corollary = VGroup(c, c1, c2, c3, dli)
        self.play(FadeOut(t), FadeOut(t1), FadeOut(t2), FadeOut(t3), FadeOut(p), FadeOut(p1), corollary.animate.shift((t1.get_center()[1] - c.get_center()[1])*UP ) )
        self.wait()



class gbf3(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/4)
        t = Tex(r"Theorem (Gauss--Bonnet Formula)", font_size = 34, color = BLUE).to_edge(UL).shift(0.7*RIGHT+ 0.2*DOWN)
        l = t.get_left()[0]
        t1 = Tex(r"Let ", r"$\Sigma$", r" be a surface, ", r"$R $", r" $\subset \Sigma $ homeomorphic to a disc, such that ",r"$\partial R$", font_size = 34).next_to(t, DOWN)
        c = Tex(r"Corollary", font_size = 34, color = BLUE)
        c1 = Tex(r"Let ", r"$\Delta$", r" be a geodesic triangle in a surface $\Sigma$ with internal angles $\alpha$, $\beta$, $\gamma$.", font_size = 34).next_to(c, DOWN)
        c2 = Tex(r"Then", font_size = 34).next_to(c1, DOWN)
        c3 = Tex(r"$\int_{\Delta} K = \alpha + \beta + \gamma - \pi $.", font_size = 34).next_to(c1, DOWN)
        t1.shift((l - t1.get_left()[0])*RIGHT)
        c.shift((l - c.get_left()[0])*RIGHT)
        c1.shift((l - c1.get_left()[0])*RIGHT)
        c2.shift((l - c2.get_left()[0])*RIGHT)
        c3.shift(c3.get_center()[0]*LEFT)
        dli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(c1[1], DOWN).shift(0.15*UP)
        corollary = VGroup(c,c1,c2,c3,dli)
        corollary.shift( (t1.get_center()[1] - c.get_center()[1])*UP  )
        plab = Tex(r"$K > 0$", font_size = 30).next_to(c3, DOWN).shift(4.2*LEFT +0.2*DOWN )
        plab2 = Tex(r"$\alpha + \beta + \gamma > \pi $", font_size = 30).next_to(plab, DOWN)
        plab3 = VGroup(plab, plab2)
        zlab = Tex(r"$K = 0$", font_size = 30).next_to(c3, DOWN).shift( 0.2*DOWN )
        zlab2 = Tex(r"$\alpha + \beta + \gamma = \pi $", font_size = 30).next_to(zlab, DOWN)
        zlab3 = VGroup(zlab, zlab2)
        nlab = Tex(r"$K < 0$", font_size = 30).next_to(c3, DOWN).shift(4.2*RIGHT +0.2*DOWN )
        nlab2 = Tex(r"$\alpha + \beta + \gamma < \pi $", font_size = 30).next_to(nlab, DOWN)
        nlab3 = VGroup(nlab, nlab2)
        boxp = SurroundingRectangle(plab3, buff = .1, color = RED)
        boxz = SurroundingRectangle(zlab3, buff = .1, color = PURPLE_B)
        boxn = SurroundingRectangle(nlab3, buff = .1, color = BLUE)
        r = 1.1
        h = -2
        ht = 3
        m = 1.5
        mn = 1.2
        op = 0.2
        q = 0.18
        def f(u,v):
            return [r*np.cos(u)*np.sin(v) + ht , r*np.sin(u)*np.sin(v) - ht , r*np.cos(v) + h  ]
        def fn(u,v):
            return [ mn * u - ht*(1.05) , mn * v + ht*(1.05) ,  ( v**2 - u **2 ) * 0.7 + h + 0.1 ]
        gpa = ParametricFunction(lambda t :  f( t , PI /2 )  , t_range = [0, PI/2], color = ORANGE)
        gpb = ParametricFunction(lambda t :  f( 0, t )  , t_range = [0, PI /2 ], color = ORANGE)
        gpc = ParametricFunction(lambda t : f( PI/2 , PI/2 -  t )  , t_range = [0, PI/2 ], color = ORANGE)
        gza = ParametricFunction(lambda t :  [ m*(t * ( 0.7 + 0.5 )  - 0.5) , m*( (1 - t ) * (-0.8)) , h  ]  , t_range = [0,1], color = ORANGE)
        gzb = ParametricFunction(lambda t :  [ m*((1 - t) * ( 0.7 + 0.5 )  - 0.5) , m*( t * (0.8)) , h  ]  , t_range = [0, 1 ], color = ORANGE)
        gzc = ParametricFunction(lambda t :  [ m*(-0.5) ,m*( -t) , h  ] , t_range = [ -0.8 , 0.8 ], color = ORANGE)
        gna = ParametricFunction(lambda t :  fn( t * ( 0.7 + 0.5 )  - 0.5 ,  (1 - t - q* np.sin(PI * t) ) * (-0.8))  , t_range = [0,1], color = ORANGE)
        gnb = ParametricFunction(lambda t :  fn( (1 - t) * ( 0.7 + 0.5 )  - 0.5 ,  (t - q* np.sin(PI * t) )  * 0.8)  , t_range = [0, 1 ], color = ORANGE)
        gnc = ParametricFunction(lambda t :  fn( -0.5 + q * np.cos(PI*t / 1.6) , - t    ) , t_range = [ -0.8 , 0.8 ], color = ORANGE)
        sp = Surface(
            lambda u, v:  f(u,v),
            u_range = [ 0, 2*PI ],
            v_range = [ 0 , 3* PI / 4 ],
            checkerboard_colors = [RED, RED],
            fill_opacity = op*0.75, 
            resolution = [ 32, 24]
        )
        deltap = Surface(
            lambda u, v:  f(u,v),
            u_range = [ 0, PI/2 ],
            v_range = [ 0, PI /2 ],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = 2*op, 
            resolution = [8,8]
        )
        sz = Surface(
            lambda u, v:  [m*u,m*v, h  ],
            u_range = [ -1 , 1  ],
            v_range = [ -1 , 1  ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op*0.75, 
            resolution = [ 16, 16]
        )
        deltaz = Surface(
            lambda u, v:  [ m*( v * ( 0.7 + 0.5 ) - 0.5) , m*(1 - v  ) * u  , h  ],
            u_range = [ -0.8 , 0.8 ],
            v_range = [ 0 , 1  ],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = 2*op, 
            resolution = [8,8]
        )
        sn = Surface(
            lambda u, v:  fn(u,v),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = op*0.75, 
            resolution = [ 16, 16]
        )
        deltan = Surface(
            lambda u, v:  fn( v * ( 0.7 + (0.5 - q * np.cos(PI*u / 1.6)) ) - (0.5 - q * np.cos(PI*u / 1.6)) , (1 - v - q* np.sin(PI * v) ) * u ),
            u_range = [ -0.8 , 0.8 ],
            v_range = [ 0 , 1  ],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = 2*op, 
            resolution = [8,8]
        )
        self.add_fixed_in_frame_mobjects( c, c1, c2, c3, dli, plab, plab2, zlab, zlab2, nlab, nlab2, boxp, boxz, boxn)
        self.remove(plab, plab2, zlab, zlab2, nlab, nlab2, boxp, boxz, boxn)
        self.wait()
        self.play(
            Create(plab3), Create(zlab3), Create(nlab3), 
            Create(boxp), Create(boxz), Create(boxn), 
            Create(sp), Create(deltap), Create(gpa), Create(gpb), Create(gpc), 
            Create(sz), Create(deltaz), Create(gza), Create(gzb), Create(gzc), 
            Create(sn), Create(deltan), Create(gna), Create(gnb), Create(gnc)
        )
        self.wait()




class gbt(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70 * DEGREES, theta=  PI/4)
        t = Tex(r"Theorem (Gauss--Bonnet)", font_size = 34, color = BLUE).to_edge(UL).shift(0.7*RIGHT+ 0.5*DOWN)
        l = t.get_left()[0]
        t1 = Tex(r"Let ", r"$\Sigma$", r" be a closed surface. Then", font_size = 34).next_to(t, DOWN)
        t2 = MathTex("\int_{\Sigma} K = 2 \pi \chi (\Sigma ).", font_size = 34).next_to(t1, DOWN)
        ec = Tex(r"$\chi (\Sigma) =  $ ", r"$V$", r" $ - $ ", r"$E$", r" $+$ ", r"$F$", font_size = 34).to_edge(UR).shift( 1.05*  DOWN + LEFT)
        p = Tex(r"Proof", font_size = 34, color = BLUE).next_to(t2, DOWN)
        p1 = Tex(r"Let ", r"$\Delta$", r" be a geodesic triangle in a surface $\Sigma$ with internal angles $\alpha$, $\beta$, $\gamma$.", font_size = 34).next_to(p, DOWN)
        t1.shift((l - t1.get_left()[0])*RIGHT)
        t2.shift(t2.get_center()[0]*LEFT)
        p.shift((l - p.get_left()[0])*RIGHT)
        p1.shift((l - p1.get_left()[0])*RIGHT)
        vli = Line([0,0,0], [0.4, 0, 0], color = BLUE).next_to(ec[1], DOWN).shift(0.15*UP)
        eli = Line([0,0,0], [0.4, 0, 0], color = ORANGE).next_to(ec[3], DOWN).shift(0.15*UP)
        fli = Line([0,0,0], [0.4, 0, 0], color = PURPLE_B).next_to(ec[5], DOWN).shift(0.15*UP)
        r = 1.35
        rs = 0.06
        h = -1.1
        ht = 2.1
        r0 = 0.5
        r1 = 1.6
        op = 0.3
        def sphere(u,v):
            return [ r*np.cos(u)*np.sin(v) + ht , r*np.sin(u)*np.sin(v) - ht  , r*np.cos(v) + h   ]
        s1 = Surface(
            lambda u, v: sphere(u,v) ,
            u_range = [ 0, 2*PI ],
            v_range = [ 0 , PI  ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op*0.75, 
            resolution = [ 16, 16]
        )
        def torus(u,v):
            return [  ( r1 + r0*np.cos(v)) * np.cos(u)    -   ht  , ( r1 + r0*np.cos(v)) * np.sin(u)  + ht  ,  r0 * np.sin(v)  +  h  ]
        s2 = Surface(
            lambda u, v:  torus(u,v),
            u_range = [ 0, 2*PI  ],
            v_range = [ 0, 2*PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op*0.75, 
            resolution = [ 16, 16]
        )
        v1 = []
        v2 = []
        e1 = []
        e2 = []
        ver = VGroup()
        edges = VGroup()
        e1.append(
            ParametricFunction(lambda t :  sphere( t , PI /2 )  , t_range = [0, 2*PI ], color = ORANGE).set_shade_in_3d(True)
        )
        for i in range(0,4):
            v1.append(Sphere(radius = rs))
            v1[i].shift( sphere( PI * i / 2 , PI/2 ) ).set_color(BLUE)
            e1.append( 
                ParametricFunction(lambda t :  sphere( PI * i / 2 , t )  , t_range = [0, PI ], color = ORANGE).set_shade_in_3d(True)
             )
            edges.add(e1[i])
            e2.append(
                ParametricFunction(lambda t :  torus( PI * i / 2 , t )  , t_range = [0, 2*PI ], color = ORANGE).set_shade_in_3d(True)
            )
            e2.append(
                ParametricFunction(lambda t :  torus( t,  PI * i / 2 )  , t_range = [0, 2*PI ], color = ORANGE).set_shade_in_3d(True)
            )
            e2.append(
                ParametricFunction(lambda t :  torus( t,  PI * i / 2 + t )  , t_range = [0, 2*PI ], color = ORANGE).set_shade_in_3d(True)
            )
        edges.add(e1[4])
        for i in range(0,12):
            edges.add(e2[i])
        v1.append(Sphere(radius = rs, resolution = [8,8]))
        v1.append(Sphere(radius = rs, resolution = [8,8]))
        v1[4].shift(sphere(0,0)).set_color(BLUE)
        v1[5].shift(sphere(0, PI)).set_color(BLUE)
        for i in range(0,6):
            ver.add(v1[i])
        for i in range(0, 4):
            v2.append([])
            for j in range(0,4):
                v2[i].append(Sphere(radius = rs, resolution = [8,8]))
                v2[i][j].shift(  torus( PI * i /2  , PI * j /2  ) ).set_color(BLUE)
                ver.add(v2[i][j])
        self.add_fixed_in_frame_mobjects( t, t1, t2, p, p1, ec, vli, eli, fli )
        self.remove( t, t1, t2, p, p1, ec, vli, eli, fli )
        self.play(Create(t), Create(t1), Create(t2), Create(s1) , Create(s2))
        self.wait()
        self.play(Create(edges), Create(ver), Write(ec), Create(vli), Create(eli), Create(fli))
        self.wait()



class gbt2(Scene):
    def construct(self):
        t = Tex(r"Theorem (Gauss--Bonnet)", font_size = 34, color = BLUE).to_edge(UL).shift(0.7*RIGHT+ 0.5*DOWN)
        l = t.get_left()[0]
        t1 = Tex(r"Let ", r"$\Sigma$", r" be a closed surface. Then", font_size = 34).next_to(t, DOWN)
        t2 = MathTex("\int_{\Sigma} K = 2 \pi \chi (\Sigma ).", font_size = 34).next_to(t1, DOWN)
        ec = Tex(r"$\chi (\Sigma) =  $ ", r"$V$", r" $ - $ ", r"$E$", r" $+$ ", r"$F$", font_size = 34).to_edge(UR).shift( 1.05* DOWN + LEFT)
        p = Tex(r"Proof", font_size = 34, color = BLUE).next_to(t2, DOWN)
        p1 = Tex(r"$\Sigma  = \bigcup _{i =1} ^F R_i$. Apply The Gauss-Bonnet Formula to each triangle to obtain", font_size = 34).next_to(p, DOWN)
        p2 = Tex(r"$\sum_{i = 1} ^F    \int_{R_i} K $" , r" $+$ ", r"$\sum_{i = 1}^F \Psi (\partial R_i ) $", r" $=$ ", r"$2 \pi F$", font_size = 34).next_to(p1, DOWN)
        p22 = Tex(r"$  \sum_{i = 1} ^F    \int_{R_i} K  $" , r" $=$ ", r"$2 \pi F$", r" $-$ ", r"$\sum_{i = 1}^F \Psi (\partial R_i )$" , font_size = 34).next_to(p1, DOWN)
        t1.shift((l - t1.get_left()[0])*RIGHT)
        t2.shift(t2.get_center()[0]*LEFT)
        p.shift((l - p.get_left()[0])*RIGHT)
        p1.shift((l - p1.get_left()[0])*RIGHT)
        p2.shift(p2.get_center()[0]*LEFT)
        p22.shift(p22.get_center()[0]*LEFT)
        p3 = Tex(r"$\int_{\Sigma} K = $ ", font_size = 34).next_to(p22, LEFT)
        eq1 = VGroup(p3, p22)
        cen = eq1.get_center()
        p3.shift(cen*LEFT)
        p4 = Tex(r"Claim:", r" $\sum _{i=1}^F \Psi (\partial R_i) = 2\pi (E - V )$.", font_size = 34).next_to(p22, DOWN)
        p5 = Tex(r"If the claim is true, then $\int_{\Sigma} K = 2 \pi (F - E + V) $ ", r"$ = 2 \pi \chi (\Sigma) $.", font_size = 34).next_to(p4, DOWN)
        p4.shift((l - p4.get_left()[0])*RIGHT)
        p5.shift((l - p5.get_left()[0])*RIGHT)
        vli = Line([0,0,0], [0.4, 0, 0], color = BLUE).next_to(ec[1], DOWN).shift(0.15*UP)
        eli = Line([0,0,0], [0.4, 0, 0], color = ORANGE).next_to(ec[3], DOWN).shift(0.15*UP)
        fli = Line([0,0,0], [0.4, 0, 0], color = PURPLE_B).next_to(ec[5], DOWN).shift(0.15*UP)
        gcli = Line([0,0,0], [1.8, 0, 0], color = TEAL).next_to(p22[4], DOWN).shift(0.15*UP + cen*LEFT )
        gcli2 = Line([0,0,0], [1, 0, 0], color = TEAL).next_to(p4[0], DOWN).shift(0.15*UP)
        fli2 = Line([0,0,0], [0.4, 0, 0], color = PURPLE_B).next_to(p2[4], DOWN).shift(0.15*UP + 0.2*RIGHT)
        fli3 = Line([0,0,0], [0.4, 0, 0], color = PURPLE_B).next_to(p22[2], DOWN).shift(0.15*UP + 0.2*RIGHT)
        fli4 = Line([0,0,0], [0.4, 0, 0], color = PURPLE_B).next_to(p5, DOWN).shift(0.2*UP + (1.2)*RIGHT)
        eli2 = Line([0,0,0], [0.4, 0, 0], color = ORANGE).next_to(p4, DOWN).shift(0.2*UP + 1.65*RIGHT)
        eli3 = Line([0,0,0], [0.4, 0, 0], color = ORANGE).next_to(p5, DOWN).shift(0.2*UP + 1.9*RIGHT)
        vli2 = Line([0,0,0], [0.4, 0, 0], color = BLUE).next_to(p4, DOWN).shift(0.2*UP + 2.35*RIGHT)
        vli3 = Line([0,0,0], [0.4, 0, 0], color = BLUE).next_to(p5, DOWN).shift(0.2*UP + 2.62*RIGHT)
        self.add(t, t1, t2, ec, vli, eli, fli)
        self.wait()
        self.play(Write(p), Write(p1), Write(p2), Create(fli2))
        self.wait()
        self.play(Transform(p2[0], p22[0]), FadeOut(p2[1]), Transform(p2[2], p22[4]), Transform(p2[3], p22[1]), Transform(p2[4], p22[2]), Write(p22[3]), Transform(fli2, fli3))
        self.wait()
        self.play(Write(p3), p22[0].animate.shift(cen*LEFT), p22[4].animate.shift(cen*LEFT), p22[1].animate.shift(cen*LEFT), p22[2].animate.shift(cen*LEFT), p22[3].animate.shift(cen*LEFT), p2[0].animate.shift(cen*LEFT), p2[2].animate.shift(cen*LEFT), p2[3].animate.shift(cen*LEFT), p2[4].animate.shift(cen*LEFT), fli2.animate.shift(cen*LEFT), fli3.animate.shift(cen*LEFT))
        self.wait()
        self.play(Create(gcli))
        self.wait()
        self.play(Write(p4), Create(gcli2), Create(eli2), Create(vli2))
        self.wait()
        self.play(Write(p5), Create(fli4), Create(eli3), Create(vli3))
        self.wait()
        sur = VGroup(p4, eli2, vli2, gcli2)
        self.play(FadeOut(p), FadeOut(t1), FadeOut(t2), FadeOut(t), FadeOut(p1), FadeOut(p3), FadeOut(p5), FadeOut(p22[0]), FadeOut(p22[4]), FadeOut(p22[1]), FadeOut(p22[2]), FadeOut(p22[3]), FadeOut(p2[0]), FadeOut(p2[2]), FadeOut(p2[3]), FadeOut(p2[4]), sur.animate.shift( ( ec.get_center()[1] - p4.get_center()[1] - 0.4 )*UP), FadeOut(gcli), FadeOut(fli2), FadeOut(fli3), FadeOut(vli3), FadeOut(eli3), FadeOut(fli4))
        self.wait()




class gbt3(Scene):
    def construct(self):
        slide = 0.3
        ec = Tex(r"$\chi (\Sigma) =  $ ", r"$V$", r" $ - $ ", r"$E$", r" $+$ ", r"$F$", font_size = 34).to_edge(UR).shift( 1.05* DOWN + LEFT)
        p4 = Tex(r"Claim:", r" $\sum _{i=1}^F \Psi (\partial R_i) = 2\pi (E - V )$.", font_size = 34).to_edge(UL)
        p4.shift(0.7*RIGHT + ( ec.get_center()[1] - p4.get_center()[1] - slide - 0.1 )*UP )
        l = p4.get_left()[0]
        p5 = Tex(r"$\sum _{i=1}^F \Psi (\partial R_i) = $ ", r"Edge contributions", r" $+$ Vertex contributions", font_size = 34).next_to(p4, DOWN)
        p5.shift((l - p5.get_left()[0])*RIGHT + slide * DOWN )
        p6 = Tex(r"Contribution at ", r"a vertex", r" $ = $", r" $ \sum_j ( \pi - \theta_j )$", font_size = 34).next_to(p5, DOWN)
        p6.shift((l - p6.get_left()[0])*RIGHT + slide * DOWN )
        p7 = Tex(r"$ =  \pi \cdot \# $(faces adjacent to the point) $ - 2 \pi $", font_size = 34).next_to(p6, DOWN)
        p8 = Tex(r"$ \sum_{i = 1 } ^ F \Psi ( \partial R_i )  = $ ", r"$ 3 \pi F$", r" $ - 2 \pi V$", r" $ = 2 \pi E - 2 \pi V$." , font_size = 34).next_to(p7, DOWN)
        p7.shift((l - p7.get_left()[0])*RIGHT + 0.1*UP )
        p8.shift((l - p8.get_left()[0])*RIGHT + slide * DOWN )
        vli = Line([0,0,0], [0.4, 0, 0], color = BLUE).next_to(ec[1], DOWN).shift(0.15*UP)
        eli = Line([0,0,0], [0.4, 0, 0], color = ORANGE).next_to(ec[3], DOWN).shift(0.15*UP)
        fli = Line([0,0,0], [0.4, 0, 0], color = PURPLE_B).next_to(ec[5], DOWN).shift(0.15*UP)
        gcli2 = Line([0,0,0], [1, 0, 0], color = TEAL).next_to(p4[0], DOWN).shift(0.15*UP)
        eli2 = Line([0,0,0], [0.4, 0, 0], color = ORANGE).next_to(p4, DOWN).shift(0.2*UP + 1.65*RIGHT)
        vli2 = Line([0,0,0], [0.4, 0, 0], color = BLUE).next_to(p4, DOWN).shift(0.2*UP + 2.35*RIGHT)
        nedli = Line([0,0,0], [3, 0.35, 0], color = RED).next_to(p5[1], DOWN).shift( 0.6 * UP )
        vli3 = Line([0,0,0], [1.3, 0, 0], color = BLUE).next_to(p6[1], DOWN).shift(0.16*UP)
        fli2 = Line([0,0,0], [0.4, 0, 0], color = PURPLE_B).next_to(p8[1], DOWN).shift(0.16*UP + 0.2*RIGHT )
        vli4 = Line([0,0,0], [0.4, 0, 0], color = BLUE).next_to(p8[2], DOWN).shift(0.16*UP + 0.3*RIGHT )
        eli3 = fli2.copy().set_color(ORANGE).shift(2.2*RIGHT)
        vli5 = vli4.copy().shift(2.3*RIGHT)
        op = 0.4
        a = PI/10
        f1 = Triangle(color = PURPLE_B, fill_opacity = op)
        f2 = f1.copy()
        v1 = Dot( [ f1.get_left()[0], f1.get_bottom()[1] , 0 ] ,   radius = 0.1 , color = BLUE )
        v2 = Dot( [ f1.get_right()[0], f1.get_bottom()[1] , 0 ] ,   radius = 0.1 , color = BLUE )
        v3 = Dot( [ f1.get_top()[0], f1.get_top()[1] , 0 ] ,   radius = 0.1 , color = BLUE )
        v4 = Dot( [ f1.get_top()[0], 2*f1.get_bottom()[1] - f1.get_top()[1] , 0 ] ,   radius = 0.1 , color = BLUE )
        e1 = Line( v1.get_center() , v2.get_center() , color = ORANGE )
        e2 = Line( v2.get_center() , v3.get_center() , color = ORANGE )
        e3 = Line( v3.get_center() , v1.get_center() , color = ORANGE )
        e4 = Line( v4.get_center() , v2.get_center() , color = ORANGE )
        e5 = Line( v1.get_center() , v4.get_center() , color = ORANGE )
        cen = [ f1.get_left()[0], f1.get_bottom()[1] , 0 ]
        fig1 = VGroup(f1, f2, e1, e2, e3, e4, e5, v1, v2, v3, v4)
        f2.rotate(angle =  - PI / 3 , about_point = cen )
        fig1.rotate(angle =  PI / 6 + a , about_point = cen ).shift( 1.6 * DOWN + 0.5 * RIGHT )
        sho = 0.1
        lon = 0.8
        a1 = Arrow( [ (1 + sho) * v3.get_center()[0] - sho * v1.get_center()[0] , (1 + sho) * v3.get_center()[1] - sho * v1.get_center()[1] , 0 ] , [ (1 - lon) * v3.get_center()[0] + lon * v1.get_center()[0] , (1 - lon) * v3.get_center()[1] + lon * v1.get_center()[1] , 0 ] , color = YELLOW , max_tip_length_to_length_ratio=0.2)
        a2 = Arrow( [ (1 + sho) * v4.get_center()[0] - sho * v2.get_center()[0] , (1 + sho) * v4.get_center()[1] - sho * v2.get_center()[1] , 0 ] , [ (1 - lon) * v4.get_center()[0] + lon * v2.get_center()[0] , (1 - lon) * v4.get_center()[1] + lon * v2.get_center()[1] , 0 ] , color = YELLOW , max_tip_length_to_length_ratio=0.2)
        d1 = [ v2.get_center()[0] - v1.get_center()[0] , v2.get_center()[1] - v1.get_center()[1] , 0 ]
        d2 = [ v3.get_center()[0] - v2.get_center()[0] , v3.get_center()[1] - v2.get_center()[1] , 0 ]
        d3 = [ v1.get_center()[0] - v3.get_center()[0] , v1.get_center()[1] - v3.get_center()[1] , 0 ]
        d4 = [ v2.get_center()[0] - v4.get_center()[0] , v2.get_center()[1] - v4.get_center()[1] , 0 ]
        d5 = [ v4.get_center()[0] - v1.get_center()[0] , v4.get_center()[1] - v1.get_center()[1] , 0 ]
        fig2 = VGroup()
        fig3 = VGroup()
        leng  = 1.6
        m = 6
        ale = 0.35 ### arrow length = 1 - ale
        ee = []
        ff = []
        aa = []
        tt = []
        vv = []
        th = []
        dd = []
        ff.append( Triangle( color = PURPLE_B, fill_opacity = op ))
        ff[0].shift( [ - ff[0].get_left()[0] , - ff[0].get_bottom()[1] , 0 ] )
        fig2.add(ff[0])
        for i in range(0,m):
            ee.append( Line( [0,0,0] , [leng,0,0] , color = ORANGE ))
            ee[i].rotate( 2 * PI * i / m, about_point = [0,0,0]  )
            vv.append( Dot ( [ 1.2* leng , 0 , 0] ) )
            vv[i].rotate( 2 * PI * i / m, about_point = [0,0,0]  )
            dd.append(  vv[i].get_center()  )
            aa.append( Arrow( [ 1.05* dd[i][0] , 1.05* dd[i][1] , 0 ] , [ ale * dd[i][0] , ale * dd[i][1] , 0 ] , color = YELLOW , max_tip_length_to_length_ratio=0.2 ) )
            if i > 0:
                ff.append(ff[0].copy())
                ff[i].rotate( 2 * PI * i / m , about_point = [0,0,0])
                fig2.add(ff[i])
                tt.append( Angle(ee[i-1], ee[i], radius = 0.5 + 0.2 * np.sin( 3 * i )) )
                fig2.add(tt[i-1])
        aa[0].set_color(MAROON_E)
        aa[1].set_color(RED)
        aa[3].set_color(GREEN)
        aa[4].set_color(BLUE_E)
        aa[5].set_color(PINK)
        tt.append( Angle( ee[5], ee[0] , radius = 0.4 ) )
        fig2.add(tt[5])
        th.append( Tex(r"$\theta_1$", font_size = 25 ).next_to(tt[0], RIGHT).shift( 0.1 * LEFT + 0.1 *UP ) )
        th.append( Tex(r"$\theta_2$", font_size = 25 ).next_to(tt[1], UP).shift( 0.1 * DOWN ))
        th.append( Tex(r"$\theta_3$", font_size = 25 ).next_to(tt[2], LEFT).shift( 0.1 * RIGHT + 0.1 * UP ))
        th.append( Tex(r"$\theta_4$", font_size = 25 ).next_to(tt[3], LEFT).shift( 0.1 * RIGHT + 0.1 * DOWN ))
        th.append( Tex(r"$\theta_5$", font_size = 25 ).next_to(tt[4], DOWN).shift( 0.1 * UP ))
        th.append( Tex(r"$\theta_6$", font_size = 25 ).next_to(tt[5], RIGHT).shift( 0.1 * LEFT + 0.1 * DOWN ))
        for i in range(0,m):
            fig2.add(th[i], ee[i])
            fig3.add(aa[i])
        vv.append( Dot( radius = 0.1 , color = BLUE ) )
        fig2.add(vv[6])
        fig2.shift( 1.5*DOWN + 4 * RIGHT )
        fig3.shift( 1.5*DOWN + 4 * RIGHT )
        self.add( ec, vli, eli, fli, p4, gcli2, eli2, vli2)
        self.wait()
        self.play( Create(fig1) , Create(p5))
        self.wait()
        self.play(Create(a1), Create(a2))
        self.wait()
        self.play( a1.animate.shift( d3 ) , a2.animate.shift( d4 ) )
        self.play(  Rotate(a1, angle = 2* PI /3, about_point = v1.get_center()) , Rotate(a2, angle = 2* PI /3, about_point = v2.get_center()) )
        self.play( a1.animate.shift( d1 ) , a2.animate.shift( [ - d1[0] , - d1[1], 0 ] ) )
        self.play(  Rotate(a1, angle = 2* PI /3, about_point = v2.get_center()) , Rotate(a2, angle = 2* PI /3, about_point = v1.get_center()) )
        self.play( a1.animate.shift( d2 ) , a2.animate.shift( d5 ) , Create(nedli) )
        self.wait()
        self.play(FadeOut(fig1), FadeOut(a1), FadeOut(a2))
        self.wait()
        self.play(Write(p6[0]), Write(p6[1]), Create(fig2), Create(vli3), Write(p6[2]))
        self.wait()
        self.play(Create(fig3))
        self.wait()
        self.play(
            aa[0].animate.shift( [ - dd[0][0] , - dd[0][1], 0 ]  ), 
            aa[1].animate.shift( [ - dd[1][0] , - dd[1][1], 0 ]  ), 
            aa[2].animate.shift( [ - dd[2][0] , - dd[2][1], 0 ]  ), 
            aa[3].animate.shift( [ - dd[3][0] , - dd[3][1], 0 ]  ), 
            aa[4].animate.shift( [ - dd[4][0] , - dd[4][1], 0 ]  ), 
            aa[5].animate.shift( [ - dd[5][0] , - dd[5][1], 0 ]  ) 
        )
        self.play(
            Rotate( aa[0] , angle = 2 * PI / 3 , about_point = vv[6].get_center() ), 
            Rotate( aa[1] , angle = 2 * PI / 3 , about_point = vv[6].get_center() ), 
            Rotate( aa[2] , angle = 2 * PI / 3 , about_point = vv[6].get_center() ), 
            Rotate( aa[3] , angle = 2 * PI / 3 , about_point = vv[6].get_center() ), 
            Rotate( aa[4] , angle = 2 * PI / 3 , about_point = vv[6].get_center() ), 
            Rotate( aa[5] , angle = 2 * PI / 3 , about_point = vv[6].get_center() ), 
            run_time = 2
        )
        self.play(
            aa[0].animate.shift( [ 0.5* dd[5][0] , 0.5*dd[5][1] , 0 ] ), 
            aa[1].animate.shift( [ 0.5* dd[0][0] , 0.5*dd[0][1] , 0 ] ), 
            aa[2].animate.shift( [ 0.5* dd[1][0] , 0.5*dd[1][1] , 0 ] ), 
            aa[3].animate.shift( [ 0.5* dd[2][0] , 0.5*dd[2][1] , 0 ] ), 
            aa[4].animate.shift( [ 0.5* dd[3][0] , 0.5*dd[3][1] , 0 ] ), 
            aa[5].animate.shift( [ 0.5* dd[4][0] , 0.5*dd[4][1] , 0 ] ) 
        )
        self.wait()
        self.play(Write(p6[3]))
        self.wait()
        self.play(Write(p7))
        self.wait()
        self.play(Write(p8[0]), Write(p8[1]), Write(p8[2]), Create(fli2), Create(vli4))
        self.wait()
        self.play(Write(p8[3]), Create(eli3), Create(vli5))
        self.wait()







