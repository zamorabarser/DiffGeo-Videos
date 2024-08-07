from this import d
from tkinter import E
from manim import *
from numpy import sqrt
import math



class ti(Scene):
    def construct(self):
        t1 = Text("Curvature of surfaces", font_size=60)
        self.play(Write(t1))
        self.wait()        


class normal0(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 45*DEGREES)
        sl = Tex(r"For a point ", r"$p$", r" in a surface ", r"$\Sigma$", r" there are two unit normal vectors.", font_size = 34).shift(2.15*UP)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE_B).next_to(sl[3], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(sl[1], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 5, y_length = 5, z_length = 4)
        def f(u,v):
            return 2**(2 - (u**2+v**2)/8)
        s = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v) - 3 ),
            u_range = [-4,4],
            v_range = [-4,4],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        a = math.log(2) #the following has length 3
        def nor1(u,v):
            return 3*a*u*f(u,v)/(4*sqrt(1 + ((a**2)*(u**2 + v**2)*(f(u,v))**2)/16 ))
        def nor2(u,v):
            return 3*a*v*f(u,v)/(4*sqrt(1 + ((a**2)*(u**2 + v**2)*(f(u,v))**2)/16 ))
        def nor3(u,v):
            return 3/sqrt(1 + ((a**2)*(u**2 + v**2)*(f(u,v))**2)/16 )
        n = Arrow( ax.c2p(0 - 0.1 * nor1(0,sqrt(2)) ,sqrt(2) - 0.1 * nor2(0,sqrt(2)) , 2**(7/4) - 3 - 0.1 * nor3(0,sqrt(2))) , ax.c2p( 0 + nor1(0,sqrt(2)), sqrt(2) + nor2(0,sqrt(2)), -3 + 2**(7/4) + nor3(0,sqrt(2))  ), color = RED)   
        n2 = Arrow( ax.c2p(0 + 0.1 * nor1(0,sqrt(2)) ,sqrt(2) + 0.1 * nor2(0,sqrt(2)) , 2**(7/4) - 3 + 0.1 * nor3(0,sqrt(2))) , ax.c2p( 0 - nor1(0,sqrt(2)), sqrt(2) - nor2(0,sqrt(2)), -3 + 2**(7/4) - nor3(0,sqrt(2))  ), color = GREEN)   
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0, sqrt(2), 2**(7/4)-3))
        self.add_fixed_in_frame_mobjects(sl, sli, pli)
        self.remove(sl, sli, pli)
        self.play(Create(s), Create(p), Write(sl), Create(sli), Create(pli))
        self.wait()
        self.play(Create(n), Create(n2))
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/4)
        self.wait(2)
        self.stop_ambient_camera_rotation()
        self.wait()
        #self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(sl), FadeOut(n2), FadeOut(sli), FadeOut(pli))
        


class normal1(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 135*DEGREES)
        sl = Tex(r"For a point ", r"$p$", r" in a surface ", r"$\Sigma$", r" there are two unit normal vectors.", font_size = 34).shift(2.15*UP)
        nl = Tex(r"If ", r"$\Sigma$", r" has two sides, there is ", r"$N$", r" $: $ ", r"$\Sigma$", r" $ \to \mathbb{S}^2$ normal vector.", font_size = 34).shift(2.15*UP)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE_B).next_to(sl[3], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(sl[1], DOWN).shift(0.15*UP)
        nli = Line([0,0,0], [0.3,0,0], color = RED).next_to(nl[3], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.4,0,0], color = PURPLE_B).next_to(nl[1], DOWN).shift(0.15*UP)
        sli3 = Line([0,0,0], [0.4,0,0], color = PURPLE_B).next_to(nl[5], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,6,10], x_length = 5, y_length = 5, z_length = 4)
        def f(u,v):
            return 2**(2 - (u**2+v**2)/8)
        s = Surface(
            lambda u, v: ax.c2p( u , v , f(u,v) - 3 ),
            u_range = [-4,4],
            v_range = [-4,4],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        a = math.log(2) #the following has length 3
        def nor1(u,v):
            return 3*a*u*f(u,v)/(4*sqrt(1 + ((a**2)*(u**2 + v**2)*(f(u,v))**2)/16 ))
        def nor2(u,v):
            return 3*a*v*f(u,v)/(4*sqrt(1 + ((a**2)*(u**2 + v**2)*(f(u,v))**2)/16 ))
        def nor3(u,v): 
            return 3/sqrt(1 + ((a**2)*(u**2 + v**2)*(f(u,v))**2)/16 )
        n = Arrow( ax.c2p(0 - 0.1 * nor1(0,sqrt(2)) ,sqrt(2) - 0.1 * nor2(0,sqrt(2)) , 2**(7/4) - 3 - 0.1 * nor3(0,sqrt(2))) , ax.c2p( 0 + nor1(0,sqrt(2)), sqrt(2) + nor2(0,sqrt(2)), -3 + 2**(7/4) + nor3(0,sqrt(2))  ), color = RED)   
        n2 = Arrow( ax.c2p(0 + 0.1 * nor1(0,sqrt(2)) ,sqrt(2) + 0.1 * nor2(0,sqrt(2)) , 2**(7/4) - 3 + 0.1 * nor3(0,sqrt(2))) , ax.c2p( 0 - nor1(0,sqrt(2)), sqrt(2) - nor2(0,sqrt(2)), -3 + 2**(7/4) - nor3(0,sqrt(2))  ), color = GREEN)   
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0, sqrt(2), 2**(7/4)-3))
        self.add_fixed_in_frame_mobjects(sl, sli, pli, nl, sli2, sli3, nli)
        self.add(p,n,n2,s)
        self.remove(sli2, sli3, nli, nl)
        self.wait()
        self.play(FadeOut(n2), FadeOut(sl), FadeOut(sli), FadeOut(pli))
        self.wait()
        self.play(Write(nl), Create(sli2), Create(nli), Create(sli3))
        #self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(sl), FadeOut(pl), FadeOut(n2), FadeOut(sli), FadeOut(pli))
        z = ValueTracker(0)
        def xco(t):
            return np.sin(t)*(4 - np.cos(t/2))*sqrt(2)/3
        def yco(t):
            return np.cos(t)*(4 -  np.cos(t/2))*sqrt(2)/3
        n.add_updater( lambda x: x.become( Arrow( ax.c2p(  xco(z.get_value()) - (0.1)*nor1( xco(z.get_value()), yco(z.get_value()) ) , yco(z.get_value()) - (0.1)*nor2( xco(z.get_value()), yco(z.get_value()) ) , -3 + f(xco(z.get_value()), yco(z.get_value())) - (0.1)*nor3( xco(z.get_value()), yco(z.get_value()) )),ax.c2p(  xco(z.get_value()) + nor1( xco(z.get_value()), yco(z.get_value()) ) , yco(z.get_value()) + nor2( xco(z.get_value()), yco(z.get_value()) ) , -3 + f(xco(z.get_value()), yco(z.get_value())) + nor3( xco(z.get_value()), yco(z.get_value()) )) , color=RED  )))
        p.add_updater( lambda x: x.become( Sphere(radius = 0.05).set_color(ORANGE).shift(ax.c2p( xco(z.get_value()), yco(z.get_value()) , -3 + f( xco(z.get_value()), yco(z.get_value()))))))
        self.play(z.animate.set_value(2*TAU), run_time = 8)
        self.wait()

class mobnor(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,2,10], x_length = 10, y_length = 10, z_length = 4 )
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        mb = Surface(
            lambda u, v: axes.c2p( np.cos(u)*(3 + np.cos(u/2)*v),np.sin(u)*(3 + np.cos(u/2)*v) , np.sin(u/2)*v ),
            u_range = [0,2*PI],
            v_range = [-1,1],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        mob = Tex(r"$\Sigma  =$  M\"obius band", font_size = 40).shift(3*UP)
        self.add_fixed_in_frame_mobjects(mob)
        self.remove(mob)
        self.play(Create(mb), Write(mob))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/5)
        self.wait(10)
        self.stop_ambient_camera_rotation()
        self.wait()
        n = Arrow(axes.c2p(3, 0, -0.1) , axes.c2p(3,0,1.5), color = RED  )
        p = Sphere(radius = 0.1 ).set_color(ORANGE).shift(axes.c2p(3,0,0))
        self.play(Create(n), Create(p))
        z = ValueTracker(0)
        n.add_updater( lambda x: x.become( Arrow( axes.c2p( 3* np.cos(z.get_value()) + 0.3*np.cos(z.get_value()) * np.sin(z.get_value()/2) , 3 *np.sin(z.get_value()) + 0.3*np.sin(z.get_value())*np.sin(z.get_value()/2) , - 0.3*np.cos(z.get_value()/2) )  , axes.c2p( 3* np.cos(z.get_value()) -1.5*np.cos(z.get_value()) * np.sin(z.get_value()/2) , 3 *np.sin(z.get_value()) -1.5*np.sin(z.get_value())*np.sin(z.get_value()/2) , + 1.5*np.cos(z.get_value()/2)   ) , color=RED  )))
        p.add_updater( lambda x: x.become( Sphere(radius = 0.1).set_color(ORANGE).shift(axes.c2p( 3* np.cos(z.get_value()) , 3* np.sin(z.get_value()), 0  ))))
        self.play(z.animate.set_value(TAU), run_time = 8)
        self.wait()


class ord(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        d = Tex(r"Definition", font_size = 34, color = BLUE).to_edge(UL).shift(0.3*RIGHT+ 0.5* DOWN)
        l = d.get_left()[0]
        d1 = Tex(r"A surface ", r"$\Sigma$", r" is ", r"orientable", r" if there is a continuous map $N : \Sigma \to \mathbb{S}^2$ such that", font_size = 34).next_to(d, DOWN)
        d1.shift((d1.get_left()[0]-l)*LEFT)
        d2 = Tex(r"$N (p ) \perp T_p \Sigma$ for all $p \in \Sigma$.", font_size = 34).next_to(d1, DOWN)
        d2.shift((d2.get_center()[0])*LEFT)
        d1[3].set_color(BLUE)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE_B).next_to(d1[1], DOWN).shift(0.15*UP)
        ax1 = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 4, y_length = 4, z_length = 2 )
        ax2 = ax1.copy()
        ax1.shift([1.7,-1.7*sqrt(3), -1 ])
        ax2.shift([-1.7,1.7*sqrt(3), -1 ])
        s = Surface(
            lambda u, v: ax1.c2p( u , v , u**2  / 8 + v**2 /8 - 1  ),
            u_range = [-PI, PI],
            v_range = [-PI, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sp = Surface(
            lambda u, v: ax2.c2p( 2.5*np.sin(u)*np.cos(v) ,2.5* np.sin(u)*np.sin(v) ,2.5* np.cos(u)),
            u_range = [0, PI],
            v_range = [-PI, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        self.add_fixed_in_frame_mobjects(d, d1, d2, sli)
        self.remove(d,d1,d2,sli)
        self.play(Create(d), Create(d1), Create(d2), Create(sli), Create(s), Create(sp))
        self.wait()
        n = Arrow(ax1.c2p(0, 0, -1.2) , ax1.c2p(0,0,1.5), color = RED  )
        n2 = Arrow(ax2.c2p(0, 0, -0.2) , ax2.c2p(0,0,2.5), color = RED  )
        p = Sphere(radius = 0.05 ).set_color(ORANGE).shift(ax1.c2p(0,0,-1))
        self.play(Create(n), Create(p), Create(n2))
        z = ValueTracker(0)
        def f(u,v):
            return u**2/8 + v**2/8 - 1
        def nor1(u,v):
            return  -2.5*u/sqrt(16 + u**2 + v**2 ) 
        def nor2(u,v):
            return  -2.5*v/sqrt(16 + u**2 + v**2)
        def nor3(u,v):
            return 10/sqrt(16 + u**2 + v**2)
        def xa(t):
            return t*np.cos(t)/5
        def y(t):
            return t*np.sin(t)/5
        n.add_updater( lambda x: x.become( Arrow( ax1.c2p( xa(z.get_value()) - 0.1*nor1(xa(z.get_value()), y(z.get_value())) , y(z.get_value()) - 0.1*nor2(xa(z.get_value()), y(z.get_value()))   , f(xa(z.get_value()),y(z.get_value())) - 0.1*nor3(xa(z.get_value()), y(z.get_value())) )  ,ax1.c2p( xa(z.get_value()) + 1.1*nor1(xa(z.get_value()), y(z.get_value())) , y(z.get_value()) + 1.1*nor2(xa(z.get_value()), y(z.get_value()))   , f(xa(z.get_value()),y(z.get_value())) + 1.1* nor3(xa(z.get_value()), y(z.get_value())) ) , color=RED  )))
        n2.add_updater( lambda x: x.become( Arrow( ax2.c2p( - 0.1*nor1(xa(z.get_value()), y(z.get_value())) , - 0.1*nor2(xa(z.get_value()), y(z.get_value()))   , - 0.1*nor3(xa(z.get_value()), y(z.get_value())) )  ,ax2.c2p( 1.1*nor1(xa(z.get_value()), y(z.get_value())) ,1.1* nor2(xa(z.get_value()), y(z.get_value())),1.1* nor3(xa(z.get_value()), y(z.get_value())) ) , color=RED  )))
        p.add_updater( lambda x: x.become( Sphere(radius = 0.05).set_color(ORANGE).shift(ax1.c2p( xa(z.get_value()), y(z.get_value()), f(xa(z.get_value()), y(z.get_value())) ))))
        self.play(z.animate.set_value(2*TAU), run_time = 8)
        self.wait()


class nsm(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        t = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift(0.2*RIGHT+ 0.3* DOWN)
        l = t.get_left()[0]
        t1 = Tex(r"$N : \Sigma \to \mathbb{S}^2$ is smooth", font_size = 34).next_to(t, DOWN)
        t1.shift((t1.get_left()[0]-l)*LEFT)
        p = Tex(r"Proof", color = BLUE, font_size = 34).next_to(t1, DOWN)
        p.shift((p.get_left()[0]-l)*LEFT)
        p1 = Tex(r"Take $\phi : U \to \Sigma$ a chart.", r" Then $N = \pm \frac{ \frac{\partial \phi }{\partial u} \times \frac{\partial \phi }{\partial v} }{ \vert \frac{\partial \phi }{\partial u} \times \frac{\partial \phi }{\partial v} \vert }$", font_size = 34).next_to(p, DOWN)
        p1.shift((p1.get_left()[0]-l)*LEFT + 0.2*UP )
        phi = CurvedArrow(3.5*LEFT + 0.7*DOWN, 1.5*LEFT + 0.25*DOWN , radius= -2.5)
        n = CurvedArrow(1.5* RIGHT+ 0.25*DOWN, 3.5*RIGHT + 0.7*DOWN, radius= -2.5)
        phil = Tex(r"$\phi$", font_size = 34).next_to(phi, UP)
        nl = Tex(r"$N$", font_size = 34).next_to(n, UP)
        ax2 = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-1,2,10], x_length = 3, y_length = 3, z_length = 1.5 , tips = False)
        ax1 = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-1,2,10], x_length = 3, y_length = 3, z_length = 0.01, tips = False)
        ax3 = ax2.copy()
        ax1.shift([2.5,-2.5*sqrt(3), -1.5 ])
        ax3.shift([-2.5,2.5*sqrt(3), -1.5 ])
        dom = Surface(
            lambda u, v: ax1.c2p( u , v , 0 ),
            u_range = [ -2 , 2 ],
            v_range = [- 2 , 2 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        s = Surface(
            lambda u, v: ax2.c2p( u , v , - u**2 /3    -  v**2 /3        ),
            u_range = [-2, 2],
            v_range = [-2, 2],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sp = Surface(
            lambda u, v: ax3.c2p( 2.5*np.sin(u)*np.cos(v) ,2.5* np.sin(u)*np.sin(v) ,2.5* np.cos(u)),
            u_range = [0, PI],
            v_range = [-PI, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        x = Sphere(radius = 0.06)
        x.set_color(ORANGE)
        y = x.copy()
        x.shift(ax1.c2p(0,1,0))
        y.shift(ax2.c2p(0,1,-1/3))
        nx = Arrow(ax2.c2p(0,1 - 0.1,-1/3 - 0.15), ax2.c2p(0,1 + (2.8)* 2/sqrt(13) ,-1/3 + (2.8)* 3/sqrt(13)), color = RED)
        ns = Arrow(ax3.c2p(0, - 0.1, - 0.15), ax3.c2p(0,(2.8)* 2/sqrt(13) , (2.8)* 3/sqrt(13)), color = RED)
        e1 = Arrow(ax1.c2p(-0.2 ,1 ,0 ), ax1.c2p( 2 , 1, 0 ), color = YELLOW)
        e2 = Arrow(ax1.c2p(0 ,1-0.2 ,0 ), ax1.c2p( 0 , 3, 0 ), color = YELLOW)
        par1 = Arrow(ax2.c2p(-0.2 ,1 ,-1/3 ), ax2.c2p( 2.5 , 1, -1/3 ), color = YELLOW)
        par2 = Arrow(ax2.c2p(0 , 1 -0.2 ,-1/3 ), ax2.c2p( 0 , 1 + 2.7*3 /sqrt(13), -1/3 - 2.7*2/sqrt(13) ), color = YELLOW)
        self.add_fixed_in_frame_mobjects(t, t1, p, p1, phi, n, phil, nl)
        self.remove(t, t1, p, p1, phi, n, phil, nl)
        self.play(Create(t), Create(t1), Create(s), Create(sp), Create(n), Create(nl))
        self.wait()
        self.play(Create(dom), Create(ax1), Create(phi), Create(phil), Write(p), Write(p1[0]))
        self.wait()
        self.play(Create(x), Create(y))
        self.wait()
        self.play(Create(nx), Create(ns))
        self.wait()
        self.play(Create(e1), Create(e2), Create(par1), Create(par2), Write(p1[1]))
        self.wait()


class gad(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        d = Tex(r"Definition", font_size = 34, color = BLUE).to_edge(UL).shift(0.3*RIGHT+  DOWN)
        l = d.get_left()[0]
        d1 = Tex(r"For an orientable surface $\Sigma$, the map $N : \Sigma \to \mathbb{S}^2$ is called the ", r"Gauss map.", font_size = 34).next_to(d, DOWN)
        d1.shift((d1.get_left()[0]-l)*LEFT)
        d1[1].set_color(BLUE)
        ax1 = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 4, y_length = 4, z_length = 2 )
        ax2 = ax1.copy()
        ax1.shift([1.7,-1.7*sqrt(3), -1 ])
        ax2.shift([-1.7,1.7*sqrt(3), -1 ])
        s = Surface(
            lambda u, v: ax1.c2p( u , v , u**2  / 8 + v**2 /8 - 1  ),
            u_range = [-PI, PI],
            v_range = [-PI, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sp = Surface(
            lambda u, v: ax2.c2p( 2.5*np.sin(u)*np.cos(v) ,2.5* np.sin(u)*np.sin(v) ,2.5* np.cos(u)),
            u_range = [0, PI],
            v_range = [-PI, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        self.add_fixed_in_frame_mobjects(d, d1)
        self.remove(d,d1)
        na = CurvedArrow(0.7*LEFT, 1.3*RIGHT, radius= -3)
        nal = Tex(r"$N$", font_size = 34).next_to(na, UP)
        self.add_fixed_in_frame_mobjects(d, d1, na, nal)
        self.remove(d,d1, na, nal)
        n = Arrow(ax1.c2p(0, 0, -1.2) , ax1.c2p(0,0,1.5), color = RED  )
        n2 = Arrow(ax2.c2p(0, 0, -0.2) , ax2.c2p(0,0,2.5), color = RED  )
        p = Sphere(radius = 0.05 ).set_color(ORANGE).shift(ax1.c2p(0,0,-1))
        self.play(Create(d), Create(d1), Create(s), Create(sp), Create(n), Create(p), Create(n2), Create(na), Create(nal))
        self.wait()




class nconst(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        pl = Tex(r"If $\Sigma$ is a plane, then $N : \Sigma \to \mathbb{S}^2$ is constant.", font_size = 34).to_edge(UL).shift(0.3*RIGHT+  1.5*DOWN)
        ax1 = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 4, y_length = 4, z_length = 2 )
        ax2 = ax1.copy()
        ax1.shift([1.7,-1.7*sqrt(3), -1 ])
        ax2.shift([-1.7,1.7*sqrt(3), -1 ])
        s = Surface(
            lambda u, v: ax1.c2p( u , v , 0 ),
            u_range = [-PI, PI],
            v_range = [-PI, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sp = Surface(
            lambda u, v: ax2.c2p( 2.5*np.sin(u)*np.cos(v) ,2.5* np.sin(u)*np.sin(v) ,2.5* np.cos(u)),
            u_range = [0, PI],
            v_range = [-PI, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        self.remove(pl)
        na = CurvedArrow(0.9*LEFT, 1.1*RIGHT, radius= -3)
        nal = Tex(r"$N$", font_size = 34).next_to(na, UP)
        n = Arrow(ax1.c2p(0.8*PI,0,-0.15), ax1.c2p(0.8*PI,0,2.7), color = RED)
        n2 = Arrow(ax2.c2p(0,0,-0.15), ax2.c2p(0,0,2.7), color = RED)
        z = ValueTracker(PI)
        n.add_updater( lambda x: x.become( Arrow( ax1.c2p( 0.8*z.get_value() , 2* np.sin( - 3*z.get_value()/2 + 3*PI/2 ) , -0.15 )  ,ax1.c2p( 0.8*z.get_value() , 2* np.sin( - 3*z.get_value()/2 + 3*PI/2 ) , 2.7 ) , color=RED  )))
        pt = Sphere(radius = 0.06)
        pt.set_color(ORANGE).shift(ax1.c2p(0.8*PI, 0, 0))
        pt.add_updater( lambda x: x.become( Sphere( radius = 0.06 ).set_color(ORANGE).shift(ax1.c2p( 0.8*z.get_value() , 2* np.sin( - 3*z.get_value()/2 + 3*PI/2 ) , 0 ))  ))
        self.add_fixed_in_frame_mobjects(pl, na, nal)
        self.remove(pl, na, nal)
        self.play(Create(pl), Create(s), Create(sp), Create(na), Create(nal), Create(pt), Create(n), Create(n2))
        self.wait()
        self.play(z.animate.set_value(-PI), run_time = 4)
        self.wait()


class tps(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        dp = Tex(r"The derivative is a linear map $d_pN : T_p \Sigma \to T_{N(p)}\mathbb{S}^2$", font_size = 34).to_edge(UL).shift(0.3*RIGHT+ 1.3* DOWN)
        ax1 = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 4, y_length = 4, z_length = 2 )
        ax2 = ax1.copy()
        ax1.shift([1.7,-1.7*sqrt(3), -1.5 ])
        ax2.shift([-1.7,1.7*sqrt(3), -1 ])
        s = Surface(
            lambda u, v: ax1.c2p( (2 +  np.cos(v))*np.cos(u) , (2 +  np.cos(v))*np.sin(u) , v ),
            u_range = [ 0, PI ],
            v_range = [ 0 , PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sp = Surface(
            lambda u, v: ax2.c2p( 2.5*np.sin(u)*np.cos(v) ,2.5* np.sin(u)*np.sin(v) ,2.5* np.cos(u)),
            u_range = [0, PI],
            v_range = [-PI, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tp1 = Surface(
            lambda u, v: ax1.c2p( sqrt(3)*u/sqrt(2) , 2 + 1/sqrt(2) - v/sqrt(2) , PI/4 + v ),
            u_range = [ -2 , 2  ],
            v_range = [-2 , 2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        tp2 = Surface(
            lambda u, v: ax2.c2p( sqrt(3)*u/sqrt(2) ,  2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) - v/sqrt(2) , 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) + v ),
            u_range = [ -2 , 2  ],
            v_range = [-2 , 2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        na = CurvedArrow(0.7*LEFT, 1.3*RIGHT, radius= -3)
        nal = Tex(r"$N$", font_size = 34).next_to(na, UP)
        self.add_fixed_in_frame_mobjects(dp, na, nal)
        self.remove(dp, na, nal)
        p = Sphere(radius = 0.07 ).set_color(ORANGE).shift(ax1.c2p(0, 2 + 1/sqrt(2), PI/4))
        n = Arrow(ax1.c2p( 0, 2 + 1/sqrt(2) - (4 + sqrt(2))/20 , PI/4 - (1 + 2 *sqrt(2))/20 ) , ax1.c2p(0, 2 + 1/sqrt(2) + 2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) , PI/4 + 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) ), color = RED  )
        n2 = Arrow(ax2.c2p(0, -(4 + sqrt(2))/20 , - (1 + 2 *sqrt(2))/20 ) , ax2.c2p( 0, 2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) , 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) ), color = RED  )
        self.play(Create(dp), Create(s), Create(sp), Create(n), Create(p), Create(n2), Create(na), Create(nal))
        self.wait()
        self.play(Create(tp1), Create(tp2))
        self.wait()

class tqs2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        dp = Tex(r"$ T_{q}\mathbb{S}^2 = \{ v \in \mathbb{R}^3 \vert v \cdot q = 0 \}$", font_size = 34).shift(2.5*UP)
        ax = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 6, y_length = 6, z_length = 3 )
        ax.shift([0,0,-1.1])
        sp = Surface(
            lambda u, v: ax.c2p( 2.5*np.sin(u)*np.cos(v) ,2.5* np.sin(u)*np.sin(v) ,2.5* np.cos(u)),
            u_range = [0, PI],
            v_range = [-PI, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tp2 = Surface(
            lambda u, v: ax.c2p( sqrt(3)*u/sqrt(2) ,  2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) - v/sqrt(2) ,  2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) + v ),
            u_range = [ -2 , 2  ],
            v_range = [-2 , 2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        q = Sphere(radius = 0.06)
        q.set_color(ORANGE).shift(ax.c2p( 0, 2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) , 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) ))
        self.add_fixed_in_frame_mobjects(dp)
        self.remove(dp)
        self.play(Create(sp), Create(tp2), Create(dp), Create(q))
        self.wait()


class tps2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        dp = Tex(r"The derivative is a linear map $d_pN : T_p \Sigma \to T_{N(p)}\mathbb{S}^2$ ", r"$ = T_p \Sigma $", font_size = 34).to_edge(UL).shift(0.3*RIGHT+ 1.3* DOWN)
        ax1 = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 4, y_length = 4, z_length = 2 )
        ax2 = ax1.copy()
        ax1.shift([1.7,-1.7*sqrt(3), -1.5 ])
        ax2.shift([-1.7,1.7*sqrt(3), -1 ])
        s = Surface(
            lambda u, v: ax1.c2p( (2 +  np.cos(v))*np.cos(u) , (2 +  np.cos(v))*np.sin(u) , v ),
            u_range = [ 0, PI ],
            v_range = [ 0 , PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sp = Surface(
            lambda u, v: ax2.c2p( 2.5*np.sin(u)*np.cos(v) ,2.5* np.sin(u)*np.sin(v) ,2.5* np.cos(u)),
            u_range = [0, PI],
            v_range = [-PI, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tp1 = Surface(
            lambda u, v: ax1.c2p( sqrt(3)*u/sqrt(2) , 2 + 1/sqrt(2) - v/sqrt(2) , PI/4 + v ),
            u_range = [ -2 , 2  ],
            v_range = [-2 , 2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        tp2 = Surface(
            lambda u, v: ax2.c2p( sqrt(3)*u/sqrt(2) ,  2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) - v/sqrt(2) , 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) + v ),
            u_range = [ -2 , 2  ],
            v_range = [-2 , 2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        na = CurvedArrow(0.7*LEFT, 1.3*RIGHT, radius= -3)
        nal = Tex(r"$N$", font_size = 34).next_to(na, UP)
        p = Sphere(radius = 0.07 ).set_color(ORANGE).shift(ax1.c2p(0, 2 + 1/sqrt(2), PI/4))
        q = Sphere(radius = 0.07 ).set_color(ORANGE).shift(ax2.c2p( 0, 2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) , 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) ))
        n = Arrow(ax1.c2p( 0, 2 + 1/sqrt(2) - (4 + sqrt(2))/20 , PI/4 - (1 + 2 *sqrt(2))/20 ) , ax1.c2p(0, 2 + 1/sqrt(2) + 2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) , PI/4 + 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) ), color = RED  )
        n2 = Arrow(ax2.c2p(0, -(4 + sqrt(2))/20 , - (1 + 2 *sqrt(2))/20 ) , ax2.c2p( 0, 2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) , 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) ), color = RED  )
        self.add_fixed_in_frame_mobjects(dp, na, nal)
        self.remove(dp[1])
        self.add(s, sp, n, p, n2, na, nal, tp1, tp2, q)
        self.wait()
        self.play(Write(dp[1]))
        self.wait()
        self.play(FadeOut(dp))

class shape(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sh = Tex(r"If $\Sigma$ is an oriented surface, its ", r"shape operator", r" at $p \in \Sigma$ is defined as", font_size = 34).to_edge(UL).shift(0.3*RIGHT+ 1.3* DOWN)
        l = sh.get_left()[0]
        sh1 = Tex(r"$S: T_p \Sigma \to T_p \Sigma$, $S(v ) = - d_pN(v)$", font_size = 34).next_to(sh, DOWN)
        sh1.shift((sh1.get_center()[0])*LEFT)
        ax1 = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 4, y_length = 4, z_length = 2 )
        ax2 = ax1.copy()
        ax1.shift([1.7,-1.7*sqrt(3), -1.5 ])
        ax2.shift([-1.7,1.7*sqrt(3), -1 ])
        s = Surface(
            lambda u, v: ax1.c2p( (2 +  np.cos(v))*np.cos(u) , (2 +  np.cos(v))*np.sin(u) , v ),
            u_range = [ 0, PI ],
            v_range = [ 0 , PI ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sp = Surface(
            lambda u, v: ax2.c2p( 2.5*np.sin(u)*np.cos(v) ,2.5* np.sin(u)*np.sin(v) ,2.5* np.cos(u)),
            u_range = [0, PI],
            v_range = [-PI, PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tp1 = Surface(
            lambda u, v: ax1.c2p( sqrt(3)*u/sqrt(2) , 2 + 1/sqrt(2) - v/sqrt(2) , PI/4 + v ),
            u_range = [ -2 , 2  ],
            v_range = [-2 , 2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        tp2 = Surface(
            lambda u, v: ax2.c2p( sqrt(3)*u/sqrt(2) ,  2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) - v/sqrt(2) , 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) + v ),
            u_range = [ -2 , 2  ],
            v_range = [-2 , 2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        na = CurvedArrow(0.7*LEFT, 1.3*RIGHT, radius= -3)
        nal = Tex(r"$N$", font_size = 34).next_to(na, UP)
        p = Sphere(radius = 0.07 ).set_color(ORANGE).shift(ax1.c2p(0, 2 + 1/sqrt(2), PI/4))
        q = Sphere(radius = 0.07 ).set_color(ORANGE).shift(ax2.c2p( 0, 2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) , 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) ))
        n = Arrow(ax1.c2p( 0, 2 + 1/sqrt(2) - (4 + sqrt(2))/20 , PI/4 - (1 + 2 *sqrt(2))/20 ) , ax1.c2p(0, 2 + 1/sqrt(2) + 2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) , PI/4 + 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) ), color = RED  )
        n2 = Arrow(ax2.c2p(0, -(4 + sqrt(2))/20 , - (1 + 2 *sqrt(2))/20 ) , ax2.c2p( 0, 2.8*(4 + sqrt(2))/sqrt(27+12*sqrt(2)) , 2.8*(1 + 2 *sqrt(2))/sqrt(27+12*sqrt(2)) ), color = RED  )
        self.add_fixed_in_frame_mobjects(sh, sh1, na, nal)
        self.add(s, sp, n, p, n2, na, nal, tp1, tp2, q)
        self.remove(sh,sh1)
        self.wait()
        self.play(Create(sh), Create(sh1))
        self.wait()


class sad(Scene): 
    def construct(self):       
        pro = Tex(r"Proposition", color = BLUE, font_size = 34).to_edge(UL).shift(0.3*RIGHT + 0.5*DOWN)
        l = pro.get_left()[0]
        pro1 = Tex(r"The shape operator $S : T_p \Sigma \to T_p \Sigma$ is self adjoint. That is, ", font_size = 34).next_to(pro, DOWN)
        pro1.shift((pro1.get_left()[0]-l)*LEFT)
        pro2 = Tex(r"$S (V) \cdot W = V \cdot S (W)$ for all $V,W \in T_p\Sigma$", font_size = 34).next_to(pro1, DOWN)
        pro2.shift((pro2.get_center()[0])*LEFT)
        pf = Tex(r"Proof", color = BLUE, font_size = 34).next_to(pro2, DOWN)
        pf.shift((pf.get_left()[0]-l)*LEFT)
        pf1 = Tex(r"For a parametrization $\phi : U \to \Sigma$, set ",  font_size = 34).next_to(pf, DOWN)
        pf1.shift((pf1.get_left()[0]-l)*LEFT)
        pf2 = Tex(r"$\phi_u =  \frac{\partial \phi }{\partial u}$, $\phi_v =  \frac{\partial \phi }{\partial v}$, $N_u = \frac{\partial N}{\partial u}$, $N_v = \frac{\partial N}{\partial v}$.", font_size = 34).next_to(pf1, DOWN)
        pf2.shift((pf2.get_center()[0])*LEFT)
        pf3 = Tex(r"$S(\phi_u) \cdot \phi_v $ ", r"$= - N_u \cdot \phi_v$ ",r"$ = - (N \cdot \phi_v)_u + N \cdot \phi_{vu}$ ",   font_size = 34).next_to(pf2, DOWN)
        pf3.shift((pf3.get_left()[0]-l-1.5)*LEFT+0.1*DOWN)
        pf4 = Tex(r"$= N \cdot \phi_{uv}$ " ,r"$ = (N \cdot \phi_{u})_v - N_v \cdot \phi_u$ ", r"$= S(\phi _v) \cdot  \phi_u$." ,  font_size = 34).next_to(pf3, DOWN)
        pf4.shift((pf4.get_left()[0]-l -3.5)*LEFT + 0.1 *DOWN)
        pf5 = Tex(r"$S(\phi_u) \cdot \phi_v $ ", r"$= S(\phi _v) \cdot  \phi_u$." ,  font_size = 34).next_to(pf4, DOWN)
        pf5.shift((pf5.get_center()[0])*LEFT + 0.1 *DOWN)
        pf6 = pf5.copy()
        pf6.shift(2.8*UP)
        pf7 = Tex(r"For $V, W \in T_p\Sigma$, we can write $V= V_1 \phi_u + V_2 \phi _v$, $W = W_1 \phi_u + W_2 \phi_v$," ,  font_size = 34).next_to(pf6, DOWN)
        pf7.shift((pf7.get_left()[0]-l)*LEFT )
        pf8 = Tex(r"$S(V) \cdot W$ ", r"$ =S(V_1 \phi_u + V_2 \phi_v) \cdot (W_1 \phi_u + W_2 \phi_v) $ ", font_size = 34).next_to(pf7, DOWN)
        pf8.shift((pf8.get_left()[0]-l)*LEFT )        
        pf9 = Tex( r"$= V_1 W_1 S(\phi_u)\cdot \phi_u + V_1W_2 S(\phi _u) \cdot \phi_v  +  V_2 W_1 S(\phi_v)\cdot \phi_u  + V_2W_2 S(\phi_v) \cdot \phi_v $" ,  font_size = 34).next_to(pf8, DOWN)
        pf9.shift((pf9.get_left()[0]-l)*LEFT )        
        pf10 = Tex( r"$= V_1 W_1 S(\phi_u)\cdot \phi_u + V_1W_2 S(\phi _v) \cdot \phi_u  +  V_2 W_1 S(\phi_u)\cdot \phi_v  + V_2W_2 S(\phi_v) \cdot \phi_v $" ,  font_size = 34).next_to(pf9, DOWN)
        pf10.shift((pf10.get_left()[0]-l)*LEFT )        
        pf11 = Tex( r"$= (V_1 \phi_u + V_2 \phi_v ) \cdot S ( W_1 \phi_u + W_2 \phi_v ) $ ", r"$ = V \cdot S(W)$." ,  font_size = 34).next_to(pf10, DOWN)
        pf11.shift((pf11.get_left()[0]-l)*LEFT )        
        cros = Line([-0.7,-1,0], [0.7,-0.7,0], color = RED)
        cros2 = Line([-0.9,-1.7,0], [0.5,-1.4,0], color = RED)
        suv = Line([0,0,0], [3.5,0,0], color = GREEN).next_to(pf6, DOWN).shift(0.15*UP)
        suv1 = Line([0,0,0], [1.4,0,0], color = GREEN).next_to(pf9, DOWN).shift(0.15*UP + 0.7*LEFT)
        suv2 = Line([0,0,0], [1.4,0,0], color = GREEN).next_to(pf9, DOWN).shift(0.15*UP + 2*RIGHT)
        suv3 = Line([0,0,0], [1.4,0,0], color = GREEN).next_to(pf10, DOWN).shift(0.15*UP + 0.7*LEFT)
        suv4 = Line([0,0,0], [1.4,0,0], color = GREEN).next_to(pf10, DOWN).shift(0.15*UP + 2*RIGHT)
        self.play(Write(pro), Write(pro1), Write(pro2))
        self.wait()
        self.play(Write(pf), Write(pf1), Write(pf2))
        self.wait()
        self.play(Write(pf3[0]))
        self.wait()
        self.play(Write(pf3[1]))
        self.wait()
        self.play(Write(pf3[2]))
        self.wait()
        self.play(Create(cros))
        self.wait()
        self.play(Write(pf4[0]))
        self.wait()
        self.play(Write(pf4[1]))
        self.wait()
        self.play(Create(cros2))
        self.wait()
        self.play(Write(pf4[2]))
        self.wait()
        self.play(Transform(pf3[0], pf5[0]), Transform(pf4[2], pf5[1]))
        self.wait()
        self.add(pf5)
        self.remove(pf3[0], pf4[2])
        self.play(FadeOut(pf1), FadeOut(pf2), FadeOut(pf3[1]), FadeOut(pf3[2]), FadeOut(pf4[0]), FadeOut(pf4[1]), FadeOut(cros), FadeOut(cros2), Transform(pf5, pf6), Create(suv))
        self.wait()
        self.play(Write(pf7))
        self.wait()
        self.play(Write(pf8[0]))
        self.wait()
        self.play(Write(pf8[1]))
        self.wait()
        self.play(Write(pf9))
        self.wait()
        self.play(Create(suv1), Create(suv2))
        self.wait()
        self.play(Write(pf10), Create(suv3), Create(suv4))
        self.wait()
        self.play(Write(pf11[0]))
        self.wait()
        self.play(Write(pf11[1]))
        self.wait()
        self.play(FadeOut(pf11), FadeOut(pf10),FadeOut(pf9),FadeOut(pf8),FadeOut(pf7),FadeOut(suv),FadeOut(suv1),FadeOut(suv2),FadeOut(suv3),FadeOut(suv4),FadeOut(pf5), FadeOut(pf))
        self.wait()
        

class spec(Scene): 
    def construct(self):       
        pro = Tex(r"Proposition", color = BLUE, font_size = 34).to_edge(UL).shift(0.5*DOWN + 0.3*RIGHT)
        l = pro.get_left()[0]
        pro1 = Tex(r"The shape operator $S : T_p \Sigma \to T_p \Sigma$ is self adjoint. That is, ", font_size = 34).next_to(pro, DOWN)
        pro1.shift((pro1.get_left()[0]-l)*LEFT)
        pro2 = Tex(r"$S (V) \cdot W = V \cdot S (W)$ for all $V,W \in T_p\Sigma$", font_size = 34).next_to(pro1, DOWN)
        pro2.shift((pro2.get_center()[0])*LEFT)
        t = Tex(r"Spectral Theorem", font_size = 34, color = BLUE).next_to(pro2, DOWN)
        t.shift((t.get_left()[0]-l)*LEFT)
        t1 = Tex(r"If $X$ is a finite dimensional vector space with a dot product and $S:X \to X$", font_size = 34).next_to(t, DOWN)
        t1.shift((t1.get_left()[0]-l)*LEFT)
        t2 = Tex(r"a self adjoint linear map, then $X$ admits an orthonormal basis of eigenvectors.", font_size = 34).next_to(t1, DOWN)
        t2.shift((t2.get_left()[0]-l)*LEFT)
        c = Tex(r"Corollary", font_size = 34, color = BLUE).next_to(t2, DOWN)
        c.shift((c.get_left()[0]-l)*LEFT)
        c1 = Tex(r"There are two orthonormal vectors $\{ E_1, E_2 \} \hookrightarrow T_p \Sigma$ for which", font_size = 34).next_to(c, DOWN)
        c1.shift((c1.get_left()[0]-l)*LEFT)
        c2 = Tex(r"$S(E_1) = k_1(E_1)$, $S(E_2) = k_2E_2$, for some $k_1, k_2 \in \mathbb{R}$", r" $(k_1 \geq k_2)$.", font_size = 34).next_to(c1, DOWN)
        c2.shift((c2.get_left()[0]-l)*LEFT)        
        self.add(pro, pro1, pro2)
        self.wait()
        self.play(Write(t), Write(t1), Write(t2))
        self.wait()
        self.play(Write(c), Write(c1), Write(c2[0]))
        self.wait()
        self.play(Write(c2[1]))
        self.wait()
        

class sff(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        kn = Tex(r"Take ", r"$\gamma$", r" a curve in ", r"$\Sigma$", r" parametrized by arc length with $\gamma ( 0 ) =$ ", r"$p$", r" and $\gamma ^{\prime}(0) = $ ", r"$V$.", font_size = 34).to_edge(UL).shift(0.2*RIGHT + 0.8*DOWN)
        l = kn.get_left()[0]
        kn1 = Tex(r"$S(V) \cdot V $ ", r"$ = -dpN(V) \cdot V  $ ", r"$= - N^{\prime} ( 0 ) \cdot \gamma ^{\prime}(0) $ ", r"$= - (N \cdot \gamma ^{\prime})^{\prime }  + $", r" $ N \cdot \gamma^{\prime \prime } (0) $."  , font_size = 34).next_to(kn, DOWN)
        kn1.shift((kn1.get_left()[0]-l)*LEFT+0.2*DOWN)
        kn2 = Tex(r"$S(V) \cdot V $ ", r"$ = $ ", r"$ N \cdot \gamma^{\prime \prime } (0) $."  , font_size = 34).next_to(kn, DOWN)
        kn2.shift((kn2.get_center()[0])*LEFT+0.2*DOWN)
        ax = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 6, y_length = 6, z_length = 3 )
        ax.shift([0,0,-0.35])
        gli = Line([0,0,0], [0.25,0,0], color=BLUE).next_to(kn[1], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.25,0,0], color=ORANGE).next_to(kn[5], DOWN).shift(0.15*UP)
        vli = Line([0,0,0], [0.25,0,0], color=YELLOW).next_to(kn[7], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.25,0,0], color=PURPLE_B).next_to(kn[3], DOWN).shift(0.15*UP)
        cros = Line([0.75,1.55,0], [2.15,1.85,0], color=RED)
        sig = Surface(
            lambda u, v: ax.c2p( u, v , - u**2/3 - v**2/3 ),
            u_range = [ -2 , 2  ],
            v_range = [-2 , 2],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift(ax.c2p( 0,0,0 ))
        n = Arrow(ax.c2p( 0,0,-0.15 ), ax.c2p(0,0,2), color = RED)
        v = Arrow(  ax.c2p( -0.1,-sqrt(3)/10,0 ) ,  ax.c2p( 1,sqrt(3),0 ), color = YELLOW )
        g = ParametricFunction(lambda t : ax.c2p( t, sqrt(3)*t , -4*t**2/3 ), t_range = [-1,1], color = BLUE)
        self.add_fixed_in_frame_mobjects(kn, kn1, kn2, vli, gli, pli, sli, cros)
        self.remove(kn, kn1, vli, gli, pli, sli, cros, kn2)
        self.play(Write(kn), Create(sig), Create(p), Create(v), Create(vli), Create(gli), Create(pli), Create(g), Create(sli), Create(n))
        self.wait()        
        self.play(Write(kn1[0]))
        self.wait()        
        self.play(Write(kn1[1]))
        self.wait()        
        self.play(Write(kn1[2]))
        self.wait()        
        self.play(Write(kn1[3]), Write(kn1[4]))
        self.wait()
        self.play(Create(cros))
        self.wait()
        self.play(Transform( kn1[0], kn2[0] ), Transform(kn1[4], kn2[2]), Write(kn2[1]), FadeOut(kn1[1]), FadeOut(kn1[2]), FadeOut(kn1[3]), FadeOut(cros))
        self.wait()


class pnc(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        ep = Tex(r"Example", font_size = 34, color = BLUE).to_edge(UL).shift(1.3*RIGHT+ 1.5* DOWN)
        l = ep.get_left()[0]
        ep1 = Tex(r"$S(V) \cdot V > 0 $ for all $V \in T_p\Sigma $ with $\vert V \vert =1$.", font_size = 34).next_to(ep, DOWN)
        ep1.shift((ep1.get_center()[0])*LEFT + 0.2 *DOWN)
        ax = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 6, y_length = 6, z_length = 3 )
        ax.shift([0,0,-2.5])
        s = Surface(
            lambda u, v: ax.c2p(  u , v, u**2/3 + v**2/3  ),
            u_range = [ -2, 2 ],
            v_range = [ -2, 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        n = Arrow(ax.c2p(0,0,-0.15 ),ax.c2p(0,0,2.5), color = RED )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift(ax.c2p( 0,0,0 ))
        g = ParametricFunction(lambda t : ax.c2p( t, sqrt(3)*t , 4*t**2/3 ), t_range = [-1,1], color = BLUE)
        g2 = ParametricFunction(lambda t : ax.c2p( 2*t, 0 , 4*t**2/3 ), t_range = [-1,1], color = BLUE)
        g3 = ParametricFunction(lambda t : ax.c2p( sqrt(2)*t, sqrt(2)*t , 4*t**2/3 ), t_range = [-1,1], color = BLUE)
        self.add_fixed_in_frame_mobjects(ep, ep1)
        self.remove(ep, ep1)
        self.play(Create(ep), Create(ep1), Create(s), Create(p), Create(n) )
        self.wait()
        self.play(Create(g))
        self.wait()
        self.play(FadeOut(g))
        self.wait()
        self.play(Create(g2))
        self.wait()
        self.play(FadeOut(g2))
        self.wait()
        self.play(Create(g3))
        self.wait()
        self.play(FadeOut(g3))
        self.wait()
        

class nnc(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        ep = Tex(r"Example", font_size = 34, color = BLUE).to_edge(UL).shift(1.3*RIGHT+  DOWN)
        l = ep.get_left()[0]
        ep1 = Tex(r"$S(V) \cdot V < 0 $ for all $V \in T_p\Sigma $ with $\vert V \vert =1$.", font_size = 34).next_to(ep, DOWN)
        ep1.shift((ep1.get_center()[0])*LEFT + 0.2 *DOWN)
        ax = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 6, y_length = 6, z_length = 3 )
        ax.shift([0,0,-0.45])
        s = Surface(
            lambda u, v: ax.c2p(  u , v, - u**2/3 -  v**2/3  ),
            u_range = [ -2, 2 ],
            v_range = [ -2, 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        n = Arrow(ax.c2p(0,0,-0.15 ),ax.c2p(0,0,2.2), color = RED )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift(ax.c2p( 0,0,0 ))
        g = ParametricFunction(lambda t : ax.c2p( t, sqrt(3)*t ,- 4*t**2/3 ), t_range = [-1,1], color = BLUE)
        g2 = ParametricFunction(lambda t : ax.c2p( 2*t, 0 , - 4*t**2/3 ), t_range = [-1,1], color = BLUE)
        g3 = ParametricFunction(lambda t : ax.c2p( sqrt(2)*t, sqrt(2)*t , - 4*t**2/3 ), t_range = [-1,1], color = BLUE)
        self.add_fixed_in_frame_mobjects(ep, ep1)
        self.remove(ep, ep1)
        self.play(Create(ep), Create(ep1), Create(s), Create(p), Create(n) )
        self.wait()
        self.play(Create(g))
        self.wait()
        self.play(FadeOut(g))
        self.wait()
        self.play(Create(g2))
        self.wait()
        self.play(FadeOut(g2))
        self.wait()
        self.play(Create(g3))
        self.wait()
        self.play(FadeOut(g3))
        self.wait()
        

class mnc(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        ep = Tex(r"Example", font_size = 34, color = BLUE).to_edge(UL).shift(1.3*RIGHT+  DOWN)
        l = ep.get_left()[0]
        ep1 = Tex(r"$S(V) \cdot V > 0 $ for some $V $ and  $S (V) \cdot V < 0 $ for some $V$.", font_size = 34).next_to(ep, DOWN)
        ep1.shift((ep1.get_center()[0])*LEFT + 0.2 *DOWN)
        ax = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 6, y_length = 6, z_length = 3 )
        ax.shift([0,0,-0.45])
        s = Surface(
            lambda u, v: ax.c2p(  u , v, - u**2/3 +  v**2/3  ),
            u_range = [ -2, 2 ],
            v_range = [ -2, 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        n = Arrow(ax.c2p(0,0,-0.15 ),ax.c2p(0,0,2.2), color = RED )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift(ax.c2p( 0,0,0 ))
        g = ParametricFunction(lambda t : ax.c2p( t, sqrt(3)*t ,  2*t**2/3 ), t_range = [-1,1], color = BLUE)
        g2 = ParametricFunction(lambda t : ax.c2p( 2*t, 0 , - 4*t**2/3 ), t_range = [-1,1], color = BLUE)
        self.add_fixed_in_frame_mobjects(ep, ep1)
        self.remove(ep, ep1)
        self.play(Create(ep), Create(ep1), Create(s), Create(p), Create(n) )
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/5)
        self.wait()
        self.play(Create(g))
        self.wait()
        self.play(FadeOut(g))
        self.wait()
        self.play(Create(g2))
        self.wait()
        self.play(FadeOut(g2))
        self.wait()
        self.stop_ambient_camera_rotation()
        self.wait()
        
        
       

class k12(Scene): 
    def construct(self):       
        k = Tex(r"If $\vert V \vert = 1$, then $V =  \cos ( \theta )E_1 +  \sin (\theta ) E_2$ for some $\theta$.", font_size = 34).to_edge(UL).shift(0.4*RIGHT + 0.7*DOWN)
        l = k.get_left()[0]
        k1 = Tex(r"$S(V) \cdot V $", r" $= S ( \cos ( \theta )E_1 +  \sin (\theta )E_2 ) \cdot (  \cos ( \theta )E_1 +  \sin (\theta )E_2 )$", font_size = 34).next_to(k, DOWN)
        k1.shift((k1.get_left()[0]-l-1.5)*LEFT + 3*DOWN)
        l2 = k1[1].get_left()[0]
        k2 = Tex(r"$= ( k_1  \cos ( \theta )E_1 +  k_2 \sin (\theta ) E_2) \cdot ( \cos ( \theta ) E_1 +  \sin (\theta ) E_2)$", font_size = 34).next_to(k1, DOWN)
        k2.shift((k2.get_left()[0]-l2)*LEFT)
        k3 = Tex(r"$= k_1 \cos ^2 (\theta ) + k_2 \sin ^2 (\theta )$", font_size = 34).next_to(k2, DOWN)
        k3.shift((k3.get_left()[0]-l2)*LEFT)
        k11 = Tex(r"$S(V) \cdot V $", r" $=   k_1 \cos^2( \theta) + k_2 \sin^2 (\theta)$", font_size = 34)
        k11.shift((k1.get_center()[1] - k11.get_center()[1])*UP)
        k21 = Tex(r"$k_1$", r" $ \geq $ ", font_size = 34).next_to(k11, LEFT)
        k31 = Tex(r" $ \geq $ ", r"$ k_2 $", font_size = 34).next_to(k11, RIGHT)
        k4 = Tex(r"$S(E_1) \cdot E_1 $", r" $ \geq \text{ } S(V) \cdot V \text{ } \geq $ ", r"$ S(E_2) \cdot E_2$", font_size = 34).next_to(k11, DOWN)
        v = Arrow([-1.2,-0.2,0], [0.5, 1.5, 0], color = YELLOW)
        e1 = Arrow([-1 - 0.2* sqrt(2),0,0], [-1 + 1.5*sqrt(2), 0, 0])
        e2 = Arrow([-1,-0.2* sqrt(2) ,0], [-1, 1.5*sqrt(2), 0])
        pd1 = Line([-1,0,0], [-1 + 1.5*sqrt(2), 0, 0])
        dir = Line([-1,0,0], [0.5, 1.5, 0])
        th = Angle( pd1, dir , radius = 0.6)
        e1l = Tex(r"$E_1$", font_size = 30).next_to(e1, RIGHT)
        e2l = Tex(r"$E_2$", font_size = 30).next_to(e2, LEFT).shift(0.15*RIGHT)
        vl = Tex(r"$V$", font_size = 30).next_to(v, RIGHT).shift(0.32*UP+0.05*LEFT)
        thl = Tex(r"$\theta$", font_size = 30).next_to(th, RIGHT).shift(0.1*UP)
        vli = Line([0,0,0], [0.2,0,0], color = YELLOW).next_to(vl, DOWN).shift(0.15*UP)
        k1li1 = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(k21[0], DOWN).shift(0.15*UP)
        k2li1 = Line([0,0,0], [0.3,0,0], color = PINK).next_to(k31[1], DOWN).shift(0.15*UP)
        k1li2 = Line([0,0,0], [0.8,0,0], color = BLUE).next_to(k4[0], DOWN).shift(0.15*UP)
        k2li2 = Line([0,0,0], [0.8,0,0], color = PINK).next_to(k4[2], DOWN).shift(0.15*UP)
        self.play(Write(k))
        self.wait()
        self.play(Create(v), Create(e1), Create(e2), Create(vl), Create(e1l), Create(e2l), Create(vli), Create(th), Create(thl))
        self.wait()
        self.play(Write(k1), Write(k2), Write(k3))
        self.wait()
        self.play(Transform(k1[0], k11[0]), Transform(k3, k11[1]), FadeOut(k1[1]), FadeOut(k2))
        self.wait()
        self.play(Write(k21), Write(k31), Create(k1li1), Create(k2li1))
        self.wait()
        self.play(Write(k4), Create(k1li2), Create(k2li2))
        self.wait()



class pdir(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70 * DEGREES, theta= 30*DEGREES)
        ax = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 6, y_length = 6, z_length = 3 )
        ax.shift([0,0,-2.2])
        pd1 = Tex(r"$E_1$", r" is the direction in which ", r"$\Sigma$", r" bends the most towards ", r"$N$", font_size = 34).to_edge(UL).shift(0.7*DOWN)
        pd1.shift(pd1.get_center()[0]*LEFT)
        l = pd1.get_left()[0]
        pd2 = Tex(r"$E_2$", r" is the direction in which ", r"$\Sigma$", r" bends the least towards ", r"$N$", font_size = 34).next_to(pd1, DOWN)
        pd2.shift((l - pd2.get_left()[0])*RIGHT)
        e1li = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(pd1[0], DOWN).shift(0.15*UP)
        sli1 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(pd1[2], DOWN).shift(0.15*UP)
        e2li = Line([0,0,0], [0.3,0,0], color = PINK).next_to(pd2[0], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(pd2[2], DOWN).shift(0.15*UP)
        nli1 = Line([0,0,0], [0.3,0,0], color = RED).next_to(pd1[4], DOWN).shift(0.15*UP)
        nli2 = Line([0,0,0], [0.3,0,0], color = RED).next_to(pd2[4], DOWN).shift(0.15*UP)
        s = Surface(
            lambda u, v: ax.c2p(  u , v,  v**2/3  ),
            u_range = [ -3, 3 ],
            v_range = [ -3, 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        n = Arrow(ax.c2p(0,0,-0.2 ),ax.c2p(0,0,2.2), color = RED )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift(ax.c2p( 0,0,0 ))
        g = ParametricFunction(lambda t : ax.c2p( 0, t ,  t**2/3 ), t_range = [-2.5,2.5], color = BLUE)
        g2 = ParametricFunction(lambda t : ax.c2p( t, 0 , 0 ), t_range = [-2.5,2.5], color = PINK)
        self.add_fixed_in_frame_mobjects(pd1, pd2, e1li, e2li, sli1, sli2, nli1, nli2)
        self.remove(pd1, pd2, e1li, e2li, sli1, sli2, nli1, nli2)
        self.play(Write(pd1), Create(s), Create(p), Create(n), Create(e1li), Create(sli1), Create(nli1))
        self.wait()
        self.play(Create(g))
        self.wait()
        self.play(FadeOut(g), Write(pd2), Create(sli2), Create(e2li), Create(nli2))
        self.wait()
        self.play(Create(g2))
        self.wait()
        self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(g2))
        self.wait()


class pdir2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=70 * DEGREES, theta= 30*DEGREES)
        ax = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 6, y_length = 6, z_length = 3 )
        ax.shift([0,0,-1.5])
        pd1 = Tex(r"$E_1$", r" is the direction in which ", r"$\Sigma$", r" bends the most towards ", r"$N$", font_size = 34).to_edge(UL).shift(0.7*DOWN)
        pd1.shift(pd1.get_center()[0]*LEFT)
        l = pd1.get_left()[0]
        pd2 = Tex(r"$E_2$", r" is the direction in which ", r"$\Sigma$", r" bends the least towards ", r"$N$", font_size = 34).next_to(pd1, DOWN)
        pd2.shift((l - pd2.get_left()[0])*RIGHT)
        e1li = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(pd1[0], DOWN).shift(0.15*UP)
        sli1 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(pd1[2], DOWN).shift(0.15*UP)
        e2li = Line([0,0,0], [0.3,0,0], color = PINK).next_to(pd2[0], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(pd2[2], DOWN).shift(0.15*UP)
        nli1 = Line([0,0,0], [0.3,0,0], color = RED).next_to(pd1[4], DOWN).shift(0.15*UP)
        nli2 = Line([0,0,0], [0.3,0,0], color = RED).next_to(pd2[4], DOWN).shift(0.15*UP)
        s = Surface(
            lambda u, v: ax.c2p(  u , v, u**2/3 - v**2/3  ),
            u_range = [ -3, 3 ],
            v_range = [ -3, 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        n = Arrow(ax.c2p(0,0,-0.2 ),ax.c2p(0,0,2.2), color = RED )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift(ax.c2p( 0,0,0 ))
        g = ParametricFunction(lambda t : ax.c2p( t , 0 , t**2/3 ), t_range = [-2.5,2.5], color = BLUE)
        g2 = ParametricFunction(lambda t : ax.c2p( 0 , t , -t**2/3 ), t_range = [-2.5,2.5], color = PINK)
        self.add_fixed_in_frame_mobjects(pd1, pd2, e1li, e2li, sli1, sli2, nli1, nli2)
        self.play(Create(s), Create(p), Create(n))
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/8)
        self.wait()
        self.play(Create(g))
        self.wait()
        self.play(FadeOut(g))
        self.wait()
        self.play(Create(g2))
        self.wait()
        self.stop_ambient_camera_rotation()
        self.wait()
        self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(g2))
        self.wait()
        



class gak(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        ax = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 6, y_length = 6, z_length = 3 )
        ax.shift([0,0,-1.5])
        de = Tex(r"Definition", color = BLUE, font_size = 34).to_edge(UL).shift(0.3*DOWN + 0.3*RIGHT)
        l = de.get_left()[0]
        de1 = Tex(r"The ", r"Gauss curvature", r" at a point ", r"$p$", r" $\in $ ", r"$ \Sigma $", r"  is defined as", font_size = 34).next_to(de, DOWN)
        de1.shift((de1.get_left()[0]-l)*LEFT)
        de1[1].set_color(BLUE)
        de2 = Tex(r"$K ( p ) : = $ det$(S)(p)$", r" $=$ det$(d_pN)$", r" $ = k_1 k_2 $", font_size = 34).next_to(de1, DOWN)
        de2.shift((de2.get_center()[0])*LEFT)    
        ex = Tex(r"$k_1, k_2 > 0 $, so  $K > 0$", font_size = 34).shift(2.6*DOWN)
        pli = Line([0,0,0], [0.3, 0, 0], color = ORANGE).next_to(de1[3], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.4, 0, 0], color = PURPLE_B).next_to(de1[5], DOWN).shift(0.15*UP)
        s = Surface(
            lambda u, v: ax.c2p(  u , v, u**2/3 + v**2/3  ),
            u_range = [ -2, 2 ],
            v_range = [ -2, 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        n = Arrow(ax.c2p(0,0,-0.15 ),ax.c2p(0,0,2.5), color = RED )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift(ax.c2p( 0,0,0 ))
        self.add_fixed_in_frame_mobjects(de, de1, de2, pli, sli, ex)
        self.remove(de, de1, de2, pli, sli, ex)
        self.play(Write(de), Write(de1), Create(pli), Create(sli), Write(de2[0]))
        self.wait()
        self.play( Write(de2[1]))
        self.wait()
        self.play( Write(de2[2]))
        self.wait()
        self.play(Create(s), Create(p), Create(n), Write(ex))
        self.wait()
        self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(ex))
        self.wait()
        

class gak2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        ax = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 6, y_length = 6, z_length = 3 )
        ax.shift([0,0,0.3])
        de = Tex(r"Definition", color = BLUE, font_size = 34).to_edge(UL).shift(0.3*DOWN + 0.3*RIGHT)
        l = de.get_left()[0]
        de1 = Tex(r"The ", r"Gauss curvature", r" at a point ", r"$p$", r" $\in $ ", r"$ \Sigma $", r"  is defined as", font_size = 34).next_to(de, DOWN)
        de1.shift((de1.get_left()[0]-l)*LEFT)
        de1[1].set_color(BLUE)
        de2 = Tex(r"$K ( p ) : = $ det$(S)(p)$", r" $=$ det$(d_pN)$", r" $ = k_1 k_2 $", font_size = 34).next_to(de1, DOWN)
        de2.shift((de2.get_center()[0])*LEFT)    
        ex = Tex(r"$k_1 , k_2  < 0 $, so  $K > 0$", font_size = 34).shift(2.75*DOWN)
        pli = Line([0,0,0], [0.3, 0, 0], color = ORANGE).next_to(de1[3], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.4, 0, 0], color = PURPLE_B).next_to(de1[5], DOWN).shift(0.15*UP)
        s = Surface(
            lambda u, v: ax.c2p(  u , v,  - u**2/3 - v**2/3  ),
            u_range = [ -2, 2 ],
            v_range = [ -2, 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        n = Arrow(ax.c2p(0,0,-0.15 ),ax.c2p(0,0,1.8), color = RED )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift(ax.c2p( 0,0,0 ))
        self.add_fixed_in_frame_mobjects(de, de1, de2, pli, sli, ex)
        self.remove(ex)
        self.play(Write(de), Write(de1), Create(pli), Create(sli), Write(de2[0]))
        self.wait()
        self.play(Create(s), Create(p), Create(n), Write(ex))
        self.wait()
        self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(ex))
        self.wait()
    

class gak3(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        ax = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 6, y_length = 6, z_length = 3 )
        ax.shift([0,0,-0.3])
        de = Tex(r"Definition", color = BLUE, font_size = 34).to_edge(UL).shift(0.3*DOWN + 0.3*RIGHT)
        l = de.get_left()[0]
        de1 = Tex(r"The ", r"Gauss curvature", r" at a point ", r"$p$", r" $\in $ ", r"$ \Sigma $", r"  is defined as", font_size = 34).next_to(de, DOWN)
        de1.shift((de1.get_left()[0]-l)*LEFT)
        de1[1].set_color(BLUE)
        de2 = Tex(r"$K ( p ) : = $ det$(S)(p)$", r" $=$ det$(d_pN)$", r" $ = k_1 k_2 $", font_size = 34).next_to(de1, DOWN)
        de2.shift((de2.get_center()[0])*LEFT)    
        ex = Tex(r"$k_1 > 0 > k_2 $, so  $K < 0$", font_size = 34).shift(2.75*DOWN)
        pli = Line([0,0,0], [0.3, 0, 0], color = ORANGE).next_to(de1[3], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.4, 0, 0], color = PURPLE_B).next_to(de1[5], DOWN).shift(0.15*UP)
        s = Surface(
            lambda u, v: ax.c2p(  u , v,  - u**2/3 + v**2/3  ),
            u_range = [ -2, 2 ],
            v_range = [ -2, 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        n = Arrow(ax.c2p(0,0,-0.15 ),ax.c2p(0,0,1.8), color = RED )
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift(ax.c2p( 0,0,0 ))
        self.add_fixed_in_frame_mobjects(de, de1, de2, pli, sli, ex)
        self.remove(ex)
        self.play(Write(de), Write(de1), Create(pli), Create(sli), Write(de2[0]))
        self.wait()
        self.play(Create(s), Create(p), Create(n), Write(ex))
        self.wait()
        self.play(FadeOut(s), FadeOut(p), FadeOut(n), FadeOut(ex))
        self.wait()
        


class gak4(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=62 * DEGREES, theta= 30*DEGREES)
        ax = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-2,2,10], x_length = 5, y_length = 5, z_length = 5/2 )
        ax.shift([0,0,-0.3])
        de = Tex(r"Definition", color = BLUE, font_size = 34).to_edge(UL).shift(0.3*DOWN + 0.3*RIGHT)
        l = de.get_left()[0]
        de1 = Tex(r"The ", r"Gauss curvature", r" at a point ", r"$p$", r" $\in $ ", r"$ \Sigma $", r"  is defined as", font_size = 34).next_to(de, DOWN)
        de1.shift((de1.get_left()[0]-l)*LEFT)
        de1[1].set_color(BLUE)
        de2 = Tex(r"$K ( p ) : = $ det$(S)(p)$", r" $=$ det$(d_pN)$", r" $ = k_1 k_2 $", font_size = 34).next_to(de1, DOWN)
        de2.shift((de2.get_center()[0])*LEFT)    
        ex = Tex(r"$K > 0 $", r" on the outside, ", r"$K < 0$", r" on the inside", font_size = 34).shift(2.75*DOWN)
        pli = Line([0,0,0], [0.3, 0, 0], color = ORANGE).next_to(de1[3], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.4, 0, 0], color = PURPLE_B).next_to(de1[5], DOWN).shift(0.15*UP)
        spli = Line([0,0,0], [0.8, 0, 0], color = BLUE).next_to(ex[0], DOWN).shift(0.15*UP)
        snli = Line([0,0,0], [0.8, 0, 0], color = GREEN).next_to(ex[2], DOWN).shift(0.15*UP)
        sp = Surface(
            lambda u, v: ax.c2p(  (3 + np.cos(v))*np.cos(u) , (3 + np.cos(v))*np.sin(u), np.sin(v) ),
            u_range = [ 0, TAU ],
            v_range = [ - PI/2 , PI/2 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        sn = Surface(
            lambda u, v: ax.c2p(  (3 + np.cos(v))*np.cos(u) , (3 + np.cos(v))*np.sin(u), np.sin(v) ),
            u_range = [ 0, TAU ],
            v_range = [  PI/2 , 3*PI/2 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        self.add_fixed_in_frame_mobjects(de, de1, de2, pli, sli, ex, spli, snli)
        self.remove(ex, spli, snli)
        self.wait()
        self.play(Create(sp), Create(sn), Write(ex), Create(spli), Create(snli))
        self.wait()
        