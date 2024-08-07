from this import d
from tkinter import E
from manim import *
from numpy import sqrt
import math



class ti(Scene):
    def construct(self):
        t1 = Text("Surfaces", font_size=60)
        self.play(Write(t1))
        self.wait()        


class sdef(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-6,6,1], y_range = [-6,6,1], z_range = [-2,3,1], x_length = 9, y_length = 9, z_length = 15/4 )
        axes.shift([0,0,-1.5])
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        smo = Tex(r"$f : U \subset \mathbb{R}^2 \to \mathbb{R}$ is smooth if $\frac{\partial ^{m + n} f}{\partial^m x \partial^n y}$ exist for all $m,n \in \mathbb{N}$.", font_size = 34)
        sudef = Tex(r"$\Sigma \subset \mathbb{R}^3$ is a smooth surface if for all $p \in \Sigma$ there is $p\in U \subset \mathbb{R}^3$", font_size = 34)
        sudef1 = Tex( r"open and a coordinate", font_size = 34)
        sudef2 = Tex(r"system in which $\Sigma \cap U$ is the graph of a smooth function with open domain.", font_size = 34)
        smo.to_corner(UL).shift(0.3*DOWN)
        sudef.next_to(smo, DOWN)
        sudef.shift((smo.get_left()-sudef.get_left())*RIGHT)
        sudef2.next_to(sudef,DOWN)
        sudef2.shift((smo.get_left()-sudef2.get_left())*RIGHT)
        sudef1.next_to(sudef,RIGHT)
        self.add(axes)
        self.wait()
        s1 = Surface(
            lambda u, v: axes.c2p(u-3*PI/2 , v - 3*PI/2 , -np.sin(u)/2 - np.sin(v)/2 - v / 9 +0.2 ),
            u_range = [0,3*PI],
            v_range = [0,3*PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        s2 = Surface(
            lambda u, v: axes.c2p(u-3*PI/2 , v - 3*PI/2 , -np.sin(u)/2 - np.sin(v)/2 - v / 9 +0.2 ),
            u_range = [7*PI/6, 11*PI/6],
            v_range = [3*PI/2,13*PI/6],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        self.play(Create(s1))
        self.wait()
        self.add_fixed_in_frame_mobjects(smo)
        self.remove(smo)
        self.play(Write(smo))
        self.wait()
        self.add_fixed_in_frame_mobjects(sudef)
        self.add_fixed_in_frame_mobjects(sudef1)
        self.add_fixed_in_frame_mobjects(sudef2)
        self.remove(sudef, sudef1, sudef2)
        l1 = Line([-1,1.9,0], [-0.7,1.9,0], color = ORANGE)
        l11 = Line([1.25,1.9,0], [1.55,1.9,0], color = ORANGE)
        l2 = Line([-6.7,1.9,0], [-6.3,1.9,0], color = PURPLE_B)
        l22 = Line([-0.4,1.9,0], [0,1.9,0], color = PURPLE_B)
        l3 = Line([-4.1,1.4,0], [-3.1,1.4,0], color = GREEN)
        l4 = Line([1.8,1.9,0], [2.2,1.9,0], color = BLUE)
        self.add_fixed_in_frame_mobjects(l1, l2, l3, l4, l11, l22)
        self.remove(l1,l2,l3,l4,l11,l22)
        self.play(Write(sudef), Write(sudef1), Write(sudef2), Create(l1), Create(l2), Create(l3), Create(l4), Create(l11), Create(l22))
        self.wait()
        p = Sphere(radius = 0.1)
        p.set_color(ORANGE)
        p.shift(axes.c2p(0,PI/3, 0.95 - PI *11/54 ))
        self.play(Create(p))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.wait()
        c1 = Surface(
            lambda u,v: axes.c2p(u,v + PI/3, 0.95 -PI/3 - PI *11/54 ),
            u_range = [-PI/3,PI/3],
            v_range = [-PI/3,PI/3],
            checkerboard_colors = [BLUE,BLUE],
            fill_opacity = 0.3
        )
        c2 = Surface(
            lambda u,v: axes.c2p(u,v + PI/3, 0.95 + PI/3 - PI *11/54 ),
            u_range = [-PI/3,PI/3],
            v_range = [-PI/3,PI/3],
            checkerboard_colors = [BLUE,BLUE],
            fill_opacity = 0.3
        )
        c3 = Surface(
            lambda u,v: axes.c2p(u, 2*PI/3 ,  v + 0.95 - PI *11/54 ),
            u_range = [-PI/3,PI/3],
            v_range = [-PI/3,PI/3],
            checkerboard_colors = [BLUE,BLUE],
            fill_opacity = 0.3
        )
        c4 = Surface(
            lambda u,v: axes.c2p(u, 0,  v + 0.95 - PI *11/54 ),
            u_range = [-PI/3,PI/3],
            v_range = [-PI/3,PI/3],
            checkerboard_colors = [BLUE,BLUE],
            fill_opacity = 0.3
        )
        c5 = Surface(
            lambda u,v: axes.c2p(PI/3, u+ PI/3,  v + 0.95 - PI *11/54 ),
            u_range = [-PI/3,PI/3],
            v_range = [-PI/3,PI/3],
            checkerboard_colors = [BLUE,BLUE],
            fill_opacity = 0.3
        )
        c6 = Surface(
            lambda u,v: axes.c2p(-PI/3, u+ PI/3,  v + 0.95 - PI *11/54 ),
            u_range = [-PI/3,PI/3],
            v_range = [-PI/3,PI/3],
            checkerboard_colors = [BLUE,BLUE],
            fill_opacity = 0.3
        )
        self.play(Create(c1), Create(c2), Create(c3), Create(c4), Create(c5), Create(c6))
        self.wait()
        axes1 = ThreeDAxes(x_range = [-3,3,1], y_range = [-3,3,1], z_range = [0,3,1], x_length = 6, y_length = 6, z_length = 3 )
        axes1.rotate_about_origin(-PI/6 , RIGHT)
        axes1.shift([0,0,-1.5])
        self.add(s2)
        self.play(FadeOut(s1), Create(axes1), FadeOut(axes))
        self.wait()
        self.play(FadeOut(c1), FadeOut(c2), FadeOut(c3), FadeOut(c4), FadeOut(c5), FadeOut(c6))
        self.wait()
        self.begin_ambient_camera_rotation(rate = - PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.wait()
        

class graphex(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-6,6,1], y_range = [-6,6,1], z_range = [-1,3,1], x_length = 9, y_length = 9, z_length = 3 )
        axes.shift([0,0,-1])
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        gra = Tex(r"$U \subset \mathbb{R}^2$ open, $f : U \to \mathbb{R}$ smooth $\Rightarrow $ Graph$(f)$ is a surface.", font_size = 34)
        gra.to_corner(UL).shift(0.5*DOWN)
        gra.shift( gra.get_center()[0]*LEFT )
        ul = Line([1.1,2.5,0], [2.5,2.5,0], color = PURPLE_B)
        self.add_fixed_in_frame_mobjects(gra, ul)
        self.remove(gra, ul)
        self.play(Create(axes), Write(gra), Create(ul))
        self.wait()
        nd = Surface(
            lambda u,v: axes.c2p(u, v, 2**((2 - u**2 - v**2)/2) ),
            u_range = [-5,5],
            v_range = [-5,5],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        self.play(Create(nd))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.wait()
        


class plaex(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-6,6,1], y_range = [-6,6,1], z_range = [-1,3,1], x_length = 9, y_length = 9, z_length = 3 )
        axes.shift([0,0,-1])
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        pla = Tex(r"$\Sigma = \{ ax + by + cz = d \}$ ", r" $= $ Graph $\left( (x,y) \mapsto  \frac{d- ax- by}{c} \right)$", font_size = 34)
        pla.to_corner(UL).shift(0.5*DOWN)
        pla.shift( pla.get_center()[0]*LEFT )
        pla[0].shift(LEFT)
        pla1 = MathTex(r"{{ax}} +{{ by}}  + {{c}}{{z}} = {{d}}", font_size = 34).next_to(pla[0], RIGHT).shift(3*RIGHT)
        pla2 = MathTex(r" {{c}}{{z}} = {{d}} -  {{ax}} - {{ by}}", font_size = 34).next_to(pla[0], RIGHT).shift(3*RIGHT)
        pla3 = MathTex(r"z = \frac{d - ax - by}{c}", font_size = 34).next_to(pla[0], RIGHT).shift(3*RIGHT) 
        self.add_fixed_in_frame_mobjects(pla, pla1, pla2, pla3)
        self.remove(pla, pla1,pla2, pla3)
        self.play(Create(axes), Write(pla[0]))
        self.wait()
        sig = Surface(
            lambda u,v: axes.c2p(u, v, 1 - u/5 - v/5 ),
            u_range = [-5,5],
            v_range = [-5,5],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        self.play(Create(sig))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/6)
        self.wait(2)
        self.stop_ambient_camera_rotation()
        self.wait()
        self.play(Write(pla1))
        self.wait()
        self.play(TransformMatchingShapes(pla1, pla2))
        self.wait()
        self.play(TransformMatchingShapes(pla2, pla3))
        self.wait()
        self.play(FadeOut(pla3))
        self.wait()
        self.play(Create(pla[1]), pla[0].animate.shift(RIGHT))
        self.wait()



class sphex(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-3,3,1], y_range = [-3,3,1], z_range = [-3,3,1], x_length = 6, y_length = 6, z_length = 6 )
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        s1 = Surface(
            lambda u, v: axes.c2p( 2*np.sin(u)*np.cos(v),2* np.sin(u)*np.sin(v),2* np.cos(u) ),
            u_range = [0,PI],
            v_range = [0,2*PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        s2 = Surface(
            lambda u, v: axes.c2p(2* np.sin(u)*np.cos(v),2* np.sin(u)*np.sin(v), 2*np.cos(u) ),
            u_range = [PI/8,3*PI/8],
            v_range = [PI/3,2*PI/3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        p = Sphere(radius = 0.1)
        p.set_color(ORANGE)
        p.shift([0,sqrt(2), sqrt(2)])
        axes1 = ThreeDAxes(x_range = [-3,3,1], y_range = [-3,3,1], z_range = [0,3,1], x_length = 6, y_length = 6, z_length = 3 )
        axes1.rotate_about_origin(-PI/4 , RIGHT)
        self.play(Create(axes), Create(s1))
        self.wait()
        self.play(Create(p))
        self.wait()
        self.add(s2)
        self.play(FadeOut(s1), FadeOut(axes))
        self.wait()
        self.play(Create(axes1))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.wait()
        s2s = Tex(r"$\mathbb{S}^2 $ is a surface", font_size = 40).shift(3*UP)
        self.add_fixed_in_frame_mobjects(s2s)
        self.remove(s2s)
        self.play(Write(s2s))
        self.wait()
        

class discex(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-5,5,1], y_range = [-5,5,1], z_range = [-2,2,1], x_length = 10, y_length = 10, z_length = 4 )
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        d = Surface(
            lambda u, v: axes.c2p(3* u*np.cos(v), 3* u*np.sin(v), 0 ),
            u_range = [0,1],
            v_range = [0,2*PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )      
        dt = Tex(r"$\Sigma  = \{ (x,y, z ) \in \mathbb{R}^3  \vert z= 0 , x^2 + y^2 \leq 1 \}$ is ", r"not ", r"a surface", font_size = 40).shift(3*UP)
        ul = Line([0,0,0], [0.6, 0,0], color = PURPLE_B).next_to(dt[1], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(dt, ul)
        self.remove(dt, ul)
        p = Sphere(radius = 0.1, color = ORANGE)
        p.set_color(ORANGE)
        p.shift([3/2, 3*sqrt(3)/2, 0])
        self.play(Create(axes), Create(d))
        self.wait()
        self.play(Create(p))
        self.wait()
        self.play(Write(dt), Create(ul))
        self.wait()
       

class mobex(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-5,5,1], y_range = [-5,5,1], z_range = [-2,2,1], x_length = 10, y_length = 10, z_length = 4 )
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        mb = Surface(
            lambda u, v: axes.c2p( np.cos(u)*(3 + np.cos(u/2)*v),np.sin(u)*(3 + np.cos(u/2)*v) , np.sin(u/2)*v ),
            u_range = [0,2*PI],
            v_range = [-1,1],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        mb2 = Surface(
            lambda u, v: axes.c2p( np.cos(u)*(3 + np.cos(u/2)*v),np.sin(u)*(3 + np.cos(u/2)*v) , np.sin(u/2)*v ),
            u_range = [7*PI/8,9*PI/8],
            v_range = [-1,1],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )   
        axes1 = ThreeDAxes(x_range = [-2,2,1], y_range = [-3,3,1], z_range = [0,5,1], x_length = 4, y_length = 6, z_length = 5 )
        axes1.rotate_about_origin(-PI/2 , UP)   
        mob = Tex(r"$\Sigma  =$  M\"obius band", font_size = 40).shift(3*UP)
        self.add_fixed_in_frame_mobjects(mob)
        p = Sphere(radius = 0.1, color = ORANGE)
        p.set_color(ORANGE)
        p.shift([-3, 0, 0])
        self.play(Create(axes), Create(mb), Write(mob))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/6)
        self.wait(15)
        self.stop_ambient_camera_rotation()
        self.wait()
        self.play(Create(p))
        self.wait()
        self.add(mb2)
        self.play(FadeOut(mb), Create(axes1), FadeOut(axes))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.wait()


class mobex2(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-5,5,1], y_range = [-5,5,1], z_range = [-2,2,1], x_length = 10, y_length = 10, z_length = 4 )
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        mb = Surface(
            lambda u, v: axes.c2p( np.cos(u)*(3 + np.cos(u/2)*v),np.sin(u)*(3 + np.cos(u/2)*v) , np.sin(u/2)*v ),
            u_range = [0,2*PI],
            v_range = [-1,1],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        edge = ParametricFunction(lambda t : axes.c2p( np.cos(t)*(3 + np.cos(t/2)),np.sin(t)*(3 + np.cos(t/2)) , np.sin(t/2) ), t_range = [0,4*PI], color = GREEN)
        mob = Tex(r"M\"obius band (without edges) is a surface", font_size = 40).shift(3*UP)
        self.add_fixed_in_frame_mobjects(mob) 
        self.add(axes)
        self.wait()
        self.play(Create(mb), Write(mob), Create(edge))
        self.wait()
        self.play(FadeOut(edge))
        self.wait()



class imphyp(ThreeDScene):
    def construct(self):
        axes1 = ThreeDAxes(x_range = [-3,3,1], y_range = [-3,3,1], z_range = [-3,3,1], x_length = 2, y_length = 2, z_length = 2, tips = False )
        axes2 = axes1.copy()
        axes3 = axes1.copy()
        axes1.shift([-4,4,0.5])
        axes2.shift([-2,2,-0.5])
        axes3.shift([-4,4,-2.5])
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        imp = Tex(r"If $f: \mathbb{R}^3 \to \mathbb{R}$ is smooth and $c \in \mathbb{R}$, then $\{ p \in \mathbb{R}^3 \vert f(p) = c \}$ is usually a surface.", font_size = 36).shift(3*UP)
        ex = Tex(r"Example", color = BLUE, font_size = 36).next_to(imp,DOWN).shift(5.7 *LEFT + (0.3)*DOWN)
        l = ex.get_left()[0]
        imp.shift((l - imp.get_left()[0])*RIGHT)
        ex1 = Tex(r"If $f(x,y,z) = x^2 + y^2 - z^2$,", font_size = 36).next_to(ex,DOWN)
        ex1.shift((l-ex1.get_left()[0])*RIGHT)
        ex2 = Tex(r"$\bullet$ If ", r"$c> 0 $", r", $f^{-1}(c)$ is a ", r"one sheeted hyperboloid.", font_size = 36).next_to(ex1,DOWN).shift(0.5*DOWN)
        ex2.shift((l-ex2.get_left()[0])*RIGHT)
        ex3 = Tex(r"$\bullet$ If ", r"$c< 0 $", r", $f^{-1}(c)$ is a ", r"two sheeted hyperboloid.", font_size = 36).next_to(ex2,DOWN).shift(0.5*DOWN)
        ex3.shift((l-ex3.get_left()[0])*RIGHT)
        ex4 = Tex(r"$\bullet$ If ", r"$c = 0 $", r", $f^{-1}(c)$ is a ", r"cone.", font_size = 36).next_to(ex3,DOWN).shift(0.5*DOWN)
        ex4.shift((l-ex4.get_left()[0])*RIGHT)
        ulp = Line([0,0,0], [1,0,0], color = PURPLE_B).next_to(ex2[1], DOWN).shift(0.15*UP)
        uln = Line([0,0,0], [1,0,0], color = GREEN).next_to(ex3[1], DOWN).shift(0.15*UP)
        ulz = Line([0,0,0], [1,0,0], color = BLUE).next_to(ex4[1], DOWN).shift(0.15*UP)
        ulpp = Line([0,0,0], [3,0,0], color = PURPLE_B).next_to(ex2[3], DOWN).shift(0.15*UP)
        ulnn = Line([0,0,0], [3,0,0], color = GREEN).next_to(ex3[3], DOWN).shift(0.15*UP)
        ulzz = Line([0,0,0], [0.8,0,0], color = BLUE).next_to(ex4[3], DOWN).shift(0.15*UP)
        hp = Surface(
            lambda u, v: axes1.c2p( np.cosh(u)*np.cos(v), np.cosh(u)*np.sin(v), np.sinh(u) ),
            u_range = [-1.8,1.8],
            v_range = [0,2*PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        hn1 = Surface(
            lambda u, v: axes2.c2p( np.sinh(u)*np.cos(v), np.sinh(u)*np.sin(v), np.cosh(u) ),
            u_range = [-1.8,1.8],
            v_range = [0,2*PI],
            checkerboard_colors = [GREEN,GREEN],
            fill_opacity = 0.7
        )
        hn2 = Surface(
            lambda u, v: axes2.c2p( np.sinh(u)*np.cos(v), np.sinh(u)*np.sin(v), - np.cosh(u) ),
            u_range = [-1.8,1.8],
            v_range = [0,2*PI],
            checkerboard_colors = [GREEN,GREEN],
            fill_opacity = 0.7
        )
        hz = Surface(
            lambda u, v: axes3.c2p( u * np.cos(v) , u*np.sin(v), u ),
            u_range = [-3,3],
            v_range = [0,2*PI],
            checkerboard_colors = [BLUE,BLUE],
            fill_opacity = 0.7
        )
        self.add_fixed_in_frame_mobjects(imp)
        self.remove(imp)
        self.play(Write(imp))
        self.wait()
        self.add_fixed_in_frame_mobjects(ex,ex1,ex2,ex3,ex4)
        self.remove(ex,ex1,ex2,ex3,ex4)
        self.play(Write(ex), Write(ex1))
        self.wait()
        self.add_fixed_in_frame_mobjects(ulp, uln, ulz, ulpp, ulnn, ulzz)
        self.play(Write(ex2), Write(ex3), Write(ex4), Create(axes1), Create(axes2), Create(axes3), Create(hp), Create(hn1), Create(hn2), Create(hz))
        self.wait()
        


class imp(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-4,4,1], y_range = [-4,4,1], z_range = [-4,4,1], x_length = 4, y_length = 4, z_length = 4 )
        self.set_camera_orientation(phi=75 * DEGREES, theta= -60*DEGREES)
        xl = axes.get_x_axis_label(Tex("$x$", font_size = 30))
        xl.rotate_about_origin(PI/2 , RIGHT)
        yl = axes.get_y_axis_label(Tex("$y$", font_size = 30))
        yl.rotate_about_origin(PI/2 , RIGHT)
        yl.shift([0,2,-2])
        zl = axes.get_z_axis_label(Tex("$z$", font_size = 30))
        xl.shift([3.5,2,-0.5])
        yl.shift([3.5,2,-0.5])
        zl.shift([3.5,2,-0.5])
        #zl.rotate_about_origin(PI/2 , RIGHT)
        axes.shift([3.5,2,-0.5])
        imp = Tex(r"If $f: \mathbb{R}^3 \to \mathbb{R}$ is smooth and $c \in \mathbb{R}$, then $\{ p \in \mathbb{R}^3 \vert f(p) = c \}$ is usually a surface.", font_size = 36).shift(3*UP)
        ex = Tex(r"Example", color = BLUE, font_size = 36).next_to(imp,DOWN).shift(5.7 *LEFT + (0.3)*DOWN)
        l = ex.get_left()[0]
        imp.shift((l - imp.get_left()[0])*RIGHT)
        thm = Tex(r"Theorem", color = BLUE, font_size = 36).next_to(imp,DOWN)
        thm.shift((l - thm.get_left()[0]) *RIGHT + (0.3)*DOWN)
        th1 = Tex(r"If $f : \mathbb{R}^3 \to \mathbb{R}$ is smooth and $c \in \mathbb{R}$ are such that", font_size = 36).next_to(thm,DOWN)
        th1.shift((l-th1.get_left()[0])*RIGHT)
        th2 = Tex(r"for all $p \in f^{-1}(c)$, one has $\nabla f (p) \neq 0$, then", font_size = 36).next_to(th1,DOWN)
        th2.shift((l-th2.get_left()[0])*RIGHT)
        th3 = Tex(r"each connected component of $f^{-1}(c)$ is a surface.", font_size = 36).next_to(th2,DOWN)
        th3.shift((l-th3.get_left()[0])*RIGHT)
        pr = Tex(r"Proof", color = BLUE, font_size = 36).next_to(th3,DOWN)
        pr.shift((l - pr.get_left()[0]) *RIGHT)
        pr1 = Tex(r"For $p$ with $f(p) = c$, w.l.o.g. $\frac{\partial f}{\partial x} (p)\neq 0$", font_size = 36).next_to(pr,DOWN)
        pr1.shift((l-pr1.get_left()[0])*RIGHT)
        pr2 = Tex(r"then by IFT, $f = c$ becomes $x = h(y,z)$ near $p$.", font_size = 36).next_to(pr1,DOWN)
        pr2.shift((l-pr2.get_left()[0])*RIGHT)
        sig = Surface(
            lambda u, v: axes.c2p( 3*np.cos(u)*np.cos(v),1.5* np.cos(u)*np.sin(v),1.5* np.sin(u) ),
            u_range = [-PI/2, PI/2],
            v_range = [0,2*PI],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sig2 = Surface(
            lambda u, v: axes.c2p( 3*np.cos(u)*np.cos(v),1.5* np.cos(u)*np.sin(v),1.5* np.sin(u) ),
            u_range = [-PI/8, PI/8],
            v_range = [-PI/8,PI/8],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE)
        p.shift([5,2,-0.5])
        self.add_fixed_in_frame_mobjects(imp)
        self.wait()
        self.add_fixed_in_frame_mobjects(thm,th1,th2,th3, pr, pr1, pr2)
        self.remove(thm,th1,th2,th3, pr, pr1, pr2)
        self.play(Create(axes), Write(thm), Write(th1), Write(th2), Write(th3), Create(sig), Create(xl), Create(yl), Create(zl))
        self.wait()
        self.play(Write(pr), Write(pr1), Write(pr2), Create(p))
        self.wait()
        self.add(sig2)
        self.play(FadeOut(sig))
        self.wait()
        self.play( 
            Rotate(axes, about_point = [3.5,2,-0.5], angle=-2*PI/3), 
            Rotate(p, about_point = [3.5,2,-0.5], angle=-2*PI/3), 
            Rotate(sig2, about_point = [3.5,2,-0.5], angle=-2*PI/3), 
            Rotate(xl, about_point = [3.5,2,-0.5], angle=-2*PI/3), 
            Rotate(yl, about_point = [3.5,2,-0.5], angle=-2*PI/3), 
            Rotate(zl, about_point = [3.5,2,-0.5], angle=-2*PI/3), 
            run_time = 4 
        )
        self.wait()


class ugly(ThreeDScene):
    def construct(self):
        imp = Tex(r"For $U \subset \mathbb{R}^2$ open, and $\varphi : U \to \mathbb{R}^3$ smooth, the image may not be a surface.", font_size = 36).shift(2.5*UP)
        self.add_fixed_in_frame_mobjects(imp)
        axes = ThreeDAxes(x_range = [-4,4,10], y_range = [-4,4,10], z_range = [-1,2,10], x_length = 6, y_length = 4.2, z_length = 3 )
        dom = ThreeDAxes(x_range = [-4,4,1], y_range = [-4,4,1], z_range = [-1,1,1], x_length = 4, y_length = 4, z_length = 0.01, tips = False )
        axes.shift([-sqrt(3)*0.9, 3*0.9, -0.5])
        dom.shift([sqrt(3)*0.9, -3*0.9,-0.5])
        ar = CurvedArrow(LEFT, RIGHT, radius= -2.5)
        phil = Tex(r"$\varphi$", font_size = 40).next_to(ar, UP)
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        self.play(Create(dom), Create(axes))
        self.wait()
        u = Surface(
            lambda u, v: dom.c2p( u,v,0 ),
            u_range = [-3, 3],
            v_range = [-3,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sig = Surface(
            lambda u, v: axes.c2p( 2*np.cos(2*u), 2*np.sin(2*u), u/3   ),
            u_range = [-3,3],
            v_range = [-3,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        self.play(Create(u))
        self.wait()
        self.add_fixed_in_frame_mobjects(ar, phil)
        self.remove(ar,phil)
        self.play(Create(ar), Write(phil))
        self.wait()
        self.play(Transform(u,sig))
        self.wait()



class regular(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        reg0 = Tex(r"For $U \subset \mathbb{R}^2$ open, and $\varphi : U \to \mathbb{R}^3$ smooth, we say $\varphi$ is regular if", font_size = 36).shift(3*UP)
        reg1 = Tex(r"$\left\{ \frac{\partial \varphi}{\partial u} ,\frac{\partial \varphi}{\partial v}   \right\}$ is linearly independent", font_size = 36).next_to(reg0, DOWN)
        th = Tex(r"Theorem", font_size = 36, color = BLUE).to_edge(UL).shift(1.8*DOWN + 0.2*RIGHT)
        l = th.get_left()[0]
        th1 = Tex(r"Let $U\subset \mathbb{R}^2$ be open connected and $\varphi : U \to \mathbb{R}^3$ be a", font_size = 36).next_to(th, DOWN)
        th1.shift((l - th1.get_left()[0])*RIGHT)
        th2 = Tex(r"smooth regular embedding, then $\varphi (U)$ is a surface.", font_size = 36).next_to(th1, DOWN)
        th2.shift((l - th2.get_left()[0])*RIGHT)
        pr = Tex(r"Proof", font_size = 36, color = BLUE).to_edge(UL).shift(1.85*DOWN + 0.2*RIGHT)
        pr1 = Tex(r"For $p = (x_0, y_0, z_0) = \varphi (u_0, v_0)$,", r"w.l.o.g.", font_size = 36).next_to(pr, DOWN).shift(0.05*UP)
        pr1.shift((l - pr1.get_left()[0])*RIGHT)
        pr2 = Tex(r"$  \begin{pmatrix}  \frac{\partial x }{\partial u} (u_0, v_0) & \frac{\partial x }{\partial v} (u_0, v_0) \\ \frac{\partial y }{\partial u} (u_0, v_0) & \frac{\partial y }{\partial v} (u_0, v_0) \end{pmatrix} $ is non-degenerate.", font_size = 36).next_to(pr1, DOWN)
        pr2.shift((l - pr2.get_left()[0])*RIGHT + 0.05*UP) 
        pr3 = Tex(r"By IFV, there are $V_1 , V_2 \subset \mathbb{R}^2$ open s.t. $(u_0, v_0) \in V_1$, $(x_0, y_0 ) \in V_2$, and", font_size = 36).to_edge(UL).shift(0.8*DOWN)
        pr3.shift((l - pr3.get_left()[0])*RIGHT)
        pr4 = Tex(r"proj $ \circ \varphi : V_1 \to V_2$ is a homeomorphism with smooth inverse $\psi : V_2 \to V_1$.", font_size = 36).next_to(pr3, DOWN)
        pr4.shift((l - pr4.get_left()[0])*RIGHT) 
        pr5 = Tex(r" Then $\varphi \circ \psi : V_2 \to \mathbb{R}^3$ respects $xy$-coordinates." , font_size = 36).next_to(pr4, DOWN)
        pr5.shift((l - pr5.get_left()[0])*RIGHT) 
        axes = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-2,3,10], x_length = 4.5, y_length = 4.5, z_length = 3.75 )
        dom = ThreeDAxes(x_range = [-2.5,2.5,1], y_range = [-4,4,1], z_range = [-1,1,2], x_length = 3.75, y_length = 6, z_length = 0.01, tips = False )
        axes.shift([-sqrt(3), 3, -1.5])
        dom.shift([sqrt(3), -3,-2.5])
        e1 = Arrow(dom.c2p( 0,PI/6,0 ), dom.c2p( 2,PI/6,0 ), color = YELLOW)
        e2 = Arrow(dom.c2p( 0,PI/6,0 ), dom.c2p( 0,PI/6 + 2,0 ), color = YELLOW)
        par1 = Arrow(axes.c2p( 0,1,sqrt(3) ), axes.c2p( 2,1,sqrt(3) ), color = YELLOW)
        par2 = Arrow(axes.c2p( 0,1,sqrt(3) ), axes.c2p( 0, 1 + 2*sqrt(3),sqrt(3) - 1 ), color = YELLOW)
        par11 = par1.copy()
        par22 = par2.copy()
        ppar1 = Arrow(axes.c2p( 0,1,0) , axes.c2p( 2,1,0 ), color = YELLOW)
        ppar2 = Arrow(axes.c2p( 0,1,0), axes.c2p( 0, 1 + 2*sqrt(3), 0 ), color = YELLOW)
        projp = DashedLine( axes.c2p( 0,1,sqrt(3) ),axes.c2p( 0,1,0 ) )
        proju = DashedLine( axes.c2p( 1.7,1,sqrt(3) ),axes.c2p( 1.7,1,0 ) )
        projv = DashedLine( axes.c2p( 0, 1 + 2*sqrt(3)*0.9  ,sqrt(3) - 0.9 ),axes.c2p( 0, 1 + 2*sqrt(3)*0.9, 0  ) )
        ar = Arrow(LEFT + 1.5*DOWN, 1.4* RIGHT + 0.9*DOWN, stroke_width=4, max_tip_length_to_length_ratio=0.2)
        phil = Tex(r"$\varphi$", font_size = 40).next_to(ar, UP).shift(0.2*DOWN )
        ari = CurvedArrow(  1.6* RIGHT + 2.3*DOWN ,  LEFT + 2.5*DOWN  , radius = -7)
        psil = Tex(r"$\psi$", font_size = 40).next_to(ari, DOWN).shift(0.2*UP+0.2*RIGHT)
        u = Surface(
            lambda u, v: dom.c2p( u,v - PI,0 ),
            u_range = [-2,2],
            v_range = [PI/6,11*PI/6],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        u0 = u.copy()
        u1 = Surface(
            lambda u, v: dom.c2p( u,v - PI,0 ),
            u_range = [-0.5, 0.5],
            v_range = [PI,4*PI/3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        u10 = u1.copy()
        sig = Surface(
            lambda u, v: axes.c2p( u , - (2- (u**2)/4)*np.sin(v), -(2- (u**2)/4)*np.cos(v)),
            u_range = [-2,2],
            v_range = [PI/6,11*PI/6],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sig1 = Surface(
            lambda u, v: axes.c2p( u , - (2- (u**2)/4)*np.sin(v), -(2- (u**2)/4)*np.cos(v)),
            u_range = [-0.5, 0.5],
            v_range = [PI,4*PI/3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        sig1p = Surface(
            lambda u, v: axes.c2p( u , - (2- (u**2)/4)*np.sin(v), 0 ),
            u_range = [-0.5, 0.5],
            v_range = [PI,4*PI/3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        q = Dot( dom.c2p( 0,PI/6,0 ), color = ORANGE )
        p = Sphere( radius = 0.05 )
        p.set_color(ORANGE)
        p.shift(axes.c2p( 0, 1, sqrt(3) ))
        p0 = p.copy()
        pp = Sphere( radius = 0.05 )
        pp.set_color(ORANGE)
        pp.shift(axes.c2p( 0, 1, 0 ))
        self.add_fixed_in_frame_mobjects(reg0, reg1, th, th1, th2, pr, pr1, pr2, ar, phil, pr3, pr4, pr5, ari, psil)
        self.remove(reg0, reg1, th, th1, th2, pr, pr1, pr2, ar, phil, pr3, pr4, pr5, ari, psil)
        self.play(Write(reg0), Write(reg1))
        self.wait()
        self.play(Write(th), Write(th1), Write(th2))
        self.wait()
        self.play(th.animate.shift(1.6*UP), th1.animate.shift(1.65*UP), th2.animate.shift(1.7*UP), FadeOut(reg0), FadeOut(reg1))
        self.wait()
        self.play(Write(pr), Create(axes), Create(dom))
        self.wait()
        self.play(Create(u))
        self.wait()
        self.add(u0)
        self.play(Transform(u0,sig), Create(ar), Create(phil))
        self.wait()
        self.play(Create(p), Create(q), Write(pr1[0]))
        self.wait()
        self.play( Create(e1), Create(e2), Create(par1), Create(par2))
        self.wait()
        self.add(par11, par22, p0)
        self.play(Write(pr1[1]), Write(pr2), Transform(par11, ppar1), Transform(par22, ppar2), Create(projp), Create(proju), Create(projv), Transform(p0,pp))
        self.wait()
        self.play( FadeOut(th), FadeOut(th1), FadeOut(th2), FadeOut(pr), FadeOut(pr1), FadeOut(pr2))
        self.wait()
        self.add(u1)
        self.play(FadeOut(u), Write(pr3), Write(pr4), FadeOut(e1), FadeOut(e2), FadeOut(par1), FadeOut(par2), FadeOut(par11), FadeOut(par22), FadeOut(projp), FadeOut(proju), FadeOut(projv), FadeOut(u0))
        self.wait()
        self.add(u10)
        self.play(Transform(u10, sig1))
        self.wait()
        self.play(Transform(u10, sig1p))
        self.wait()
        self.add(sig1p)
        self.play(Transform(sig1p,u1), Write(pr5), Write(psil), Create(ari))
        self.play(Transform(sig1p, sig1))
        self.wait()


class charts(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        th = Tex(r"Theorem", font_size = 36, color = BLUE).to_edge(UL).shift(0.2*RIGHT)
        l = th.get_left()[0]
        th1 = Tex(r"Let $U\subset \mathbb{R}^2$ be open connected and $\varphi : U \to \mathbb{R}^3$ be a", font_size = 36).next_to(th, DOWN).shift(0.1*UP)
        th1.shift((l - th1.get_left()[0])*RIGHT)
        th2 = Tex(r"smooth regular embedding, then $\varphi (U)$ is a surface.", font_size = 36).next_to(th1, DOWN).shift(0.1*UP)
        th2.shift((l - th2.get_left()[0])*RIGHT)
        co = Tex(r"Corollary", font_size = 36, color = BLUE).next_to(th2, DOWN)
        co.shift((l - co.get_left()[0])*RIGHT)
        co1 = Tex(r"A connected subset ", r"$\Sigma$", r" $ \subset \mathbb{R}^3$ is a surface if and only if for each", font_size = 36).next_to(co, DOWN).shift(0.1*UP)
        co1.shift((l - co1.get_left()[0])*RIGHT)      
        co2 = Tex(r"$p$", r" $\in \Sigma$ there are $ p \in$ ", r"$ U $", r" $ \subset \mathbb{R}^3$, ", r"$W $ ", r"$ \subset \mathbb{R}^2$ open, and a smooth", font_size = 36).next_to(co1, DOWN).shift(0.1*UP)
        co2.shift((l - co2.get_left()[0])*RIGHT)      
        co3 = Tex(r"regular embedding $\varphi : W \to \mathbb{R}^3$ with $\varphi (W)  = \Sigma \cap U$.", font_size = 36).next_to(co2, DOWN).shift(0.1*UP)
        co3.shift((l - co3.get_left()[0])*RIGHT)   
        axes = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-2,3,10], x_length = 4.5, y_length = 4.5, z_length = 3.75 )
        dom = ThreeDAxes(x_range = [-2.5,2.5,10], y_range = [-4,4,10], z_range = [-1,1,2], x_length = 3.75, y_length = 6, z_length = 0.01, tips = False )
        axes.shift([-sqrt(3), 3, -2.5])
        dom.shift([sqrt(3), -3,-2.5])
        uls = Line([0,0,0], [0.5, 0,0], color = PURPLE_B).next_to(co1[1], DOWN).shift(0.15*UP)
        ulp = Line([0,0,0], [0.4, 0,0], color = ORANGE).next_to(co2[0], DOWN).shift(0.15*UP)
        ulu = Line([0,0,0], [0.5, 0,0], color = BLUE).next_to(co2[2], DOWN).shift(0.15*UP)
        ulw = Line([0,0,0], [0.5, 0,0], color = GREEN).next_to(co2[4], DOWN).shift(0.15*UP)
        u1 = Surface(
            lambda u, v: dom.c2p( 2*u,3*(v - PI),0 ),
            u_range = [-0.5, 0.5],
            v_range = [7*PI/8,9*PI/8],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        u10 = u1.copy()
        sig = Surface(
            lambda u, v: axes.c2p( u , - (2- (u**2)/4)*np.sin(v), -(2- (u**2)/4)*np.cos(v)),
            u_range = [-2,2],
            v_range = [PI/6,11*PI/6],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sig1 = Surface(
            lambda u, v: axes.c2p( u , - (2- (u**2)/4)*np.sin(v), -(2- (u**2)/4)*np.cos(v)),
            u_range = [-0.5, 0.5],
            v_range = [7*PI/8,9*PI/8],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        b1 = Surface(
            lambda u, v: axes.c2p( u , - (2- (u**2)/4)*np.sin(v), -(2- (u**2)/4)*np.cos(v) + 1/2 ),
            u_range = [-0.5, 0.5],
            v_range = [7*PI/8,9*PI/8],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        b2 = Surface(
            lambda u, v: axes.c2p( u , - (2- (u**2)/4)*np.sin(v), -(2- (u**2)/4)*np.cos(v) - 1/2 ),
            u_range = [-0.5, 0.5],
            v_range = [7*PI/8,9*PI/8],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        b3 = Surface(
            lambda u, v: axes.c2p( 1/2 , - ( 31/16)*np.sin(v), -(31/16)*np.cos(v) +u ),
            u_range = [-0.5, 0.5],
            v_range = [7*PI/8,9*PI/8],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        b4 = Surface(
            lambda u, v: axes.c2p( -1/2 , - ( 31/16)*np.sin(v), -(31/16)*np.cos(v) +u ),
            u_range = [-0.5, 0.5],
            v_range = [7*PI/8,9*PI/8],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        b5 = Surface(
            lambda u, v: axes.c2p( u , - (2- (u**2)/4)*np.sin(PI/8), (2- (u**2)/4)*np.cos(PI/8) +v ),
            u_range = [-0.5, 0.5],
            v_range = [-0.5, 0.5],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        b6 = Surface(
            lambda u, v: axes.c2p( u ,  (2- (u**2)/4)*np.sin(PI/8), (2- (u**2)/4)*np.cos(PI/8) +v ),
            u_range = [-0.5, 0.5],
            v_range = [-0.5, 0.5],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        p = Sphere( radius = 0.05 )
        p.set_color(ORANGE)
        p.shift(axes.c2p( 0, 0 , 2 ))
        ar = CurvedArrow(1.2*LEFT + 2*DOWN, 1.2* RIGHT + 1.2*DOWN, radius = -7)
        phil = Tex(r"$\varphi$", font_size = 40).next_to(ar, UP)
        self.add_fixed_in_frame_mobjects(th, th1, th2, co, co1, co2, co3, ar, phil, uls, ulp, ulu, ulw)
        self.remove(th, th1, th2, co, co1, co2, co3, ar, phil, uls, ulp, ulu, ulw)
        self.play(Write(th), Write(th1), Write(th2))
        self.wait()
        self.play(Write(co), Write(co1), Write(co2), Write(co3), Create(uls), Create(ulp), Create(uls), Create(ulu), Create(ulw))
        self.wait()
        self.play(Create(axes), Create(dom), Create(sig))
        self.wait()
        self.play(Create(p))
        self.wait()
        self.play(Create(u1), Create(b1), Create(b2), Create(b3), Create(b4), Create(b5), Create(b6))
        self.wait()
        self.add(sig1)
        self.play(FadeOut(sig))
        self.wait()
        self.play(FadeOut(b1), FadeOut(b2), FadeOut(b3), FadeOut(b4), FadeOut(b5), FadeOut(b6))
        self.wait()
        self.add(u10)
        self.play(Transform(u10, sig1), Create(ar), Write(phil))
        self.wait()


class cgraph(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        axes = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-1,3,10], x_length = 5, y_length = 5, z_length = 20/6 )
        dom = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-1,1,2], x_length = 5, y_length = 5, z_length = 0.01, tips = False )
        axes.shift([-sqrt(3), 3, -1.5])
        dom.shift([sqrt(3), -3,-1.5])
        u = Surface(
            lambda u, v: dom.c2p( u,v,0 ),
            u_range = [-2.5, 2.5],
            v_range = [-2.5,2.5],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7
        )
        u0 = u.copy()
        sig = Surface(
            lambda u, v: axes.c2p( u,v, 1 + np.sin(2*u)/3 + np.sin(2*v)/3 ),
            u_range = [-2.5, 2.5],
            v_range = [-2.5,2.5],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7
        )
        cap = Tex(r"If $\Sigma$ is the graph of $f: U \to \mathbb{R}$, then one can use", font_size = 36).to_edge(UL).shift(RIGHT + 0.2*DOWN)
        l = cap.get_left()[0]
        cap1 = Tex(r"$\varphi : U \to \Sigma$ given by  $\varphi (u,v) = (u,v, f(u,v))$ as a chart.", font_size = 36).next_to(cap, DOWN)
        cap1.shift((l - cap1.get_left()[0])*RIGHT)
        self.add_fixed_in_frame_mobjects(cap, cap1)
        self.remove(cap, cap1)
        self.play(Write(cap), Write(cap1), Create(axes), Create(dom), Create(sig))  
        self.wait()
        self.play(Create(u))
        self.wait()
        self.add(u0)
        self.play(Transform(u0, sig))
        self.wait()    
        self.remove(u0)
        self.play(axes.animate.shift([sqrt(3), -3, 0]) , sig.animate.shift([sqrt(3), -3, 0]), dom.animate.shift([-sqrt(3), 3, 0]), u.animate.shift([-sqrt(3), 3, 0]) )
        self.remove(dom)
        self.wait()
        u00 = u.copy()
        self.add(u00)
        self.play(Transform(u00, sig))
        self.wait()  
        


class csph(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        axes = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-1.5,1.5,10], x_length = 5, y_length = 5, z_length = 5 )
        axes.shift([0,0, -1.6])
        sig = Surface(
            lambda u, v: axes.c2p( np.sin(u)*np.cos(v), np.sin(u)*np.sin(v), np.cos(u) ),
            u_range = [0, PI],
            v_range = [0, 2*PI],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7
        )
        u = Surface(
            lambda u, v: axes.c2p(np.sin(u)* np.cos(v), np.sin(u)*np.sin(v), 0 ),
            u_range = [0,PI/2],
            v_range = [0, 2*PI],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        u0 = u.copy()
        u00 = Surface(
            lambda u, v: axes.c2p(np.sin(u)* np.cos(v), np.sin(u)*np.sin(v), 0 ),
            u_range = [0,PI/2],
            v_range = [0, 2*PI],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        s1 = Surface(
            lambda u, v: axes.c2p(np.sin(u)* np.cos(v), np.sin(u)*np.sin(v), np.cos(u) ),
            u_range = [0,PI/2],
            v_range = [0, 2*PI],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        s2 = Surface(
            lambda u, v: axes.c2p(np.sin(u)* np.cos(v), np.sin(u)*np.sin(v), - np.cos(u) ),
            u_range = [0,PI/2],
            v_range = [0, 2*PI],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        cap = Tex(r"For $\Sigma = \mathbb{S}^2$, the chart $\varphi_1 : U \to \Sigma$ given by", font_size = 36).to_edge(UL).shift( 0.2*DOWN)
        l = cap.get_left()[0]
        cap1 = Tex(r"$\varphi_1 (u,v) = (u,v, \sqrt{1 - u^2 - v^2 } )$", r" covers the upper hemisphere.", font_size = 36).next_to(cap, DOWN)
        cap1.shift((l - cap1.get_left()[0])*RIGHT)
        cap2 = Tex(r"$\varphi_2 (u,v) = (u,v,  - \sqrt{1 - u^2 - v^2 } )$", r" covers the lower hemisphere.", font_size = 36).next_to(cap1, DOWN)
        cap2.shift((l - cap2.get_left()[0])*RIGHT)
        ul1 = Line([0,0,0], [3,0,0], color = GREEN).next_to(cap1[0], DOWN).shift(0.15*UP)
        ul2 = Line([0,0,0], [3,0,0], color = BLUE).next_to(cap2[0], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(cap, cap1, cap2, ul1, ul2)
        self.remove(cap, cap1, cap2, ul1, ul2)
        self.play(Write(cap), Write(cap1), Create(axes), Create(sig), Create(ul1))
        self.wait()
        self.play(Create(u), FadeOut(sig))
        self.wait()
        self.play(Transform(u, s1))
        self.wait()
        self.play(FadeOut(u))
        self.wait()
        self.play(Create(u00), Create(cap2), Create(ul2))
        self.wait()
        self.play(Transform(u00, s2))
        self.wait()
        self.play(FadeOut(u00), FadeOut(cap), FadeOut(cap1), FadeOut(cap2))
        self.wait()



class csph2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        axes = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-1.5,1.5,10], x_length = 5, y_length = 5, z_length = 5 )
        axes.shift([0,0, -1.6])
        sig = Surface(
            lambda u, v: axes.c2p( np.sin(u)*np.cos(v), np.sin(u)*np.sin(v), np.cos(u) ),
            u_range = [0, PI],
            v_range = [0, 2*PI],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7
        )
        u = Surface(
            lambda u, v: axes.c2p(np.sin(u)* np.cos(v),0,  np.sin(u)*np.sin(v) ),
            u_range = [0,PI/2],
            v_range = [0, 2*PI],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        u00 = Surface(
            lambda u, v: axes.c2p(np.sin(u)* np.cos(v),0,  np.sin(u)*np.sin(v) ),
            u_range = [0,PI/2],
            v_range = [0, 2*PI],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        s1 = Surface(
            lambda u, v: axes.c2p(np.sin(u)* np.cos(v),np.cos(u), np.sin(u)*np.sin(v) ),
            u_range = [0,PI/2],
            v_range = [0, 2*PI],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        s2 = Surface(
            lambda u, v: axes.c2p(np.sin(u)* np.cos(v),-np.cos(u),  np.sin(u)*np.sin(v)),
            u_range = [0,PI/2],
            v_range = [0, 2*PI],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        cap = Tex(r"$\varphi_3 (u,v) = (u, \sqrt{1 - u^2 - v^2 } , v)$", r" covers the right hemisphere.", font_size = 36).to_edge(UL).shift( 0.2*DOWN)
        l = cap.get_left()[0]
        cap1 = Tex(r"$\varphi_4 (u,v) = (u, -  \sqrt{1 - u^2 - v^2 },v )$", r" covers the left hemisphere.", font_size = 36).next_to(cap, DOWN)
        cap1.shift((l - cap1.get_left()[0])*RIGHT)
        ul1 = Line([0,0,0], [3,0,0], color = GREEN).next_to(cap[0], DOWN).shift(0.15*UP)
        ul2 = Line([0,0,0], [3,0,0], color = BLUE).next_to(cap1[0], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(cap, cap1, ul1, ul2)
        self.remove(cap, cap1, ul1, ul2)
        self.add(axes)
        self.play(Write(cap), Write(cap1), Create(u), Create(ul1), Create(ul2))
        self.wait()
        self.play(Transform(u, s1))
        self.wait()
        self.play(FadeOut(u))
        self.wait()
        self.play(Create(u00))
        self.wait()
        self.play(Transform(u00, s2))
        self.wait()


class ste(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        axes = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-1.5,1.5,10], x_length = 6, y_length = 6, z_length = 6 )
        axes.shift([0,0, -1.5])
        sig = Surface(
            lambda u, v: axes.c2p( np.sin(u)*np.cos(v), np.sin(u)*np.sin(v), np.cos(u) ),
            u_range = [0, PI],
            v_range = [0, 2*PI],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.3
        )
        eq = ParametricFunction(lambda t : axes.c2p(np.cos(t) ,np.sin(t) ,0 ), t_range = [0,2*PI])
        sou = Sphere(radius = 0.1).set_color(YELLOW).shift(axes.c2p(0,0,-1))
        lin = DashedLine(axes.c2p(0,0,-1), axes.c2p(0,1,1) )
        p = Sphere(radius = 0.1)
        p.set_color(ORANGE)
        p.shift(axes.c2p( 0, 0.5, 0))
        q = Sphere(radius = 0.1)
        q.set_color(GREEN)
        q.shift(axes.c2p( 0, 4/5, 3/5))
        def g(t):
            return [t*np.cos(2*PI*t - PI/2), t*np.sin(2*PI*t - PI/2)   ,0] 
        def g2(t):
            nor = g(t)[0]**2 + g(t)[1]**2
            return [  2*g(t)[0]/(nor + 1) , 2*g(t)[1]/(nor + 1) , (1 - nor) / (nor + 1) ] 
        z = ValueTracker(1/2)
        p.add_updater( lambda x: x.become( Sphere(radius = 0.1).set_color(ORANGE).shift(axes.c2p(  g(z.get_value())[0],g(z.get_value())[1] , 0 )))  )
        q.add_updater( lambda x: x.become( Sphere(radius = 0.1).set_color(GREEN).shift(axes.c2p(  g2(z.get_value())[0], g2(z.get_value())[1], g2(z.get_value())[2] )))  )
        lin.add_updater( lambda x: x.become( DashedLine( axes.c2p(0,0,-1), axes.c2p(2*g(z.get_value())[0], 2* g(z.get_value())[1] ,1) ))  )
        cap = Tex(r"For ", r"$p$ ", r"$ \in \mathbb{R}^2$, ", r"$\varphi (p)$ ", r"$ \in \mathbb{S}^2$ is the point collinear with $p$ and the ", r"south pole.", font_size = 36).to_edge(UL).shift( 0.2*DOWN)
        cap.shift(( cap.get_center()[0])*LEFT)
        cap2 = Tex(r"$\varphi (u,v) = \left( \frac{2u}{u^2 + v^2 + 1}, \frac{2v}{u^2 + v^2 + 1}, \frac{1 - u^2 - v^2}{u^2 + v^2 + 1}    \right).$", font_size = 36).next_to(cap, DOWN)
        cap2.shift(( cap2.get_center()[0])*LEFT)
        ulp = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(cap[1], DOWN).shift(0.15*UP)
        ulpp = Line([0,0,0], [0.5,0,0], color = GREEN).next_to(cap[3], DOWN).shift(0.15*UP)
        ulsp = Line([0,0,0], [1.5,0,0], color = YELLOW).next_to(cap[5], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(cap, cap2, ulp, ulpp, ulsp)
        self.remove(cap, cap2, ulp, ulpp, ulsp)
        self.play(Create(eq), Create(axes), Create(sig), Create(sou))
        self.wait()
        self.play(Create(p), Write(cap), Create(ulp), Create(ulpp), Create(ulsp))
        self.wait()
        self.play(Create(lin), Create(q))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.wait()
        self.play(z.animate.set_value(3/2), run_time = 6)
        self.wait()
        self.play( Write(cap2))
        self.wait()
        self.begin_ambient_camera_rotation(rate = -PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.wait()


class ste2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        axes = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-1.5,1.5,10], x_length = 6, y_length = 6, z_length = 6 )
        axes.shift([0,0, -1.5])
        sig = Surface(
            lambda u, v: axes.c2p( np.sin(u)*np.cos(v), np.sin(u)*np.sin(v), np.cos(u) ),
            u_range = [0, PI],
            v_range = [0, 2*PI],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.3
        )
        eq = ParametricFunction(lambda t : axes.c2p(np.cos(t) ,np.sin(t) ,0 ), t_range = [0,2*PI])
        nor = Sphere(radius = 0.1).set_color(YELLOW).shift(axes.c2p(0,0,1))
        lin = DashedLine(axes.c2p(0,0,1), axes.c2p(0,1,-1) )
        p = Sphere(radius = 0.1)
        p.set_color(ORANGE)
        p.shift(axes.c2p( 0, 0.5, 0))
        q = Sphere(radius = 0.1)
        q.set_color(GREEN)
        q.shift(axes.c2p( 0, 4/5, - 3/5))
        def g(t):
            return [t*np.cos(2*PI*t - PI/2), t*np.sin(2*PI*t - PI/2)   ,0] 
        def g2(t):
            nor = g(t)[0]**2 + g(t)[1]**2
            return [  2*g(t)[0]/(nor + 1) , 2*g(t)[1]/(nor + 1) , ( nor - 1) / (nor + 1) ] 
        z = ValueTracker(1/2)
        p.add_updater( lambda x: x.become( Sphere(radius = 0.1).set_color(ORANGE).shift(axes.c2p(  g(z.get_value())[0],g(z.get_value())[1] , 0 )))  )
        q.add_updater( lambda x: x.become( Sphere(radius = 0.1).set_color(GREEN).shift(axes.c2p(  g2(z.get_value())[0], g2(z.get_value())[1], g2(z.get_value())[2] )))  )
        lin.add_updater( lambda x: x.become( DashedLine( axes.c2p(0,0,1), axes.c2p(2*g(z.get_value())[0], 2* g(z.get_value())[1] ,- 1) ))  )
        cap = Tex(r"For ", r"$p$ ", r"$ \in \mathbb{R}^2$, ", r"$\psi (p)$ ", r"$ \in \mathbb{S}^2$ is the point collinear with $p$ and the ", r"north pole.", font_size = 36).to_edge(UL).shift( 0.2*DOWN)
        cap.shift(( cap.get_center()[0])*LEFT)
        cap2 = Tex(r"$\psi (u,v) = \left( \frac{2u}{u^2 + v^2 + 1}, \frac{2v}{u^2 + v^2 + 1}, \frac{ u^2 + v^2 - 1}{u^2 + v^2 + 1}    \right).$", font_size = 36).next_to(cap, DOWN)
        cap2.shift(( cap2.get_center()[0])*LEFT)
        ulp = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(cap[1], DOWN).shift(0.15*UP)
        ulpp = Line([0,0,0], [0.5,0,0], color = GREEN).next_to(cap[3], DOWN).shift(0.15*UP)
        ulsp = Line([0,0,0], [1.5,0,0], color = YELLOW).next_to(cap[5], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(cap, cap2, ulp, ulpp, ulsp)
        self.remove(cap, cap2, ulp, ulpp, ulsp)
        self.play(Create(eq), Create(axes), Create(sig), Create(nor))
        self.wait()
        self.play(Create(p), Write(cap), Create(ulp), Create(ulpp), Create(ulsp), Write(cap2))
        self.wait()
        self.play(Create(lin), Create(q))
        self.wait()
        self.play(z.animate.set_value(3/2), run_time = 6)
        self.wait()


class cmob(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-5,5,1], y_range = [-5,5,1], z_range = [-2,2,1], x_length = 5, y_length = 5, z_length = 2 )
        dom = ThreeDAxes(x_range = [-2,2,10], y_range = [-5,5,10], z_range = [-1,1,2], x_length = 2, y_length = 5, z_length = 0.01, tips = False )
        axes.shift([-sqrt(3), 3, -1.5])
        dom.shift([sqrt(3), -3,-1.5])
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        mb = Surface(
            lambda u, v: axes.c2p( np.cos(u)*(3 + np.cos(u/2)*v),np.sin(u)*(3 + np.cos(u/2)*v) , np.sin(u/2)*v ),
            u_range = [0,2*PI],
            v_range = [-1,1],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        m1 = Surface(
            lambda u, v: axes.c2p( np.cos(u)*(3 + np.cos(u/2)*v),np.sin(u)*(3 + np.cos(u/2)*v) , np.sin(u/2)*v ),
            u_range = [-5*PI / 4,PI/4],
            v_range = [-1,1],
            checkerboard_colors = [GREEN,GREEN],
            fill_opacity = 0.7
        )
        m2 = Surface(
            lambda u, v: axes.c2p( np.cos(u)*(3 + np.cos(u/2)*v),np.sin(u)*(3 + np.cos(u/2)*v) , np.sin(u/2)*v ),
            u_range = [-PI/4, 5*PI/4],
            v_range = [-1,1],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7
        )
        u1 = Surface(
            lambda u, v: dom.c2p( v , u ,0 ),
            u_range = [-5*PI / 4,PI/4],
            v_range = [-1,1],
            checkerboard_colors = [GREEN,GREEN],
            fill_opacity = 0.7
        )
        u2 = Surface(
            lambda u, v: dom.c2p( v , u , 0 ),
            u_range = [-PI/4, 5*PI/4],
            v_range = [-1,1],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7
        )
        u11 = u1.copy()
        u22 = u2.copy()
        mob = Tex(r"$\Sigma  =$  M\"obius band", font_size = 40).to_edge(UL)
        mob1 = Tex(r"$\varphi_1 : [-1,1]\times \left[ -\frac{5 \pi}{4}, \frac{\pi}{4} \right] \to \Sigma$, ", r" $\varphi_2 : [-1,1]\times \left[ -\frac{ \pi}{4}, \frac{5 \pi}{4} \right] \to \Sigma$", font_size = 40).next_to(mob, DOWN)
        mob2 = Tex(r"$\varphi _ 1 (u,v) $", r" $=  \left( \cos (v) (3 + \cos (\frac{v}{2} ) u) , \sin (v) (3 + \cos (\frac{v}{2} ) u) , \sin (\frac{v}{2}) u \right).$", font_size = 40).next_to(mob1, DOWN)
        mob2.shift((mob2.get_center()[0])*LEFT)
        mob3 = Tex(r"$\varphi _ 2 (u,v) $", r" $= \left( \cos (v) (3 + \cos (\frac{v}{2} ) u) , \sin (v) (3 + \cos (\frac{v}{2} ) u) , \sin (\frac{v}{2}) u \right).$", font_size = 40).next_to(mob2, DOWN)
        mob3.shift((mob3.get_center()[0])*LEFT)
        mob.shift((mob3.get_left()[0] - mob.get_left()[0])*RIGHT)
        mob1.shift((mob3.get_left()[0] - mob1.get_left()[0])*RIGHT)
        ul0 = Line([0,0,0], [3,0,0], color = PURPLE_B).next_to(mob, DOWN).shift(0.15*UP)
        ul1 = Line([0,0,0], [5,0,0], color = GREEN).next_to(mob1[0], DOWN).shift(0.15*UP)
        ul2 = Line([0,0,0], [5,0,0], color = ORANGE).next_to(mob1[1], DOWN).shift(0.15*UP)
        ul11 = Line([0,0,0], [1.5,0,0], color = GREEN).next_to(mob2[0], DOWN).shift(0.15*UP)
        ul22 = Line([0,0,0], [1.5,0,0], color = ORANGE).next_to(mob3[0], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(mob, mob1, mob2, mob3, ul1, ul2, ul0, ul11, ul22)
        self.remove(mob, mob1, mob2, mob3, ul0, ul1, ul2, ul11, ul22)
        self.play(Create(axes), Create(dom), Create(mb), Write(mob), Write(mob1), Write(mob2), Write(mob3), Create(ul0), Create(ul1), Create(ul2), Create(ul11), Create(ul22))
        self.wait()
        self.play(Create(u1), FadeOut(mb))
        self.wait()
        self.add(u11)
        self.play(Transform(u11, m1))
        self.wait()
        self.play(FadeOut(u1), FadeOut(u11), FadeIn(u2))
        self.wait()
        self.add(u22)
        self.play(Transform(u22, m2))
        self.wait()













class tracktest(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        axes = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-1.5,1.5,10], x_length = 6, y_length = 6, z_length = 6 )
        p = Sphere(radius = 0.1).shift(axes.c2p(0,1/2,0)).set_color(ORANGE)
        self.play(Create(p), Create(axes))
        self.wait()
        z = ValueTracker(1/2)
        p.add_updater( lambda x: x.become( Sphere( radius = 0.1).set_color(ORANGE).shift(axes.c2p(np.sin(z.get_value()-1/2) ,z.get_value(),0))))
        self.play(z.animate.set_value(PI + 1/2), run_time = 2)
        self.wait()

class uglytest(ThreeDScene):
    def construct(self):
        imp = Tex(r"For $U \subset \mathbb{R}^2$ open, and $\varphi : U \to \mathbb{R}^3$ smooth, the image may not be a surface.", font_size = 36).shift(3.5*UP)
        self.add_fixed_in_frame_mobjects(imp)


class raptest(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-4,4,1], y_range = [-4,4,1], z_range = [-4,4,1], x_length = 4, y_length = 4, z_length = 4)
        self.set_camera_orientation(phi=75 * DEGREES, theta= -60*DEGREES)
        axes.shift([0,4,0])
        p = Sphere(radius = 0.1)
        self.play(Create(axes), Create(p))
        self.wait()
        self.play(
            Rotate(axes, about_point = [0,4,0], angle=PI/2), Rotate(p, about_point = [0,4,0], angle=PI/2), run_time=2 
        )
        self.wait()

class axtest(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes(x_range = [-4,4,1], y_range = [-4,4,1], z_range = [-4,4,1], x_length = 4, y_length = 4, z_length = 4)
        label = axes.get_x_axis_label(Tex("$x$", font_size = 30))
        label.rotate_about_origin(PI/2 , RIGHT)
        self.set_camera_orientation(phi=75 * DEGREES, theta= -60*DEGREES)
        self.play(Create(axes), Create(label))
        self.wait()

class eqtest(Scene):
    def construct(self):
        eq2 = MathTex(r" {{ a }} = {{ c }} - {{ b }}")
        eq3 = MathTex(r"a = c - 38x + 21 - {{b}}")
        self.add(eq2)
        self.wait()
        self.play(TransformMatchingShapes(eq2, eq3))
        self.wait()


class txttest(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        smo = Tex(r"Sample text.", font_size = 34)
        smo.to_corner(UL).shift(0.3*DOWN)
        self.add_fixed_in_frame_mobjects(smo)
        axes1 = ThreeDAxes(x_range = [-3,3,1], y_range = [-3,3,1], z_range = [0,3,1], x_length = 6, y_length = 6, z_length = 3 )
        self.play(Create(axes1))
        self.wait()
        self.remove(smo)
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/6)
        self.wait(2)
        self.stop_ambient_camera_rotation()
        self.play(Create(smo))
        self.wait()



class ultest(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        smo = Tex(r"$f : U \subset \mathbb{R}^2 \to \mathbb{R}$ is smooth if $\frac{\partial ^{m + n} f}{\partial^m x \partial^n y}$ exist for all $m,n \in \mathbb{N}$.", font_size = 34)
        sudef = Tex(r"$\Sigma \subset \mathbb{R}^3$ is a smooth surface if for all $p \in \Sigma$ there is $p \in U \subset \mathbb{R}^3$", font_size = 34)
        sudef1 = Tex( r"open and a coordinate", font_size = 34)
        sudef2 = Tex(r"system in which $\Sigma \cap U$ is the graph of a smooth function with open domain.", font_size = 34)
        smo.to_corner(UL).shift(0.3*DOWN)
        sudef.next_to(smo, DOWN)
        sudef.shift((smo.get_left()-sudef.get_left())*RIGHT)
        sudef2.next_to(sudef,DOWN)
        sudef2.shift((smo.get_left()-sudef2.get_left())*RIGHT)
        sudef1.next_to(sudef,RIGHT)
        self.add_fixed_in_frame_mobjects(smo)
        self.add_fixed_in_frame_mobjects(sudef)
        self.add_fixed_in_frame_mobjects(sudef1)
        self.add_fixed_in_frame_mobjects(sudef2)
        self.wait()
        l1 = Line([-1,1.9,0], [-0.7,1.9,0], color = ORANGE)
        l11 = Line([1.25,1.9,0], [1.55,1.9,0], color = ORANGE)
        l2 = Line([-6.7,1.9,0], [-6.3,1.9,0], color = PURPLE_B)
        l22 = Line([-0.4,1.9,0], [0,1.9,0], color = PURPLE_B)
        l3 = Line([-4.1,1.4,0], [-3.1,1.4,0], color = GREEN)
        l4 = Line([1.8,1.9,0], [2.2,1.9,0], color = BLUE)
        self.add_fixed_in_frame_mobjects(l1, l2, l3, l4, l11, l22)
        self.wait()

class rottest(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        axes1 = ThreeDAxes(x_range = [-3,3,1], y_range = [-3,3,1], z_range = [0,3,1], x_length = 6, y_length = 6, z_length = 3 )
        self.play(Create(axes1))
        x = Dot([1,0,0], color = RED)
        y = Dot([0,1,0], color = BLUE)
        z = Dot([0,0,1], color = GREEN)
        self.add(x,y,z)
        self.wait()
        self.play(axes1.animate.rotate_about_origin(1, RIGHT))
        self.wait()
        self.begin_ambient_camera_rotation(rate = PI/2)
        self.wait(2)
        self.stop_ambient_camera_rotation()
        self.wait()


class tan(Scene):
    def construct(self):
        a = Dot()
        self.play(Create(a))


class texttest(ThreeDScene):
    def construct(self):
        cap = Tex(r"For ", r"$p$ ", r"$ \in \mathbb{R}^2$, ", r"$\varphi (p)$ ", r"$ \in \mathbb{S}^2$ is the point collinear with $p$ and the ", r"south pole.", font_size = 36).to_edge(UL).shift( 0.2*DOWN)
        cap.shift(( cap.get_center()[0])*LEFT)
        self.add_fixed_in_frame_mobjects(cap)
        self.wait()

