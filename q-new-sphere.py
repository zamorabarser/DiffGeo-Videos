from fnmatch import translate
from this import d
from tkinter import E
from manim import *
from numpy import sqrt




class tis(Scene):
    def construct(self):
        t1 = Text("Spherical curves", font_size=60)
        self.play(Write(t1))
        self.wait(3) 
        self.play(FadeOut(t1))
        self.wait()       



 

class spdef2(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES, distance = 3)
        axes = ThreeDAxes()
        param = Tex(r"$\gamma : [a,b] \to \mathbb{S}^2$", font_size=50)
        param.to_edge(UL).shift(RIGHT+0.3*DOWN)
        self.add_fixed_in_frame_mobjects(param)
        self.add(axes)
        self.wait()
        s22 = Surface(
            lambda u, v: np.array([
                1.5* np.cos(u) * np.cos(v),
                1.5* np.cos(u) * np.sin(v),
                1.5* np.sin(u)
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            fill_color=BLUE, resolution=(20, 20), fill_opacity=0.4
        )
        gam = ParametricFunction(lambda t : [1.5* np.cos((t-PI)/1.8)*np.cos(PI*np.sin(2*t)/4),1.5*np.sin((t-PI)/1.8)*np.cos(PI*np.sin(2*t)/4),1.5*np.sin(PI*np.sin(2*t)/4) ], t_range = [0 , TAU ] , color = YELLOW)
        self.play(Create(s22))
        self.wait()
        self.play(Create(gam), run_time=3)
        self.wait()
        self.begin_ambient_camera_rotation(rate=TAU/10)
        self.wait(10)
        self.stop_ambient_camera_rotation()
        self.wait()


class fre2d(Scene):
    def construct(self):
        ti = Tex(r"Hemisphere Lemma",font_size=40, color =BLUE ).to_edge(UL)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is closed and has length $<2 \pi$,", font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"then it is strictly contained in a hemisphere.", font_size=40).next_to(t1,DOWN)
        l = ti.get_left()[0]
        t1.shift((t1.get_left()[0]-l)*LEFT)
        t2.shift((t2.get_left()[0]-l)*LEFT)
        self.play(Write(ti), Write(t1), Write(t2))
        self.wait()


class fre(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES, distance = 8)
        ti = Tex(r"Hemisphere Lemma",font_size=40, color =BLUE ).to_edge(UL)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is closed and has length $<2 \pi$,", font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"then it is strictly contained in a hemisphere.", font_size=40).next_to(t1,DOWN)
        l = ti.get_left()[0]
        t1.shift((t1.get_left()[0]-l)*LEFT)
        t2.shift((t2.get_left()[0]-l)*LEFT)
        self.add_fixed_in_frame_mobjects(ti)
        self.add_fixed_in_frame_mobjects(t1)
        self.add_fixed_in_frame_mobjects(t2)
        self.wait()
        s2 = Surface(
            lambda u, v: np.array([
                2* np.cos(u) * np.cos(v),
                2* np.cos(u) * np.sin(v),
                2* np.sin(u)-3/2
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            fill_color=BLUE, resolution=(20, 20), fill_opacity=0.3
        )
        gam = ParametricFunction(lambda t : [ 2* np.cos(t) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10),2* np.sin(t) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10), 2*np.cos(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10) -3/2  ], t_range = [0 , TAU ] , color = YELLOW)
        self.play(Create(s2))
        self.play(Create(gam), run_time = 3)
        self.wait()
        self.begin_ambient_camera_rotation(rate=TAU/5)
        self.wait(5)
        self.stop_ambient_camera_rotation()
        self.wait()
        t3 = Tex(r"Assume $\gamma : [0,L] \to \mathbb{S}^2$ is parametrized by arc length, with $L < 2 \pi$.", font_size=40).next_to(ti,DOWN)
        t3.shift((t3.get_left()[0]-l)*LEFT+0.25* UP)
        self.play(FadeOut(ti), FadeOut(t1), FadeOut(t2))
        self.add_fixed_in_frame_mobjects(t3)
        self.wait()


        
class fre1(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES, distance = 8)
        ti = Tex(r"Hemisphere Lemma",font_size=40, color =BLUE ).to_edge(UL)
        l = ti.get_left()[0]
        t3 = Tex(r"Assume $\gamma : [0,L] \to \mathbb{S}^2$ is parametrized by arc length, with $L < 2 \pi$", font_size=40).next_to(ti,DOWN)
        t3.shift((t3.get_left()[0]-l)*LEFT+0.25* UP)
        pq1 = Tex(r"$p = \gamma (0)$", font_size=40).next_to(t3,DOWN)
        pq1.shift((pq1.get_left()[0]-l -0.2)*LEFT+0.1*DOWN)
        pq1b = SurroundingRectangle(pq1, buff = .1, color = RED)
        pq2 = Tex(r"$q = \gamma (L/2)$", font_size=40).next_to(pq1,RIGHT).shift(0.2*RIGHT)
        pq2b = SurroundingRectangle(pq2, buff = .1, color = RED)
        g1 = Tex(r"$\gamma_1=$ second half of $\gamma$", font_size=40).next_to(pq2,RIGHT).shift(0.2*RIGHT)
        g1b = SurroundingRectangle(g1, buff = .1, color = YELLOW)
        x1 = Tex(r"$x = $ midpoint between $p$ and $q$", font_size=40).next_to(pq1,DOWN)
        x1.shift((x1.get_left()[0]-l -0.2)*LEFT + 0.1*DOWN)
        x1b = SurroundingRectangle(x1, buff = .1, color = ORANGE)
        g2 = Tex(r"$\gamma_2=$ rotation of $\gamma_1$", font_size=40).next_to(x1,RIGHT).shift(0.2*RIGHT)
        g2b = SurroundingRectangle(g2, buff = .1, color = PINK)
        self.add_fixed_in_frame_mobjects(t3, pq1, pq2, pq1b, pq2b)
        s2 = Surface(
            lambda u, v: np.array([
                2* np.cos(u) * np.cos(v),
                2* np.cos(u) * np.sin(v),
                2* np.sin(u)-3/2
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            fill_color=BLUE, resolution=(20, 20), fill_opacity=0.3
        )
        gam = ParametricFunction(lambda t : [ 2* np.cos(t) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10),2* np.sin(t) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10), 2*np.cos(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10) -3/2  ], t_range = [0 , TAU ] , color = YELLOW)
        self.add(s2, gam)
        p = Sphere([2*np.sin(PI/6), 0, 2*np.cos(PI/6)- 3/2], color = RED , radius = 0.1)
        p.set_color(RED)
        q = Sphere([-2*np.sin(PI/6), 0, 2*np.cos(PI/6)-3/2], color = RED , radius = 0.1 )
        q.set_color(RED)
        x = Sphere([0, 0, 1/2], color = ORANGE , radius = 0.1 )
        x.set_color(ORANGE)
        self.play(Create(p), Create(q))
        self.wait()
        gam1 = ParametricFunction(lambda t : [ 2* np.cos(t) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10),2* np.sin(t) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10), 2*np.cos(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10) -3/2  ], t_range = [PI , TAU ] , color = YELLOW)
        gam2 = gam1.copy()
        gam2.set_color(PINK)
        self.add(gam1)
        self.add_fixed_in_frame_mobjects(g1, g1b)
        self.play(FadeOut(gam))
        self.wait()
        self.add_fixed_in_frame_mobjects(x1, x1b)
        self.play(Create(x))
        self.wait()
        self.add(gam2)
        self.add_fixed_in_frame_mobjects(g2,g2b)
        s = ValueTracker(0)
        gam2.add_updater(  lambda x: x.become( ParametricFunction(lambda t : [ 2* np.cos(t+s.get_value()) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10),2* np.sin(t+s.get_value()) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10), 2*np.cos(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10) -3/2  ], t_range = [PI,TAU ] , color = PINK)))
        self.play(s.animate.set_value(PI), run_time = 3)
        self.wait()
        

  
class fre2(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES, distance = 8)
        ti = Tex(r"Hemisphere Lemma",font_size=40, color =BLUE ).to_edge(UL)
        l = ti.get_left()[0]
        t3 = Tex(r"Assume $\gamma : [0,L] \to \mathbb{S}^2$ is parametrized by arc length, with $L < 2 \pi$.", font_size=40).next_to(ti,DOWN)
        t3.shift((t3.get_left()[0]-l)*LEFT+0.25* UP)
        pq1 = Tex(r"$p = \gamma (0)$", font_size=40).next_to(t3,DOWN)
        pq1.shift((pq1.get_left()[0]-l -0.2)*LEFT+ 0.1*DOWN)
        pq1b = SurroundingRectangle(pq1, buff = .1, color = RED)
        pq2 = Tex(r"$q = \gamma (L/2)$", font_size=40).next_to(pq1,RIGHT).shift(0.2*RIGHT)
        pq2b = SurroundingRectangle(pq2, buff = .1, color = RED)
        g1 = Tex(r"$\gamma_1=$ second half of $\gamma$", font_size=40).next_to(pq2,RIGHT).shift(0.2*RIGHT)
        g1b = SurroundingRectangle(g1, buff = .1, color = YELLOW)
        x1 = Tex(r"$x = $ midpoint between $p$ and $q$", font_size=40).next_to(pq1,DOWN)
        x1.shift((x1.get_left()[0]-l -0.2)*LEFT+ 0.1*DOWN)
        x1b = SurroundingRectangle(x1, buff = .1, color = ORANGE)
        g2 = Tex(r"$\gamma_2=$ rotation of $\gamma_1$", font_size=40).next_to(x1,RIGHT).shift(0.2*RIGHT)
        g2b = SurroundingRectangle(g2, buff = .1, color = PINK)
        y1 = Tex(r"$y \in \gamma _1 $, $\langle y, x \rangle = 0$", font_size=40).next_to(x1,DOWN)
        y1.shift((y1.get_left()[0]-l -0.2)*LEFT+ 0.1*DOWN)
        y1b = SurroundingRectangle(y1, buff = .1, color = GREEN)
        z1 = Tex(r"$z \in \gamma _2 $, $\langle z, x \rangle = 0$", font_size=40).next_to(y1,DOWN)
        z1.shift((y1.get_left()[0]-l -0.2)*LEFT+ 0.1*DOWN)
        z1b = SurroundingRectangle(z1, buff = .1, color = GREEN)
        self.add_fixed_in_frame_mobjects(t3, pq1, pq2, pq1b, pq2b, g1, g1b, g2, g2b, x1, x1b)
        s2 = Surface(
            lambda u, v: np.array([
                2* np.cos(u) * np.cos(v),
                2* np.cos(u) * np.sin(v),
                2* np.sin(u)-3/2
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            fill_color=BLUE, resolution=(20, 20), fill_opacity=0.3
        )
        gam = ParametricFunction(lambda t : [ 2* np.cos(t) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10),2* np.sin(t) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10), 2*np.cos(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10) -3/2  ], t_range = [PI , TAU ] , color = YELLOW)
        gam2 = ParametricFunction(lambda t : [ 2* np.cos(t+PI) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10),2* np.sin(t+PI) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10), 2*np.cos(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)/10) -3/2  ], t_range = [PI , TAU ] , color = PINK)
        p = Sphere([2*np.sin(PI/6), 0, 2*np.cos(PI/6)- 3/2], color = RED , radius = 0.1)
        p.set_color(RED)
        q = Sphere([-2*np.sin(PI/6), 0, 2*np.cos(PI/6)-3/2], color = RED , radius = 0.1 )
        q.set_color(RED)
        x = Sphere([0, 0, 1/2], color = ORANGE , radius = 0.1 )
        x.set_color(ORANGE)
        self.add(s2, gam, gam2, p, q, x)
        self.wait()
        eq = ParametricFunction(lambda t : [ 2* np.cos(t) , 2*np.sin(t), -3/2  ], t_range = [0 , TAU ] , color = BLACK)
        self.play(Create(eq))
        self.wait()
        y = Sphere([ 2* np.cos(4.41) * np.sin(PI/3 + PI * np.cos(8.82 - PI)/6 - np.sin(4.41)/10), 2* np.sin(4.41) * np.sin(PI/3 + PI * np.cos(8.82 - PI)/6 - np.sin(4.41)/10), 2*np.cos(PI/3 + PI * np.cos(8.82 - PI)/6 - np.sin(4.41)/10) -3/2 ], radius = 0.1)
        y.set_color(GREEN)
        z = Sphere([ - 2* np.cos(4.41) * np.sin(PI/3 + PI * np.cos(8.82 - PI)/6 - np.sin(4.41)/10), - 2* np.sin(4.41) * np.sin(PI/3 + PI * np.cos(8.82 - PI)/6 - np.sin(4.41)/10), 2*np.cos(PI/3 + PI * np.cos(8.82 - PI)/6 - np.sin(4.41)/10) -3/2 ], radius = 0.1)
        z.set_color(GREEN)
        self.add_fixed_in_frame_mobjects(y1,y1b,z1,z1b)
        self.play(Create(y), Create(z))
        self.wait()
        self.begin_ambient_camera_rotation(rate=TAU/5)
        self.wait(5)
        self.stop_ambient_camera_rotation()
        
        






class frec(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES, distance = 8)
        ti = Tex(r"Hemisphere Lemma",font_size=40, color =BLUE ).to_edge(UL)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is closed and has length $<2 \pi$,", font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"then it is strictly contained in a hemisphere.", font_size=40).next_to(t1,DOWN)
        tv = Tex(r"There is $v \in \mathbb{S}^2$ with $ \langle v, \gamma (t) \rangle > 0$ for all $t$.", font_size=40).next_to(t2,DOWN)
        vul = Line([-5.07,1.1,0],[-4.2,1.1,0],color=ORANGE)
        l = ti.get_left()[0]
        t1.shift((t1.get_left()[0]-l)*LEFT)
        t2.shift((t2.get_left()[0]-l)*LEFT)
        tv.shift((tv.get_left()[0]-l)*LEFT)
        self.add_fixed_in_frame_mobjects(ti)
        self.add_fixed_in_frame_mobjects(t1)
        self.add_fixed_in_frame_mobjects(t2)
        self.wait()
        s2 = Surface(
            lambda u, v: np.array([
                2* np.cos(u) * np.cos(v),
                2* np.cos(u) * np.sin(v),
                2* np.sin(u)-3/2
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            fill_color=BLUE, resolution=(20, 20), fill_opacity=0.3
        )
        gam = ParametricFunction(lambda t : [ 2* np.cos(t) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)**2/10),2* np.sin(t) * np.sin(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)**2/10), 2*np.cos(PI/3 + PI * np.cos(2*t - PI)/6 - np.sin(t)**2/10) -3/2  ], t_range = [0 , TAU ] , color = YELLOW)
        v = Arrow([0,0,-3/2],[0,0,0.8], color=ORANGE)
        c = ParametricFunction(lambda t : [ 2* np.cos(t) , 2*np.sin(t),-3/2  ], t_range = [0 , TAU ] , color = BLACK)
        self.play(Create(s2))
        self.play(Create(gam), run_time = 3)
        self.wait()
        self.add_fixed_in_frame_mobjects(tv, vul)
        self.play(Create(c), Create(v))
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/5)
        self.wait(5)
        self.stop_ambient_camera_rotation()


class frec1(Scene):
    def construct(self):
        ti = Tex(r"Hemisphere Lemma",font_size=40, color =BLUE ).to_edge(UL)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is closed and has length $<2 \pi$,", font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"then it is strictly contained in a hemisphere.", font_size=40).next_to(t1,DOWN)
        tv = Tex(r"There is $v \in \mathbb{S}^2$ with $ \langle v, \gamma (t) \rangle > 0$ for all $t$.", font_size=40).next_to(t2,DOWN)
        l = ti.get_left()[0]
        t1.shift((t1.get_left()[0]-l)*LEFT)
        t2.shift((t2.get_left()[0]-l)*LEFT)
        tv.shift((tv.get_left()[0]-l)*LEFT)
        self.add(ti,t1,t2,tv)
        co = Tex(r"Fenchel Theorem",font_size=40, color =BLUE ).next_to(tv,DOWN)
        co.shift((co.get_left()[0]-l)*LEFT+0.3*DOWN)
        c1 = Tex(r"If $\gamma:[a,b ] \to \mathbb{R}^3$ is closed and piecewise smooth regular, then $\Phi (\gamma ) \geq 2 \pi$.", font_size=40).next_to(co,DOWN)
        c1.shift((c1.get_left()[0]-l)*LEFT)
        self.play(Create(co), Create(c1))
        self.wait()        
        p = Tex(r"Proof:", font_size=40, color=BLUE).next_to(c1,DOWN)
        p.shift((p.get_left()[0]-l)*LEFT)
        p1 = Tex(r"If not, $\Phi (\gamma) = $ length$(T) < 2 \pi$",r", so there is", font_size=40).next_to(p,DOWN)
        p2 = Tex(r"$v \in \mathbb{S}^2$ with $\langle v , T (t) \rangle > 0 $ for all $t$.", font_size=40).next_to(p1,DOWN)
        p3 = Tex(r"$ 0 = \langle v , \gamma ( b) - \gamma ( a ) \rangle = \int_a^b \langle v, \gamma ^{\prime} (t) \rangle dt > 0 . $", font_size=40).shift(3*DOWN)
        p1.shift((p1.get_left()[0]-l)*LEFT)
        p2.shift((p2.get_left()[0]-l)*LEFT)
        self.play(Create(p), Write(p1[0]))
        self.wait()
        self.play(Write(p1[1]), Write(p2))
        self.wait()
        self.play(Create(p3))




class ch(Scene):
    def construct(self):
        ti = Tex(r"Chord Lemma",font_size=40, color =BLUE ).to_edge(UL).shift(0.5*DOWN+RIGHT)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{R}^3$ is piecewise smooth regular,", font_size=40).next_to(ti,DOWN)
        t2 = Tex(r" $\mathrm{W} = \gamma (b) - \gamma (a) $, $\alpha = \measuredangle ( \mathrm{W} , \gamma ^{\prime} (a )  ) $, $\beta = \measuredangle (\mathrm{W} , \gamma ^{\prime}(b))$", font_size=40).next_to(t1,DOWN)
        t3 = Tex(r"Then $\Phi (\gamma) \geq \alpha + \beta$.", font_size=40).next_to(t2,DOWN)
        l = ti.get_left()[0]
        t1.shift((t1.get_left()[0]-l)*LEFT)
        t2.shift((t2.get_center()[0])*LEFT)
        t3.shift((t3.get_left()[0]-l)*LEFT)
        self.play(Write(ti),Write(t1),Write(t2),Write(t3))
        p = Tex(r"Proof:", font_size=40, color=BLUE).next_to(t3,DOWN)
        p1 = Tex(r"By Fenchel Theorem, $\Phi (\gamma) + (\pi - \alpha) + (\pi - \beta ) \geq 2 \pi$.", font_size=40).next_to(p,DOWN)
        p.shift((p.get_left()[0]-l)*LEFT)
        p1.shift((p1.get_left()[0]-l)*LEFT)
        ga = Tex(r"$\gamma$", font_size=40).shift(0.6*DOWN)
        gab = SurroundingRectangle(ga, buff = .1, color = YELLOW)
        g = ParametricFunction(lambda t : [  4*(t-PI/2) , 5*np.sin(t) + np.sin( 6*(t-PI/3) )/12 -6,0 ], t_range = [PI/3 , 2*PI/3 ] , color = YELLOW)
        w = Line([-2*PI/3 , 5* sqrt(3)/2 - 6 , 0],[2*PI/3 , 5* sqrt(3)/2 - 6 ,0], color = ORANGE)
        wr = Line([2*PI/3 , 5* sqrt(3)/2 - 6 ,0],[-2*PI/3 , 5* sqrt(3)/2 - 6 , 0])
        wl = Tex(r"$\mathrm{W}$", font_size=40).shift(2.1*DOWN)
        wlb = SurroundingRectangle(wl, buff = .1, color = ORANGE)
        gpa = Line([-2*PI/3 , 5* sqrt(3)/2 - 6 , 0],[-2*PI/3 + 4 , 5* sqrt(3)/2 - 3.2 , 0])
        gpb = Line([2*PI/3 , 5* sqrt(3)/2 - 6 ,0], [2*PI/3 - 4 , 5* sqrt(3)/2 - 4.3 ,0])
        al = Tex(r"$\alpha$", font_size=40).shift(1.5*LEFT + 2.1*DOWN)
        alb = SurroundingRectangle(al, buff = .1, color = TEAL)
        bl = Tex(r"$\beta$", font_size=40).shift(1.5*RIGHT + 2.1*DOWN)
        blb = SurroundingRectangle(bl, buff = .1, color = GREEN)
        alp = Angle(w,gpa, radius=0.7, color = TEAL)
        bet = Angle(gpb,wr, radius=0.7, color = GREEN)
        all = Group(g,ga,gab,w,wl,wlb,al,alb,bl,blb,alp,bet)
        all.shift(0.8*DOWN)
        self.play(Write(p))
        self.play(Create(g))
        self.play(Create(ga), Create(gab))
        self.play(Create(w))
        self.play(Create(wl), Create(wlb))
        self.play(Create(al), Create(alb), Create(bl), Create(blb), Create(alp), Create(bet))
        self.wait()
        self.play(Write(p1))
        self.wait()
        self.play(FadeOut(ga), FadeOut(gab), FadeOut(wl), FadeOut(wlb), FadeOut(al), FadeOut(alb), FadeOut(bl), FadeOut(blb))
        self.wait()
        
        


class ch1(Scene):
    def construct(self):
        ti = Tex(r"Chord Lemma",font_size=40, color =BLUE ).to_edge(UL).shift(0.5*DOWN+RIGHT)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{R}^3$ is piecewise smooth regular,", font_size=40).next_to(ti,DOWN)
        t2 = Tex(r" $\mathrm{W} = \gamma (b) - \gamma (a) $, $\alpha = \measuredangle ( \mathrm{W} , \gamma ^{\prime} (a )  ) $, $\beta = \measuredangle (\mathrm{W} , \gamma ^{\prime}(b))$", font_size=40).next_to(t1,DOWN)
        t3 = Tex(r"Then $\Phi (\gamma) \geq \alpha + \beta$.", font_size=40).next_to(t2,DOWN)
        l = ti.get_left()[0]
        t1.shift((t1.get_left()[0]-l)*LEFT)
        t2.shift((t2.get_center()[0])*LEFT)
        t3.shift((t3.get_left()[0]-l)*LEFT)
        self.play(Write(ti),Write(t1),Write(t2),Write(t3))
        p = Tex(r"Proof:", font_size=40, color=BLUE).next_to(t3,DOWN)
        p1 = Tex(r"By Fenchel Theorem, $\Phi (\gamma) + (\pi - \alpha) + (\pi - \beta ) \geq 2 \pi$.", font_size=40).next_to(p,DOWN)
        p.shift((p.get_left()[0]-l)*LEFT)
        p1.shift((p1.get_left()[0]-l)*LEFT)
        g = ParametricFunction(lambda t : [  4*(t-PI/2) , 5*np.sin(t) + np.sin( 6*t )/12 -6,0 ], t_range = [PI/3 , 2*PI/3 ] , color = YELLOW)
        w = Line([-2*PI/3 , 5* sqrt(3)/2 - 6 , 0],[2*PI/3 , 5* sqrt(3)/2 - 6 ,0], color = ORANGE)
        wr = Line([2*PI/3 , 5* sqrt(3)/2 - 6 ,0],[-2*PI/3 , 5* sqrt(3)/2 - 6 , 0])
        gpa = Line([-2*PI/3 , 5* sqrt(3)/2 - 6 , 0],[-2*PI/3 + 4 , 5* sqrt(3)/2 - 3.2 , 0])
        gpb = Line([2*PI/3 , 5* sqrt(3)/2 - 6 ,0], [2*PI/3 - 4 , 5* sqrt(3)/2 - 4.3 ,0])
        alp = Angle(w,gpa, radius=0.7, color = TEAL)
        bet = Angle(gpb,wr, radius=0.7, color = GREEN)
        all = Group(g,w,alp,bet)
        all.shift(0.8*DOWN)
        self.add(all,p,p1)
        self.wait()
        tt = Arrow([-2*PI/3 , 5* sqrt(3)/2 - 6.8 , 0], [-2*PI/3 +4/3 , 5* sqrt(3)/2 - 6.8 + 1  , 0], color = PINK)
        self.play(Create(tt))
        self.wait()
        s = ValueTracker(PI/3)
        tt.add_updater( lambda x: x.become( Arrow([ 4* s.get_value() - 2 * PI , 5*np.sin(s.get_value()) + np.sin(6*s.get_value())/12 -6.8, 0  ], [4* s.get_value() - 2 * PI + 4/3 ,  5*np.sin(s.get_value()) + np.sin(6*s.get_value())/12 -6.8 + 5*np.cos(s.get_value())/3 + np.cos(6*s.get_value())/6, 0 ], color = PINK) ))
        self.play(s.animate.set_value(2*PI/3), run_time = 2)
        self.wait()
        ttc = Arrow([8*PI/3 - 2*PI , 5 * sqrt(3)/2 - 6.8 , 0], [8*PI/3 - 2*PI + 4/3 , 5 * sqrt(3)/2 - 6.8 -2/3 , 0], color = PINK)
        self.remove(tt)
        self.add(ttc)
        self.play(Rotate(ttc, about_point = [8*PI/3 - 2*PI , 5 * sqrt(3)/2 - 6.8 , 0], angle =  np.arctan(1/2) - PI ) )
        self.wait()
        self.play(ttc.animate.shift(4*PI/3 *LEFT))
        self.wait()
        self.play(Rotate(ttc, about_point = [-2*PI/3 , 5* sqrt(3)/2 - 6.8 , 0], angle =  np.arctan(3/4) - PI ) )
        self.wait()
        

class dc(Scene):
    def construct(self):
        ti = Tex(r"Discrete Chord Lemma",font_size=40, color =BLUE ).to_edge(UL).shift(0.5*RIGHT)
        t1 = Tex(r"If $a,b,c,d,x \in \mathbb{R}^3$ are distinct, then", font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"$\Phi ( abcd ) \leq \Phi (abxcd)$.", font_size=40).next_to(t1,DOWN)
        t22 = Tex(r"$\Phi ( abcd ) \leq \Phi (abxcd)$", font_size=40).shift(2*UP)
        l = ti.get_left()[0]
        t1.shift((t1.get_left()[0]-l)*LEFT)
        t2.shift((t2.get_center()[0])*LEFT)
        self.play(Write(ti),Write(t1),Write(t2))
        p = Tex(r"Proof:", font_size=40, color=BLUE).next_to(t2,DOWN)
        p1 = Tex(r"Let $y$ be between $a$ and $b$, and $z$ between $c$ and $d$.",  font_size=40).next_to(p,DOWN)
        p11 = Tex(r"Let $y$ be between $a$ and $b$, and $z$ between $c$ and $d$.", r" Then", font_size=40).to_edge(UL).shift(0.5*RIGHT+0.5*DOWN)
        p2 = Tex(r"$\Phi (abxcd) = \Phi (ybxcz)$", r" $ \geq \measuredangle zyb + \measuredangle yzc $",r" $ =  2 \pi - \measuredangle ayz + \measuredangle yzd $.", font_size=40).next_to(p11,DOWN)
        p3 = Tex(r"As $y \to b$, and $z \to c$, we have", font_size=40).next_to(p2,DOWN)
        p4 = Tex(r"$\measuredangle ayz + \measuredangle yzd \to \measuredangle abc + \measuredangle bcd$", r" $ = 2 \pi  -  \Phi (abcd) $", font_size=40).next_to(p3,DOWN)
        p.shift((p.get_left()[0]-l)*LEFT)
        p1.shift((p1.get_left()[0]-l)*LEFT)
        p2.shift((p2.get_center()[0])*LEFT)
        p3.shift((p3.get_left()[0]-l)*LEFT)
        p4.shift((p4.get_center()[0])*LEFT)
        a = Dot([-2.5,-3.5,0], color = RED)
        b = Dot([-2,-1.5,0], color = YELLOW)
        c = Dot([2,-1.5,0], color = BLUE)
        d = Dot([2.5,-3.5,0], color = PINK)
        x = Dot([0,-0.5,0], color = GREEN)
        y = Dot([-2.125, -2,0 ], color = ORANGE)
        z = Dot([2.125, -2,0 ], color = PURPLE)
        dal = DashedLine(y,z)
        al = Tex(r"$a$", font_size=40).next_to(a,LEFT)
        bl = Tex(r"$b$", font_size=40).next_to(b,LEFT)
        cl = Tex(r"$c$", font_size=40).next_to(c,RIGHT)
        dl = Tex(r"$d$", font_size=40).next_to(d,RIGHT)
        xl = Tex(r"$x$", font_size=40).next_to(x,UP).shift(0.5*RIGHT)
        yl = Tex(r"$y$", font_size=40).next_to(y,RIGHT).shift(0.5*DOWN)
        zl = Tex(r"$z$", font_size=40).next_to(z,LEFT).shift(0.5*DOWN)
        alb = SurroundingRectangle(al, buff = .1, color = RED)
        blb = SurroundingRectangle(bl, buff = .1, color = YELLOW)
        clb = SurroundingRectangle(cl, buff = .1, color = BLUE)
        dlb = SurroundingRectangle(dl, buff = .1, color = PINK)
        xlb = SurroundingRectangle(xl, buff = .1, color = GREEN)
        ylb = SurroundingRectangle(yl, buff = .1, color = ORANGE)
        zlb = SurroundingRectangle(zl, buff = .1, color = PURPLE)
        l1 = Line([-2.5,-3.5,0], [-2,-1.5,0])
        l2 = Line([-2,-1.5,0], [0,-0.5,0])
        l3 = Line([0,-0.5,0], [2,-1.5,0])
        l4 = Line([2,-1.5,0],[2.5,-3.5,0])
        self.play(Create(l1),Create(l2), Create(l3), Create(l4) )
        self.play(Create(a), Create(b), Create(c), Create(d), Create(x), Write(al), Write(bl), Write(cl), Write(dl), Write(xl), Create(alb), Create(blb), Create(clb), Create(dlb), Create(xlb))
        self.wait()
        self.play(Write(p))
        self.play(Write(p1[0]))
        self.play(Create(y), Create(z), Create(yl), Create(zl), Create(ylb), Create(zlb))
        self.wait()
        self.play(p1.animate.shift((p1.get_center()[1] - p11.get_center()[1])*DOWN), FadeOut(ti), FadeOut(t1), FadeOut(t2), FadeOut(p))
        self.remove(p1)
        self.add(p11[0])
        self.wait()
        self.play(Write(p11[1]), Write(p2[0]))
        self.wait()
        self.play(Write(p2[1]), Create(dal))
        self.wait()
        self.play(Write(p2[2]))
        self.wait()
        self.play(Write(p3), Write(p4[0]))
        self.wait()
        s = ValueTracker(3/4)
        y.add_updater(  lambda x: x.become( Dot((1- s.get_value())*a.get_center( ) + s.get_value()* b.get_center() , color = ORANGE )    ))
        z.add_updater(  lambda x: x.become( Dot((1- s.get_value())*d.get_center( ) + s.get_value()* c.get_center() , color = PURPLE )    ))
        yl.add_updater(  lambda x: x.become(  Tex(r"$y$", font_size=40).next_to(y,RIGHT).shift(0.5*DOWN) ))
        zl.add_updater(  lambda x: x.become( Tex(r"$z$", font_size=40).next_to(z,LEFT).shift(0.5*DOWN)))
        ylb.add_updater(  lambda x: x.become(  SurroundingRectangle(yl, buff = .1, color = ORANGE ) ))
        zlb.add_updater(  lambda x: x.become(   SurroundingRectangle(zl, buff = .1, color = PURPLE ) ))
        dal.add_updater(  lambda x: x.become(DashedLine(y,z)))
        self.play(s.animate.set_value(1), run_time = 3)
        self.wait()
        self.play(Write(p4[1]))
        self.wait()
        final = SurroundingRectangle(t22, buff = .3, color = TEAL )
        self.play(FadeOut(p11), FadeOut(p1), FadeOut(p2), FadeOut(p3), FadeOut(p4))
        self.play( Write(t22), Create(final))
        self.wait()


class tc(Scene):
    def construct(self):
        ti = Tex(r"Theorem",font_size=40, color =BLUE ).to_edge(UL).shift(0.5*RIGHT)
        l =  ti.get_left()[0]
        t1 = Tex(r"For a piecewise smooth regular curve $\gamma : [a,b] \to \mathbb{R}^3$,", font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"$\Phi ( \gamma  )  = \sup \Phi (  p_0 , \ldots, p_n   )$.", font_size=40).next_to(t1,DOWN).shift(0.5*DOWN)
        t1.shift((t1.get_left()[0]-l)*LEFT)
        t2.shift((t2.get_center()[0])*LEFT)
        f0 = SurroundingRectangle(t2, buff = .15, color = PINK )
        self.play(Write(ti), Write(t1), Write(t2), Create(f0))
        pi = Tex(r"Proof",font_size=40, color =BLUE ).next_to(t2, DOWN).shift(0.2*DOWN)
        p1 = Tex(r"By the chord lemma, ", font_size=40).next_to(pi,DOWN)
        p2 = Tex(r"$\Phi ( \gamma  )  \geq $",r" $ \Phi (  p_0 , \ldots, p_n   )$.", font_size=40).next_to(p1,DOWN).shift(0.5*DOWN)
        p22 = Tex(r"$\Phi ( \gamma  )  \geq $",r" $ \sup $", r" $\Phi (  p_0 , \ldots, p_n   )$.", font_size=40).next_to(p1,DOWN).shift(0.5*DOWN)
        pi.shift((pi.get_left()[0]-l)*LEFT)
        p1.shift((p1.get_left()[0]-l)*LEFT)
        p2.shift((p2.get_center()[0])*LEFT)
        p22.shift((p22.get_center()[0])*LEFT)
        self.play(Write(pi), Write(p1), Write(p2))
        self.wait()
        self.play(p2[0].animate.become(p22[0]), Write(p22[1]), p2[1].animate.become(p22[2]))
        self.wait()
        p3 = Tex(r"If $\gamma_n$ are piecewise linear inscribed in $\gamma$ and $\gamma_n \to \gamma$, then:", font_size=40).next_to(p2,DOWN)
        p4 = Tex(r"$\lim \Phi ( \gamma _n ) \geq \Phi (\gamma )$.", font_size=40).next_to(p3,DOWN)
        p3.shift((p3.get_left()[0]-l)*LEFT)
        p4.shift((p4.get_center()[0])*LEFT)
        self.play(Write(p3), Write(p4))
        self.wait()
        self.play(FadeOut(pi), FadeOut(p1), FadeOut(p22), FadeOut(p2), FadeOut(p4), FadeOut(p3))        
        ci = Tex(r"Theorem",font_size=40, color =BLUE ).to_edge(UL).shift(0.5*RIGHT + 3*DOWN)
        c1 = Tex(r"If $\gamma : [a,b] \to \mathbb{R}^3$ is closed and $\vert \gamma (t) \vert \leq 1$ for all $t$, then", font_size=40).next_to(ci,DOWN)
        c2 = Tex(r"$\Phi ( \gamma  )  \geq $ length$(\gamma)$.", font_size=40).next_to(c1,DOWN).shift(0.5*DOWN)
        c1.shift((c1.get_left()[0]-l)*LEFT)
        c2.shift((c2.get_center()[0])*LEFT)
        final = SurroundingRectangle(c2, buff = .15, color = PINK )
        self.play(Write(ci), Write(c1), Write(c2), Create(final))
        self.wait()





class pr(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=30 * DEGREES, distance = 8)
        ti = Tex(r"For $x,v \in \mathbb{S}^2$, $x_v^{\ast} : = x_{v^{\perp}} / \vert x_{v^{\perp}} \vert  $",font_size=40 ).to_edge(UL)
        lx = Line([-5.55,3,0], [-5.85,3,0], color = ORANGE)
        lv = Line([-5.45,3,0], [-5.15,3,0], color = YELLOW)
        lp = Line([-4.05,3,0], [-3.6,3,0], color = TEAL)
        self.add_fixed_in_frame_mobjects(ti, lx, lv, lp)
        s2 = Surface(
            lambda u, v: np.array([
                2* np.cos(u) * np.cos(v),
                2* np.cos(u) * np.sin(v),
                2* np.sin(u)
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            fill_color=BLUE, resolution=(20, 20), fill_opacity=0.3, checkerboard_colors=[MAROON_B, MAROON_B]
        )
        g = ParametricFunction(lambda t : [ 2*np.cos(t), 2*np.sin(t),0 ], t_range = [ 0 , TAU ] , color = BLACK)
        self.play(Create(s2))
        x = Sphere([1, 0, sqrt(3)], radius = 0.1)
        v = Sphere([0,0,2], radius = 0.1)
        x.set_color(ORANGE)
        v.set_color(YELLOW)
        p = x.copy()        
        self.play(Create(x), Create(v), Create(p), Create(g))
        p.set_color(TEAL)
        s = ValueTracker(PI/6)
        p.add_updater(  lambda x: x.become(  Sphere([2*np.sin(s.get_value()),0,2*np.cos(s.get_value())], radius = 0.1).set_color(TEAL) ))
        self.play(s.animate.set_value(PI/2), run_time = 2)
        self.play(FadeOut(p), FadeOut(x))
        x = Sphere([0, sqrt(2), -sqrt(2)], radius = 0.1)
        x.set_color(ORANGE)
        p1 = x.copy().set_color(TEAL)
        self.play(Create(x))
        self.add(p1)
        s = ValueTracker(3*PI/4)
        p1.add_updater(  lambda x: x.become(  Sphere([0,2*np.sin(s.get_value()),2*np.cos(s.get_value())], radius = 0.1).set_color(TEAL) ))
        self.play(s.animate.set_value(PI/2), run_time = 2)
        self.play(FadeOut(p1), FadeOut(x))
        x = Sphere([sqrt(3/2), sqrt(3/2), 1], radius = 0.1)
        x.set_color(ORANGE)
        p2 = x.copy().set_color(TEAL)
        self.play(Create(x))
        self.add(p2)
        s = ValueTracker(PI/3)
        p2.add_updater(  lambda x: x.become(  Sphere([ sqrt(2)*np.sin(s.get_value()),sqrt(2)*np.sin(s.get_value()),2*np.cos(s.get_value())], radius = 0.1).set_color(TEAL) ))
        self.play(s.animate.set_value(PI/2), run_time = 2)
        self.wait()



class cr(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=-90 * DEGREES, distance = 8)
        ti = Tex(r"Theorem",font_size=40, color=BLUE).to_edge(UL)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is piecewise smooth regular, then",font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"length$(\gamma ) = \frac{1}{4 \pi} \int _{ \mathbb{S}^2} $ length$(\gamma_v^{\ast})dv $",font_size=40).next_to(t1,DOWN)
        l = ti.get_left()
        t1.shift((t1.get_left()-l)*LEFT)
        t2.shift((t2.get_center())*LEFT)
        self.add_fixed_in_frame_mobjects(ti, t1,t2)
        self.wait()
        s2 = Surface(
            lambda u, v: np.array([
                2* np.cos(u) * np.cos(v),
                2* np.cos(u) * np.sin(v),
                2* np.sin(u)-3/2
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            fill_color=BLUE, resolution=(20, 20), fill_opacity=0.3, checkerboard_colors=[MAROON_B, MAROON_B]
        )
        g = ParametricFunction(lambda t : [ 2*np.cos(t)*np.sin( PI/2 - np.sin(5*PI/6 - t)), 2*np.sin(t)*np.sin( PI/2 - np.sin(5*PI/6 - t)),2*np.cos( PI/2 - np.sin(5*PI/6 - t))-3/2 ], t_range = [ -PI/2 , PI ] , color = ORANGE)
        v1 = Sphere([0,0,1/2], radius = 0.1)
        v1.set_color(YELLOW)
        e1 = ParametricFunction(lambda t : [ 2*np.cos(t), 2*np.sin(t),-3/2 ], t_range = [ 0 , TAU ] , color = BLACK)
        self.play(Create(s2))
        p1 = g.copy()        
        self.play( Create(g))
        self.play(Create(v1), Create(e1))
        gl = Tex(r"$\gamma$",font_size=40).shift(4*LEFT)
        vl = Tex(r"$v$",font_size=40).next_to(gl, DOWN).shift(0.3*DOWN)
        gvl = Tex(r"$\gamma_v^{\ast}$",font_size=40).next_to(vl, DOWN).shift(0.3*DOWN)
        glb = SurroundingRectangle(gl, buff = .15, color = ORANGE )
        vlb = SurroundingRectangle(vl, buff = .15, color = YELLOW )
        gvlb = SurroundingRectangle(gvl, buff = .15, color = TEAL )
        self.add_fixed_in_frame_mobjects(gl, glb, vl, vlb)
        self.begin_ambient_camera_rotation(rate=PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.add(p1)
        p1.set_color(TEAL)
        s = ValueTracker(1)
        p1.add_updater(  lambda x: x.become( ParametricFunction(lambda t : [ 2*np.cos(t)*np.sin( PI/2 - s.get_value()*( np.sin(5*PI/6 - t))), 2*np.sin(t)*np.sin( PI/2 - s.get_value()*( np.sin(5*PI/6 - t))),2*np.cos( PI/2 -s.get_value()*( np.sin(5*PI/6 - t)) ) -3/2], t_range = [ -PI/2 , PI ] , color = TEAL) ))
        self.play(s.animate.set_value(0), run_time = 2)
        self.add_fixed_in_frame_mobjects(gvl, gvlb)
        self.wait()
        



class cr1(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=30 * DEGREES, distance = 8)
        ti = Tex(r"Theorem",font_size=40, color=BLUE).to_edge(UL)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is piecewise smooth regular, then",font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"length$(\gamma ) = \frac{1}{4 \pi} \int _{ \mathbb{S}^2} $ length$(\gamma_v^{\ast})dv $",font_size=40).next_to(t1,DOWN)
        l = ti.get_left()
        t1.shift((t1.get_left()-l)*LEFT)
        t2.shift((t2.get_center())*LEFT)
        gl = Tex(r"$\gamma$",font_size=40).shift(4*LEFT)
        vl = Tex(r"$v$",font_size=40).next_to(gl, DOWN).shift(0.3*DOWN)
        gvl = Tex(r"$\gamma_v^{\ast}$",font_size=40).next_to(vl, DOWN).shift(0.3*DOWN)
        glb = SurroundingRectangle(gl, buff = .15, color = ORANGE )
        vlb = SurroundingRectangle(vl, buff = .15, color = YELLOW )
        gvlb = SurroundingRectangle(gvl, buff = .15, color = TEAL )
        self.add_fixed_in_frame_mobjects(ti, t1,t2, gl, vl, gvl, glb, vlb, gvlb)
        s2 = Surface(
            lambda u, v: np.array([
                2* np.cos(u) * np.cos(v),
                2* np.cos(u) * np.sin(v),
                2* np.sin(u)-3/2
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            fill_color=BLUE, resolution=(20, 20), fill_opacity=0.3, checkerboard_colors=[MAROON_B, MAROON_B]
        )
        def ga(t):
            return [ np.cos(t)*np.sin( PI/2 - np.sin(5*PI/6 - t)) , np.sin(t)*np.sin( PI/2 - np.sin(5*PI/6 - t)) ,   np.cos( PI/2 - np.sin(5*PI/6 - t)) ]   
        g = ParametricFunction(lambda t :  [2*ga(t)[0],2*ga(t)[1], 2*ga(t)[2]-  3/2] , t_range = [ -PI/2 , PI ] , color = ORANGE)
        v1 = Sphere([0,0,1/2], radius = 0.1)
        v1.set_color(YELLOW)
        e1 = ParametricFunction(lambda t : [ 2*np.cos(t), 2*np.sin(t),-3/2 ], t_range = [ 0 , TAU ] , color = BLACK)
        p1 =  ParametricFunction(lambda t : [ 2*np.cos(t)*np.sin( PI/2 ), 2*np.sin(t)*np.sin( PI/2 ),-3/2 ], t_range = [ -PI/2 , PI ] , color = TEAL)      
        self.add( g, s2, v1, e1, p1)
        self.wait()
        self.play(FadeOut(v1), FadeOut(e1), FadeOut(p1))
        v2 = Sphere([0,sqrt(2),-3/2 - sqrt(2)], radius = 0.1).set_color(YELLOW)
        e2 = ParametricFunction(lambda t : [ 2*np.cos(t), sqrt(2)*np.sin(t),-3/2 +  sqrt(2)*np.sin(t) ], t_range = [ 0 , TAU ] , color = BLACK)
        self.play(Create(v2), Create(e2))
        s = ValueTracker(0)
        p2 = g.copy()
        p2.set_color(TEAL)
        self.add(p2)
        def proy( x,y,z,a,b,c ):
            return [(x*a+ y*b + z*c)*a, (x*a+ y*b + z*c)*b , (x*a+ y*b + z*c)*c]
        def gs(a,b,c):
            return [a / sqrt( a**2 + b**2 + c**2 ), b / sqrt( a**2 + b**2 + c**2 ), c / sqrt( a**2 + b**2 + c**2 )]
        def past(x,y,z,a,b,c,t):
            return gs( x - t*proy(x,y,z,a,b,c)[0], y- t*proy(x,y,z,a,b,c)[1] , z - t*proy(x,y,z,a,b,c)[2])
        p2.add_updater(  lambda x: x.become( ParametricFunction(lambda t : [ 2* past( ga(t)[0], ga(t)[1], ga(t)[2], 0, 1/sqrt(2), -1/sqrt(2), s.get_value() )[0], 2*past( ga(t)[0], ga(t)[1], ga(t)[2], 0, 1/sqrt(2), -1/sqrt(2), s.get_value() )[1],2*past( ga(t)[0], ga(t)[1], ga(t)[2], 0, 1/sqrt(2), -1/sqrt(2) , s.get_value() )[2] - 3/2], t_range = [ -PI/2 , PI ] , color = TEAL) ))
        self.play(s.animate.set_value(1), run_time = 2)
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.play(FadeOut(p2), FadeOut(v2), FadeOut(e2))
        v3 = Sphere([2,0,-3/2], radius = 0.1).set_color(YELLOW)
        e3 = ParametricFunction(lambda t : [ 0, 2*np.cos(t), 2*np.sin(t)-3/2], t_range = [ 0 , TAU ] , color = BLACK)
        self.play(Create(v3), Create(e3))
        s = ValueTracker(0)
        p3 = g.copy()
        p3.set_color(TEAL)
        self.add(p3)
        p3.add_updater(  lambda x: x.become( ParametricFunction(lambda t : [ 2* past( ga(t)[0], ga(t)[1], ga(t)[2], 1, 0,0, s.get_value() )[0], 2*past( ga(t)[0], ga(t)[1], ga(t)[2], 1, 0, 0, s.get_value() )[1],2*past( ga(t)[0], ga(t)[1], ga(t)[2],  1,0,0  , s.get_value() )[2] - 3/2], t_range = [ -PI/2 , PI ] , color = TEAL) ))
        self.play(s.animate.set_value(1), run_time = 2)
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.wait()
        
   


class cr1p(Scene):
    def construct(self):
        ti = Tex(r"Theorem",font_size=40, color=BLUE).to_edge(UL).shift(0.2*DOWN)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is piecewise smooth regular, then",font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"length$(\gamma ) = \frac{1}{4 \pi} \int _{ \mathbb{S}^2} $ length$(\gamma_v^{\ast})dv. $",font_size=40).next_to(t1,DOWN)
        l = ti.get_left()
        t1.shift((t1.get_left()-l)*LEFT)
        t2.shift((t2.get_center())*LEFT)
        self.add(ti, t1,t2)
        self.wait()
        tp = Tex(r"Proof",font_size=40, color=BLUE).next_to(t2, DOWN).shift(0.2*DOWN)
        p1 = Tex(r"If $\gamma_n \to \gamma $ uniformly, then",font_size=40).next_to(tp,DOWN)
        p2 = Tex(r" $(\gamma_n)_v^{\ast}  \to \gamma_v^{\ast}$ for almost all $v \in \mathbb{S}^2$. ",font_size=40).next_to(p1,DOWN)
        p3 = Tex(r" Then",font_size=40).next_to(p2,DOWN)
        p4 = Tex(r"$ \int_{\mathbb{S}^2}$ length$(\gamma_v^{\ast}) dv $", r" $ =$", r" $\int_{\mathbb{S}^2} \liminf $ length$((\gamma_n)_v^{\ast}) dv$",font_size=40).next_to(p3,DOWN)
        p5 = Tex(r"$ \leq $", r" $ \liminf \int_{\mathbb{S}^2}$ length$((\gamma_n)_v^{\ast}) dv .  $ ",font_size=40).next_to(p4,DOWN)
        tp.shift((tp.get_left()-l)*LEFT)
        p1.shift((p1.get_left()-l)*LEFT)
        p2.shift(p2.get_center()*LEFT)
        p3.shift((p3.get_left()-l)*LEFT)
        p4.shift(p4.get_center()*LEFT)
        p5.shift((p5[0].get_left()-p4[1].get_left())*LEFT)
        self.play(Write(tp), Write(p1), Write(p2))
        self.wait()
        self.play(Write(p3), Write(p4), Write(p5))
        self.wait()


class theor(Scene):
    def construct(self):
       prop = Tex("Theorem", font_size=60, color = LIGHT_PINK).to_edge(UL).shift(0.2*RIGHT)
       v0 = prop.get_left()[0]
       intro = Tex(r"Let $\lambda : \{ $", r" $Curves$", r" $ \} \to \mathbb{R} \cup \{ \infty \}$ satisfy:", font_size=40).next_to(prop,DOWN).shift((0.1)*DOWN)
       intro.shift((v0 - intro.get_left()[0])*RIGHT)
       p0 = Tex(r"$\bullet$ If $\gamma$ is a ", r"line, ", r" $\lambda (\gamma ) = \vert \gamma (b) - \gamma (a) \vert$ ", r" .", font_size=40).next_to(intro,DOWN).shift((0.2)*DOWN)
       p0.shift((v0 - p0.get_left()[0])*RIGHT)
       b0 = SurroundingRectangle(p0[2], buff = .1, color = RED)
       p1 = Tex(r"$\bullet$ ", r"$ \lambda (\gamma_1 ) = \lambda (\gamma_2)$", r" when $\gamma_1, \gamma_2$ are reparametrizations", font_size=40).next_to(p0,DOWN).shift((0.1)*DOWN)
       p11 = Tex(r"of each other.", font_size=44).next_to(p1,DOWN).shift((0.5)*RIGHT)
       p1.shift((v0 - p1.get_left()[0])*RIGHT)
       p11.shift((v0 - p11.get_left()[0]+1/2)*RIGHT)
       b1 = SurroundingRectangle(p1[1], buff = .1, color = YELLOW)
       p2 = Tex(r"$\bullet$ ", r" $\lambda (\gamma_1 \ast \gamma_2 ) = \lambda (\gamma_1) +  \lambda (\gamma _2)$ ", r" for concatenations.", font_size=40).next_to(p11,DOWN).shift((0.1)*DOWN)
       p2.shift((v0 - p2.get_left()[0])*RIGHT)
       p2[1].shift(0.1*RIGHT)
       p2[2].shift(0.2*RIGHT)
       b2 = SurroundingRectangle(p2[1], buff = .1, color = GREEN)
       b21 = SurroundingRectangle(p2[1], buff = .2, color = BLUE)
       p3 = Tex(r"$\bullet$ ", r"$ \lambda (T \circ \gamma ) = \lambda (\gamma) $ ", r" for $T$ any isometry.", font_size=40).next_to(p2,DOWN).shift((0.1)*DOWN)
       p3.shift((v0 - p3.get_left()[0])*RIGHT)
       b3 = SurroundingRectangle(p3[1], buff = .1, color = BLUE)
       p4 = Tex(r"$\bullet$ If $\gamma_n \to \gamma$", " pointwise, ", r" $ \lambda( \gamma ) \leq  \liminf  \lambda (\gamma_n) $ ", r" .", font_size=40).next_to(p3,DOWN).shift((0.1)*DOWN)
       p4.shift((v0 - p4.get_left()[0])*RIGHT)
       b4 = SurroundingRectangle(p4[2], buff = .1, color = ORANGE)
       con = Tex(r"Then ", r"$ \lambda (\gamma ) = $ length$(\gamma )$", r" for all ", r"curve", r" $\gamma$.", font_size=40).next_to(p4,DOWN).shift((0.2)*DOWN)
       con.shift((v0 - con.get_left()[0])*RIGHT)
       bf = SurroundingRectangle(con[1], buff = .1, color = LIGHT_PINK)
       self.play(Write(prop),Write(intro), Write(p0), Create(b0), Write(p1),Write(p11), Create(b1), Write(p2), Create(b2), Create(b21), Write(p3), Create(b3),Write(p4), Create(b4), Write(con), Create(bf) )
       self.wait()
       intros = Tex(r"Let $\lambda : \{ $", r"$pwsr$-$Spherical$-$Curves$", r" $ \} \to \mathbb{R} \cup \{ \infty \}$ satisfy:", font_size=40).next_to(prop,DOWN).shift((0.1)*DOWN)
       cons = Tex(r"Then ", r"$ \lambda (\gamma ) = $ length$(\gamma )$", r" for all ", r"pwsr spherical curve", r" $\gamma$.", font_size=40).next_to(p4,DOWN).shift((0.2)*DOWN)
       intros.shift((v0 - intros.get_left()[0])*RIGHT)
       cons.shift((v0 - cons.get_left()[0])*RIGHT)
       self.play(intro[1].animate.become(intros[1]), intro[2].animate.become(intros[2]), con[3].animate.become(cons[3]), con[4].animate.become(cons[4]) )
       self.wait()
       p4s = Tex(r"$\bullet$ If $\gamma_n \to \gamma$", " uniformly, ", r" $ \lambda( \gamma ) \leq  \liminf  \lambda (\gamma_n) $ ", r" .", font_size=40).next_to(p3,DOWN).shift((0.1)*DOWN)
       p4s.shift((v0 - p4s.get_left()[0])*RIGHT)
       self.wait()
       p0s = Tex(r"$\bullet$ If $\gamma$ is a ", r"``line'', ", r" $\lambda (\gamma ) = $ length$(\gamma)$ ", r" .", font_size=40).next_to(intro,DOWN).shift((0.2)*DOWN)
       p0s.shift((v0 - p0s.get_left()[0])*RIGHT)
       b0s = SurroundingRectangle(p0s[2], buff = .1, color = RED)
       self.play(p0[1].animate.become(p0s[1]),p0[2].animate.become(p0s[2]),p0[3].animate.become(p0s[3]), b0.animate.become(b0s))
       self.wait()
       b4s = SurroundingRectangle(p4s[2], buff = .1, color = ORANGE)
       self.play(p4[1].animate.become(p4s[1]), p4[2].animate.become(p4s[2]), p4[3].animate.become(p4s[3]), b4.animate.become(b4s))
       self.wait()





class thln(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=0 * DEGREES, distance = 8)
        prop = Tex("Theorem", font_size=60, color = LIGHT_PINK).to_edge(UL).shift(0.2*RIGHT)
        v0 = prop.get_left()[0]
        intros = Tex(r"Let $\lambda : \{ $", r"$pwsr$-$Spherical$-$Curves$", r" $ \} \to \mathbb{R} \cup \{ \infty \}$ satisfy:", font_size=40).next_to(prop,DOWN).shift((0.1)*DOWN)
        intros.shift((v0 - intros.get_left()[0])*RIGHT)
        p0s = Tex(r"$\bullet$ If $\gamma$ is a ", r"``line'', ", r" $\lambda (\gamma ) = $ length$(\gamma)$ ", r" .", font_size=40).next_to(intros,DOWN).shift((0.2)*DOWN)
        p0s.shift((v0 - p0s.get_left()[0])*RIGHT)
        b0s = SurroundingRectangle(p0s[2], buff = .1, color = RED)
        p1 = Tex(r"$\bullet$ ", r"$ \lambda (\gamma_1 ) = \lambda (\gamma_2)$", r" when $\gamma_1, \gamma_2$ are reparametrizations", font_size=40).next_to(p0s,DOWN).shift((0.1)*DOWN)
        p11 = Tex(r"of each other.", font_size=44).next_to(p1,DOWN).shift((0.5)*RIGHT)
        p1.shift((v0 - p1.get_left()[0])*RIGHT)
        p11.shift((v0 - p11.get_left()[0]+1/2)*RIGHT)
        b1 = SurroundingRectangle(p1[1], buff = .1, color = YELLOW)
        p2 = Tex(r"$\bullet$ ", r" $\lambda (\gamma_1 \ast \gamma_2 ) = \lambda (\gamma_1) +  \lambda (\gamma _2)$ ", r" for concatenations.", font_size=40).next_to(p11,DOWN).shift((0.1)*DOWN)
        p2.shift((v0 - p2.get_left()[0])*RIGHT)
        p2[1].shift(0.1*RIGHT)
        p2[2].shift(0.2*RIGHT)
        b2 = SurroundingRectangle(p2[1], buff = .1, color = GREEN)
        b21 = SurroundingRectangle(p2[1], buff = .2, color = BLUE)
        p3 = Tex(r"$\bullet$ ", r"$ \lambda (T \circ \gamma ) = \lambda (\gamma) $ ", r" for $T$ any isometry.", font_size=40).next_to(p2,DOWN).shift((0.1)*DOWN)
        p3.shift((v0 - p3.get_left()[0])*RIGHT)
        b3 = SurroundingRectangle(p3[1], buff = .1, color = BLUE)
        p4 = Tex(r"$\bullet$ If $\gamma_n \to \gamma$", " pointwise, ", r" $ \lambda( \gamma ) \leq  \liminf  \lambda (\gamma_n) $ ", r" .", font_size=40).next_to(p3,DOWN).shift((0.1)*DOWN)
        p4.shift((v0 - p4.get_left()[0])*RIGHT)
        b4 = SurroundingRectangle(p4[2], buff = .1, color = ORANGE)
        cons = Tex(r"Then ", r"$ \lambda (\gamma ) = $ length$(\gamma )$", r" for all ", r"pwsr spherical curve", r" $\gamma$.", font_size=40).next_to(p4,DOWN).shift((0.2)*DOWN)
        cons.shift((v0 - cons.get_left()[0])*RIGHT)
        bf = SurroundingRectangle(cons[1], buff = .1, color = LIGHT_PINK)
        self.add_fixed_in_frame_mobjects(prop,intros,p0s,b0s,p1,p11,b1,p2,b2,b21,p3, b3,p4,b4,cons,bf )
        self.wait()
        s2 = Surface(
            lambda u, v: np.array([
                3*np.cos(u) * np.cos(v)/4,
                3*np.cos(u) * np.sin(v)/4 + 4.5,
                3*np.sin(u)/4
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            fill_color=BLUE, resolution=(20, 20), fill_opacity=0.6, checkerboard_colors=[MAROON_B, MAROON_B]
        )
        gam = ParametricFunction(lambda t : [  3 *np.sin (t)/ sqrt(32), 4.5 - 3*np.sin (t)/sqrt(32), 3*np.cos(t)/4 ], t_range = [ 0 , 2*PI ] , color = YELLOW)
        self.play(Create(s2))
        self.play(Create(gam))
        self.wait()
        self.play(FadeOut(gam))
        gam1 = ParametricFunction(lambda t : [ 3*np.cos(t)/(4*sqrt(2)) + 3*np.sin(t)/(4*sqrt(6)), -3*np.cos(t)/(4*sqrt(2)) + 3*np.sin(t)/(4*sqrt(6)) + 4.5, sqrt(6) *np.sin(t)/4  ], t_range = [ -PI/3 , 5*PI/3 ] , color = YELLOW)
        self.play(Create(gam1))
        self.wait()
        self.play(FadeOut(gam1), FadeOut(s2))
        
              


class cr2pre(Scene):
    def construct(self):
        ti = Tex(r"Theorem",font_size=40, color=BLUE).to_edge(UL)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is piecewise smooth regular, then",font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"length$(\gamma ) = \frac{1}{4}  \int _{ \mathbb{S}^2}  \#  ( \gamma \cap v^{\perp} )  dv $",font_size=40).next_to(t1,DOWN)
        l = ti.get_left()
        t1.shift((t1.get_left()-l)*LEFT)
        t2.shift((t2.get_center())*LEFT)
        self.play(Create(ti),Create(t1), Create(t2))
        self.wait()

class cr2(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=30 * DEGREES, distance = 8)
        ti = Tex(r"Theorem",font_size=40, color=BLUE).to_edge(UL)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is piecewise smooth regular, then",font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"length$(\gamma ) = \frac{1}{4}  \int _{ \mathbb{S}^2}  \#   ( \gamma \cap v^{\perp} )  dv $",font_size=40).next_to(t1,DOWN)
        l = ti.get_left()
        t1.shift((t1.get_left()-l)*LEFT)
        t2.shift((t2.get_center())*LEFT)
        self.add_fixed_in_frame_mobjects(ti, t1,t2)
        s2 = Surface(
            lambda u, v: np.array([
                2* np.cos(u) * np.cos(v),
                2* np.cos(u) * np.sin(v),
                2* np.sin(u)-3/2
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            fill_color=BLUE, resolution=(20, 20), fill_opacity=0.3, checkerboard_colors=[MAROON_B, MAROON_B]
        )
        g = ParametricFunction(lambda t : [ 2*np.cos(t)*np.sin(PI/2 - np.sin(2*t)/4), 2*np.sin(t)*np.sin(PI/2 - np.sin(2*t)/4),-3/2 + 2*np.cos(PI/2 - np.sin(2*t)/4) ], t_range = [ 0 , PI ] , color = ORANGE)
        self.play(Create(s2), Create(g))
        gl = Tex(r"$\gamma$",font_size=40).shift(4*LEFT)
        vl = Tex(r"$v$",font_size=40).next_to(gl, DOWN).shift(0.3*DOWN)
        gvl = Tex(r"$\gamma \cap v^{\perp}$",font_size=40).next_to(vl, DOWN).shift(0.3*DOWN)
        glb = SurroundingRectangle(gl, buff = .15, color = ORANGE )
        vlb = SurroundingRectangle(vl, buff = .15, color = YELLOW )
        gvlb = SurroundingRectangle(gvl, buff = .075, color = PINK )
        self.add_fixed_in_frame_mobjects(ti, t1,t2, gl, vl, gvl, glb, vlb, gvlb)
        v1 =  Sphere([0,0,1/2], radius = 0.1).set_color(YELLOW)
        e1 = ParametricFunction(lambda t : [ 2*np.cos(t), 2*np.sin(t), -3/2 ], t_range = [ 0 , TAU ] , color = BLACK)
        self.play(Create(v1), Create(e1))
        self.wait()
        p = Sphere([2,0,-3/2], radius = 0.1).set_color(PINK)
        q = Sphere([0,2,-3/2], radius = 0.1).set_color(PINK)
        r = Sphere([-2,0,-3/2], radius = 0.1).set_color(PINK)
        self.play(Create(p), Create(q), Create(r))
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.play(FadeOut(v1), FadeOut(e1), FadeOut(p), FadeOut(q), FadeOut(r))
        self.wait()
        v2 =  Sphere([sqrt(2),0,sqrt(2) -3/2], radius = 0.1).set_color(YELLOW)
        e2 = ParametricFunction(lambda t : [ sqrt(2)*np.cos(t), 2* np.sin(t), -3/2 - sqrt(2)*np.cos(t) ], t_range = [ 0 , TAU ] , color = BLACK)
        self.play(Create(v2), Create(e2))
        self.wait()
        p = Sphere([0,2,-3/2], radius = 0.1).set_color(PINK)
        self.play(Create(p))
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/6)
        self.wait(4)
        self.stop_ambient_camera_rotation()
        self.wait()


class cr2p(Scene):
    def construct(self):
        ti = Tex(r"Theorem",font_size=40, color=BLUE).to_edge(UL)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is piecewise smooth regular, then",font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"length$(\gamma ) = \frac{1}{4}  \int _{ \mathbb{S}^2}  \#   ( \gamma \cap v^{\perp} )  dv $",font_size=40).next_to(t1,DOWN)
        l = ti.get_left()
        t1.shift((t1.get_left()-l)*LEFT)
        t2.shift((t2.get_center())*LEFT)
        self.add(ti, t1,t2)
        all = Group(ti,t1,t2)
        self.play(all.animate.shift(0.2*DOWN))
        self.wait()
        tp = Tex(r"Proof",font_size=40, color=BLUE).next_to(t2, DOWN).shift(0.2*DOWN)
        p1 = Tex(r"Pick a sequence $\gamma _ n : [a,b] \to \mathbb{S}^2$ of inscribed ``piecewise linear'' curves ",font_size=40).next_to(tp,DOWN)
        p11 = Tex(r"with $\gamma_n \to \gamma$, and such that $\gamma_{n+1} $ refines $\gamma_n$.",font_size=40).next_to(p1,DOWN)
        p2 = Tex(r"By Sard's Theorem, $\gamma$ and $\gamma_n$ are transversal to $v^{\perp} $ for almost all $v \in \mathbb{S}^2$.",font_size=40).next_to(p11,DOWN)
        p3 = Tex(r"For those $v$, $\# (\gamma _ n \cap v^{\perp} ) \nearrow \# (\gamma \cap v^{\perp})$.",font_size=40).next_to(p2,DOWN)
        p4 = Tex(r"By Levy monotone convence theorem, ",font_size=40).next_to(p3,DOWN)
        p5 = Tex(r"$\int _{ \mathbb{S}^2}  \#   ( \gamma_n \cap v^{\perp} )  dv \to \int _{ \mathbb{S}^2}  \#   ( \gamma \cap v^{\perp} )  dv $.",font_size=40).next_to(p4,DOWN)
        tp.shift((tp.get_left()-l)*LEFT)
        p1.shift((p1.get_left()-l)*LEFT)
        p11.shift((p11.get_left()-l)*LEFT)
        p2.shift((p2.get_left()-l)*LEFT)
        p3.shift((p3.get_left()-l)*LEFT)
        p4.shift((p4.get_left()-l)*LEFT)
        p5.shift(p5.get_center()*LEFT)
        self.play(Write(tp), Write(p1), Write(p11))
        self.wait()
        self.play(Write(p2))
        self.wait()
        self.play(Write(p3))
        self.wait()
        self.play(Write(p4), Write(p5))
        self.wait()



class cr12(Scene):
    def construct(self):
        ti = Tex(r"Theorem",font_size=40, color=BLUE).to_edge(UL).shift(0.8*DOWN)
        t1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is piecewise smooth regular, then",font_size=40).next_to(ti,DOWN)
        t2 = Tex(r"length$(\gamma ) = \frac{1}{4 \pi} \int _{ \mathbb{S}^2} $ length$(\gamma_v^{\ast})dv. $",font_size=40).next_to(t1,DOWN)
        l = ti.get_left()
        t1.shift((t1.get_left()-l)*LEFT)
        t2.shift((t2.get_center())*LEFT)
        pi = Tex(r"Theorem",font_size=40, color=BLUE).next_to(t2, DOWN).shift(0.3*DOWN)
        p1 = Tex(r"If $\gamma : [a,b] \to \mathbb{S}^2$ is piecewise smooth regular, then",font_size=40).next_to(pi,DOWN)
        p2 = Tex(r"length$(\gamma ) = \frac{1}{4}  \int _{ \mathbb{S}^2}  \#  ( \gamma \cap v^{\perp} )  dv $",font_size=40).next_to(p1,DOWN)
        pi.shift((pi.get_left()-l)*LEFT)
        p1.shift((p1.get_left()-l)*LEFT)
        p2.shift((p2.get_center())*LEFT)
        self.play(Write(ti),Write(t1), Write(t2), Write(pi),Write(p1), Write(p2))
        self.wait()


        
        
        

class ftest(Scene):
    def construct(self):
        def proy( x,y,z,a,b,c ):
            return [(x*a+ y*b + z*c)*a, (x*a+ y*b + z*c)*b , (x*a+ y*b + z*c)*c]
        def gs(a,b,c):
            return [a / sqrt( a^2 + b^2 + c^2 ), b / sqrt( a^2 + b^2 + c^2 ), c / sqrt( a^2 + b^2 + c^2 )]
        Line1 = Line([0,0,0], [1,0,0])
        Line2 = Line([0,1,0], gs(1,1,0))
        self.play(Create(Line1), Create(Line2))
        


class centertest(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES, distance = 8)
        p = Sphere([0,0,0], radius = 1)
        self.play(Create(p))
        self.play(p.animate.shift(RIGHT).set_color(GREEN))
        self.play(p.animate.set_center([0,0,0]))














class test(ThreeDScene):
    def construct(self):
        self.renderer.camera.light_source.move_to(0) # changes the source of the light
        self.set_camera_orientation(phi=70 * DEGREES, theta=45 * DEGREES, distance = 3)
        axes = ThreeDAxes()
        s22 = Surface(
            lambda u, v: np.array([
                1.5* np.cos(u) * np.cos(v),
                1.5* np.cos(u) * np.sin(v),
                1.5* np.sin(u)
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            fill_color=BLUE, resolution=(20, 20), fill_opacity=0.4
        )
        gam = ParametricFunction(lambda t : [sqrt(1.25)*np.cos(t),sqrt(1.25)*np.sin(t),1 ], t_range = [0 , TAU ] , color = YELLOW)
        self.add(axes)
        self.play(Create(s22))
        self.play(Create(gam))
        self.wait(2)