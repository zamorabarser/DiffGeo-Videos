from this import d
from tkinter import E
from manim import *
from numpy import sqrt
import math



class ti(Scene):
    def construct(self):
        t1 = Text("Calculus on surfaces", font_size=60)
        self.play(Write(t1))
        self.wait()        


class obs(ThreeDScene):
    def construct(self):
        thm = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift(0.4*RIGHT + 0.2*DOWN)
        l = thm.get_left()[0]
        thm1 = Tex(r"If $\phi : $ ", r"$U $ ", r"$ \to $ ", r"$\Sigma $", r", $\psi : $ ", r"$V$ ", r"$ \to $ ", r"$\Sigma$", r" are charts, then", font_size = 34).next_to(thm, DOWN)
        thm1.shift((l-thm1.get_left()[0])*RIGHT)
        thm2 = Tex(r"$(\psi^{-1} \circ \phi ) :$ ", r"$ \phi^{-1}( \psi (V)) $ ", r"$\to $ ", r"$\psi^{-1}(\phi(U))$ ", r"$  \subset V$ is smooth.", font_size = 34).next_to(thm1, DOWN)
        thm2.shift((l-thm2.get_left()[0])*RIGHT)
        uli = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(thm1[1], DOWN).shift(0.15*UP)
        sli1 = Line([0,0,0], [0.3,0,0], color = PURPLE).next_to(thm1[3], DOWN).shift(0.15*UP)
        vli = Line([0,0,0], [0.3,0,0], color = GREEN).next_to(thm1[5], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.3,0,0], color = PURPLE).next_to(thm1[7], DOWN).shift(0.15*UP)
        phs = Line([0,0,0], [1.4,0,0], color = YELLOW).next_to(thm2[1], DOWN).shift(0.15*UP)
        psh = Line([0,0,0], [1.4,0,0], color = ORANGE).next_to(thm2[3], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-1,3,10], x_length = 5, y_length = 5, z_length = 2 )
        ax.shift([0,0,0])
        axu = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-2,3,10], x_length = 3, y_length = 3, z_length = 0.1 , tips = False )
        axu.shift([2,-2*sqrt(3),-2.5])
        axv = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-2,3,10], x_length = 3, y_length = 3, z_length = 0.1 , tips = False  )
        axv.shift([-2,2*sqrt(3),-2.5])
        phi = Arrow([-3,-1.5,0], [-1.5,-0.5,0])
        phil = Tex(r"$\phi$", font_size = 34).next_to(phi, LEFT).shift(0.2*UP)
        psi = Arrow([3,-1.5,0], [1.5,-0.5,0])
        psil = Tex(r"$\psi$", font_size = 34).next_to(psi, RIGHT).shift(0.2*UP)
        cc = Arrow([-1.5,-2.5,0], [1.5,-2.5,0])
        ccl = Tex(r"$\psi ^{-1} \circ \phi$", font_size = 34).next_to(cc, DOWN).shift(0.2*UP)
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( u , v , 2**( -(u**2 + v**2)/4)),
            u_range = [-3,3],
            v_range = [-3,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        sigu = Surface(
            lambda u, v: ax.c2p( u , v , 2**( -(u**2 + v**2)/4)),
            u_range = [-1,2],
            v_range = [-2,1],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        sigv = Surface(
            lambda u, v: ax.c2p( u , v , 2**( -(u**2 + v**2)/4)),
            u_range = [-2,1],
            v_range = [-1,2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        siguv = Surface(
            lambda u, v: ax.c2p( u , v , 2**( -(u**2 + v**2)/4)),
            u_range = [-1,1],
            v_range = [-1,1],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = 0.7
        )
        u = Surface(
            lambda u, v: axu.c2p( u , v , 0),
            u_range = [-1,2],
            v_range = [-2,1],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        v = Surface(
            lambda u, v: axv.c2p( u , v , 0),
            u_range = [-2,1],
            v_range = [-1,2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        uuv =  Surface(
            lambda u, v: axu.c2p( u , v , 0),
            u_range = [-1,1],
            v_range = [-1,1],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = 0.7
        )
        vvu =  Surface(
            lambda u, v: axv.c2p( u , v , 0),
            u_range = [-1,1],
            v_range = [-1,1],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7
        )
        siguv2 = siguv.copy()
        ul = Line(axu.c2p(-1,-2,0), axu.c2p(2,-2,0), color = BLUE).append_points(Line(axu.c2p(2,-2,0), axu.c2p(2,1,0), color = BLUE).points )
        ul.append_points(Line(axu.c2p(2,1,0), axu.c2p(-1,1,0), color = BLUE).points).append_points(Line(axu.c2p(-1,1,0), axu.c2p(-1,-2,0), color = BLUE).points)
        vl = Line(axv.c2p(1,2,0), axv.c2p(-2,2,0), color = GREEN).append_points(Line(axv.c2p(-2,2,0), axv.c2p(-2,-1,0), color = GREEN).points )
        vl.append_points(Line(axv.c2p(-2,-1,0), axv.c2p(1,-1,0), color = GREEN).points).append_points(Line(axv.c2p(1,-1,0), axv.c2p(1,2,0), color = GREEN).points)
        uu = u.copy()
        vv = v.copy()
        self.add_fixed_in_frame_mobjects(thm, thm1, thm2, phi, psi, cc, phil, psil, ccl, uli, vli, sli1, sli2, phs, psh)
        self.remove(thm, thm1, thm2, phi, psi, cc, phil, psil, ccl, uli, vli, sli1, sli2, phs, psh)
        self.play(Create(thm), Create(thm1), Create(thm2), Create(uli), Create(vli), Create(sli1), Create(sli2), Create(phs), Create(psh))
        self.wait()
        self.play(Create(ax), Create(axu), Create(axv))
        self.wait()
        self.play(Create(sig), Create(u), Create(v), Create(phi), Create(phil), Create(psi), Create(psil), Create(cc), Create(ccl))
        self.wait()
        self.add(uu)
        self.play(Transform(uu,sigu))
        self.wait()
        self.add(vv)
        self.play(Transform(vv,sigv))
        self.wait()
        self.play(FadeOut(sig))
        self.wait()
        self.play(Create(ul), FadeOut(u), Create(vl), FadeOut(v))
        self.wait()
        self.add(siguv)
        self.play(Transform(siguv, uuv))
        self.wait()
        self.add(uuv)
        self.play(Transform(uuv, siguv2))
        self.play(Transform(uuv, vvu))
        self.wait()


class difdef(ThreeDScene):
    def construct(self):
        dif = Tex(r"A continuous function $f:$ ",r"$ \Sigma $ ", r"$ \to \mathbb{R}$ is ", r"smooth", r" if for any chart", font_size = 34).to_edge(UL).shift(0.4*RIGHT + 1.2*DOWN)
        l = dif.get_left()[0]
        dif[3].set_color(BLUE)
        dif1 = Tex(r"$\phi : $ ", r"$ U $ ", r"$ \to $ ", r"$\Sigma $, ", r"the composition $f \circ \phi : $ ", r"$U $ ", r"$ \to \mathbb{R}$ is smooth.", font_size = 34).next_to(dif, DOWN)
        dif1.shift((l-dif1.get_left()[0])*RIGHT)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(dif[1], DOWN).shift(0.15*UP)
        uli = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(dif1[1], DOWN).shift(0.15*UP)
        sli1 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(dif1[3], DOWN).shift(0.15*UP)
        uli1 = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(dif1[5], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-1,3,10], x_length = 4.5, y_length = 4.5, z_length = 1.8 )
        ax.shift([-1/5,sqrt(3)/5,-1])
        axu = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-2,3,10], x_length = 3, y_length = 3, z_length = 0.1 , tips = False )
        axu.shift([2*(1.2),-2*sqrt(3)*(1.2),-1])
        phi = Arrow([-3.5,-0.5,0], [-1.5,-0.5,0])
        phil = Tex(r"$\phi$", font_size = 34).next_to(phi, UP).shift(0.2*DOWN)
        f = Arrow([2.5,-0.5,0], [4.5,-0.5,0])
        fl = Tex(r"$f$", font_size = 34).next_to(f, UP).shift(0.2*DOWN)
        r = ThreeDAxes(x_range = [-1,1,10], y_range = [-1,1,10], z_range = [-5,5,1], x_length = 0.1, y_length = 0.1, z_length = 3, tips = False )
        r.shift([-2*(1.2),2*sqrt(3)*(1.2),-1])
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( u , v , np.sin(2*u)/2 + np.cos(2*v)/2),
            u_range = [-3,3],
            v_range = [-3,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        u = Surface(
            lambda u, v: axu.c2p( u , v , 0),
            u_range = [-2,2],
            v_range = [-2,2],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        self.add_fixed_in_frame_mobjects(dif, dif1, phi, phil, uli, sli, uli1, sli1, f, fl)
        self.remove(dif, dif1, phi, phil, uli, sli, uli1, sli1, f, fl)
        self.play(Create(uli),  Create(sli), Create(dif), Create(dif1), Create(uli1), Create(sli1))
        self.wait()
        self.play(Create(sig), Create(f), Create(fl), Create(ax), Create(r))
        self.wait()
        self.play(Create(axu), Create(u), Create(phi), Create(phil))
        self.wait()
        self.play(FadeOut(uli), FadeOut(sli),FadeOut(dif), FadeOut(dif1),FadeOut(uli1), FadeOut(sli1))
        self.wait()

class difwell(ThreeDScene):
    def construct(self):
        dif = Tex(r"If $\phi : $ ", r"$ U $ ", r"$\to  $ ", r"$\Sigma$, ", r" $\psi :$ ", r"$ V$ ", r"$ \to $ ", r"$ \Sigma$ ", r" are charts with $\phi(U) = \psi(V)$, then", font_size = 34).to_edge(UL).shift(0.4*RIGHT + 0.2*DOWN)
        l = dif.get_left()[0]
        dif1 = Tex(r"$\bullet$ $f \circ \phi = (f \circ \psi) \circ (\psi^{-1} \circ \phi ) :$ ", r"$ U$ ", r"$ \to \mathbb{R}$.", font_size = 34).next_to(dif, DOWN)
        dif1.shift((l-dif1.get_left()[0])*RIGHT)
        dif2 = Tex(r"$\bullet$ $f \circ \psi = (f \circ \phi ) \circ ( \phi^{-1} \circ \psi ): $ ", r"$ V$ ", r"$ \to \mathbb{R}$.", font_size = 34).next_to(dif1, DOWN)
        dif2.shift((l-dif2.get_left()[0])*RIGHT)
        ul1 = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(dif[1], DOWN).shift(0.15*UP)
        ul2 = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(dif1[1], DOWN).shift(0.15*UP)
        vl1 = Line([0,0,0], [0.3,0,0], color = GREEN).next_to(dif[5], DOWN).shift(0.15*UP)
        vl2 = Line([0,0,0], [0.3,0,0], color = GREEN).next_to(dif2[1], DOWN).shift(0.15*UP)
        sl1 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(dif[3], DOWN).shift(0.15*UP)      
        sl2 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(dif[7], DOWN).shift(0.15*UP)        
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-1,3,10], x_length = 4.5, y_length = 4.5, z_length = 1.8 )
        ax.shift([-1/5,sqrt(3)/5,-1])
        axu = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-2,3,10], x_length = 3, y_length = 3, z_length = 0.1 , tips = False )
        axu.shift([2*(1.2),-2*sqrt(3)*(1.2),-1])
        axv = axu.copy()
        axv.shift([0,0,-1])
        axv = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-2,3,10], x_length = 3, y_length = 3, z_length = 0.1 , tips = False )
        axv.shift([2*(1.2),-2*sqrt(3)*(1.2),-2])
        phi = Arrow([-3.5,-0.5,0], [-1.5,-0.5,0])
        phip = Arrow( [-3.5,0,0], [-1.5,-0.5,0])
        psi = Arrow([-3.5,-1.5,0], [-1.5,-1,0])
        phil = Tex(r"$\phi$", font_size = 34).next_to(phi, UP).shift(0.2*DOWN)
        philp = Tex(r"$\phi$", font_size = 34).next_to(phip, UP).shift(0.2*DOWN)
        psil = Tex(r"$\psi$", font_size = 34).next_to(psi, DOWN).shift(0.2*UP)
        f = Arrow([2.5,-0.5,0], [4.5,-0.5,0])
        fl = Tex(r"$f$", font_size = 34).next_to(f, UP).shift(0.2*DOWN)
        r = ThreeDAxes(x_range = [-1,1,10], y_range = [-1,1,10], z_range = [-5,5,1], x_length = 0.1, y_length = 0.1, z_length = 3, tips = False )
        r.shift([-2*(1.2),2*sqrt(3)*(1.2),-1])
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( u , v , np.sin(2*u)/2 + np.cos(2*v)/2),
            u_range = [-3,3],
            v_range = [-3,3],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        u = Surface(
            lambda u, v: axu.c2p( u , v , 0),
            u_range = [-2,2],
            v_range = [-2,2],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        v = Surface(
            lambda u, v: axv.c2p( u , v , 0),
            u_range = [-2,2],
            v_range = [-2,2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        self.add(ax, sig, r, axu, u, dif2)
        self.add_fixed_in_frame_mobjects(dif, dif1, phi, phil, f, fl, psil, psi, phip, philp, dif2, ul1, ul2, vl1, vl2, sl1, sl2)
        self.remove(dif, dif1, dif2, psil, psi, phip, philp, ul1, ul2, vl1, vl2, sl1, sl2)
        self.play( Create(dif), Create(dif1), Create(dif2), Create(ul1), Create(ul2), Create(vl1), Create(vl2), Create(sl1), Create(sl2) )
        self.wait()
        self.play(u.animate.shift([0,0,1.1]), axu.animate.shift([0,0,1.1]), Create(axv), Create(v), Create(psi), Create(psil), Transform(phi, phip), Transform(phil, philp))
        self.wait()


class tp(ThreeDScene):
    def construct(self):
        tan = Tex(r"For $\phi:$ ",r"$ U $ ", r"$ \to $ ", r"$\Sigma$ ", r"a chart, and ", r"$p$ ", r"$ = \phi ($ ", r"$q$ ", r"$)$,", font_size = 34).to_edge(UL).shift(0.4*RIGHT + 1.2*DOWN)
        l = tan.get_left()[0]
        tan1 = Tex(r"$ T_p\Sigma$ ", r"$ : =  \langle  $ ", r"$\frac{ \partial \phi }{ \partial u} (q) $ ,", r"$\frac{\partial \phi }{\partial v} (q) $ ", r"$ \rangle$.", font_size = 34).next_to(tan, DOWN)
        tan1.shift((l-tan1.get_left()[0])*RIGHT)
        uli = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(tan[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(tan[3], DOWN).shift(0.15*UP)
        pu = Line([0,0,0], [0.6,0,0], color = YELLOW).next_to(tan1[2], DOWN).shift(0.15*UP)
        pv = Line([0,0,0], [0.6,0,0], color = YELLOW).next_to(tan1[3], DOWN).shift(0.15*UP)
        pla = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(tan[5], DOWN).shift(0.15*UP)
        qla = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(tan[7], DOWN).shift(0.15*UP)
        tp = Line([0,0,0], [0.6,0,0], color = GREEN).next_to(tan1[0], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-2/5,6/5,10], x_length = 3, y_length = 3, z_length = 1.6 , tips = False )
        ax.shift([-(1.8),(1.8)*sqrt(3),-1])
        axu = ThreeDAxes(x_range = [-2,2,10], y_range = [-2,2,10], z_range = [-2,3,10], x_length = 3, y_length = 3, z_length = 0.1 , tips = False )
        axu.shift([(1.8),-sqrt(3)*(1.8),-1])
        phi = CurvedArrow([-1.5,-0.3,0], [1.5,-0.3,0], radius = -5)
        phil = Tex(r"$\phi$", font_size = 34).next_to(phi, UP)
        q = Sphere(radius = 0.05)
        q.set_color(ORANGE).shift(axu.c2p(0,PI/6,0))
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0,1/2,sqrt(3)/2 ))
        e1 = Arrow( axu.c2p(-0.2,PI/6,0 ), axu.c2p(1.2,PI/6,0), color = YELLOW )
        e2 = Arrow( axu.c2p(0,PI/6-0.2,0 ), axu.c2p(0,1.2+PI/6,0), color = YELLOW )
        e3 = Arrow( ax.c2p(-0.2,1/2,sqrt(3)/2 ), ax.c2p(1.2,1/2,sqrt(3)/2 ), color = YELLOW )
        e4 = Arrow( ax.c2p(0,1/2 - sqrt(3)/10 ,sqrt(3)/2 + 1/10 ), ax.c2p(0,1/2 + sqrt(3)*(1.2)/2,sqrt(3)/2 - 0.6), color = YELLOW )
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( np.cos(v)*np.sin(u), np.sin(v) , np.cos(v)*np.cos(u)  ),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        u = Surface(
            lambda u, v: axu.c2p( u , v , 0),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        tps = Surface(
            lambda u, v: ax.c2p( u , sqrt(3)*v/2 + 1/2 , sqrt(3)/2 - v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        self.add_fixed_in_frame_mobjects(tan, tan1, phi, phil, uli, sli, pu, pv, pla, qla, tp)
        self.remove(tan, tan1, phi, phil, uli, sli, pu, pv, pla, qla, tp)
        self.play(Create(uli),  Create(sli), Create(tan), Create(tan1), Create(pu), Create(pv), Create(pla), Create(qla), Create(tp))
        self.wait()
        self.play(Create(ax), Create(axu), Create(sig), Create(u), Create(phi), Create(phil), Create(q), Create(p))
        self.wait()
        self.play(Create(e1), Create(e2), Create(e3), Create(e4))
        self.wait()
        self.play(Create(tps))
        self.wait()





class balance(ThreeDScene):
    def construct(self):
        ax = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-2/5,6/5,10], x_length = 6, y_length = 6, z_length = 3.2 , tips = False )
        axs = ax.copy()
        ax.shift([(7/4),-(7/4)*sqrt(3),-1])
        axs.shift([-(7/4),(7/4)*sqrt(3),-1])
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0,1/2,sqrt(3)/2 ))
        zero = Sphere(radius = 0.05)
        zero.set_color(ORANGE).shift(axs.c2p(0,0,0))
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( np.cos(v)*np.sin(u), np.sin(v) , np.cos(v)*np.cos(u)  ),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tp = Surface(
            lambda u, v: ax.c2p( u , sqrt(3)*v/2 + 1/2 , sqrt(3)/2 - v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        tpp = Surface(
            lambda u, v: axs.c2p( u , sqrt(3)*v/2 ,- v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        self.play(Create(ax), Create(sig), Create(p), Create(tp))
        self.wait()
        self.play(Create(axs), Create(tpp), Create(zero))
        self.wait()
        t = ValueTracker(0)
        tp.add_updater( lambda x: x.become( 
            Surface(
                lambda u, v: ax.c2p( np.cos(t.get_value())*u + np.sin(t.get_value())*(sqrt(3)*v/2 + 1/2), np.cos(t.get_value())*(sqrt(3)*v/2 + 1/2) - np.sin(t.get_value())*u , sqrt(3)/2 - v/2 ),
                u_range = [ -1, 1 ],
                v_range = [ -1, 1 ],
                checkerboard_colors = [GREEN, GREEN],
                fill_opacity = 0.7
            )   
        ))
        tpp.add_updater( lambda x: x.become( 
            Surface(
                lambda u, v: axs.c2p( np.cos(t.get_value())*u + np.sin(t.get_value())*(sqrt(3)*v/2), np.cos(t.get_value())*(sqrt(3)*v/2) - np.sin(t.get_value())*u , - v/2 ),
                u_range = [ -1, 1 ],
                v_range = [ -1, 1 ],
                checkerboard_colors = [GREEN, GREEN],
                fill_opacity = 0.7
            )   
        ))
        p.add_updater(lambda x: x.become(Sphere(radius = 0.05).set_color(ORANGE).shift(ax.c2p(np.sin(t.get_value())/2,np.cos(t.get_value())/2,sqrt(3)/2 ))))
        self.play(t.animate.set_value(PI/2), rate_func=rate_functions.ease_in_sine)
        self.wait()


class balance2(ThreeDScene):
    def construct(self):
        ax = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-2/5,6/5,10], x_length = 6, y_length = 6, z_length = 3.2 , tips = False )
        axs = ax.copy()
        ax.shift([(7/4),-(7/4)*sqrt(3),-1])
        axs.shift([-(7/4),(7/4)*sqrt(3),-1])
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0,1/2,sqrt(3)/2 ))
        zero = Sphere(radius = 0.05)
        zero.set_color(ORANGE).shift(axs.c2p(0,0,0))
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( np.cos(v)*np.sin(u), np.sin(v) , np.cos(v)*np.cos(u)  ),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tp = Surface(
            lambda u, v: ax.c2p( u , sqrt(3)*v/2 + 1/2 , sqrt(3)/2 - v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        tpp = Surface(
            lambda u, v: axs.c2p( u , sqrt(3)*v/2 ,- v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        t = ValueTracker(PI/2)
        tp.add_updater( lambda x: x.become( 
            Surface(
                lambda u, v: ax.c2p( np.cos(t.get_value())*u + np.sin(t.get_value())*(sqrt(3)*v/2 + 1/2), np.cos(t.get_value())*(sqrt(3)*v/2 + 1/2) - np.sin(t.get_value())*u , sqrt(3)/2 - v/2 ),
                u_range = [ -1, 1 ],
                v_range = [ -1, 1 ],
                checkerboard_colors = [GREEN, GREEN],
                fill_opacity = 0.7
            )   
        ))
        tpp.add_updater( lambda x: x.become( 
            Surface(
                lambda u, v: axs.c2p( np.cos(t.get_value())*u + np.sin(t.get_value())*(sqrt(3)*v/2), np.cos(t.get_value())*(sqrt(3)*v/2) - np.sin(t.get_value())*u , - v/2 ),
                u_range = [ -1, 1 ],
                v_range = [ -1, 1 ],
                checkerboard_colors = [GREEN, GREEN],
                fill_opacity = 0.7
            )   
        ))
        p.add_updater(lambda x: x.become(Sphere(radius = 0.05).set_color(ORANGE).shift(ax.c2p(np.sin(t.get_value())/2,np.cos(t.get_value())/2,sqrt(3)/2 ))))
        self.add(axs,tpp,zero, ax, sig, p, tp)
        self.wait()
        self.play(t.animate.set_value(PI), rate_func=rate_functions.linear)
        self.wait()


class balance3(ThreeDScene):
    def construct(self):
        ax = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-2/5,6/5,10], x_length = 6, y_length = 6, z_length = 3.2 , tips = False )
        axs = ax.copy()
        ax.shift([(7/4),-(7/4)*sqrt(3),-1])
        axs.shift([-(7/4),(7/4)*sqrt(3),-1])
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0,1/2,sqrt(3)/2 ))
        zero = Sphere(radius = 0.05)
        zero.set_color(ORANGE).shift(axs.c2p(0,0,0))
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( np.cos(v)*np.sin(u), np.sin(v) , np.cos(v)*np.cos(u)  ),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tp = Surface(
            lambda u, v: ax.c2p( u , sqrt(3)*v/2 + 1/2 , sqrt(3)/2 - v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        tpp = Surface(
            lambda u, v: axs.c2p( u , sqrt(3)*v/2 ,- v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        t = ValueTracker(PI)
        tp.add_updater( lambda x: x.become( 
            Surface(
                lambda u, v: ax.c2p( np.cos(t.get_value())*u + np.sin(t.get_value())*(sqrt(3)*v/2 + 1/2), np.cos(t.get_value())*(sqrt(3)*v/2 + 1/2) - np.sin(t.get_value())*u , sqrt(3)/2 - v/2 ),
                u_range = [ -1, 1 ],
                v_range = [ -1, 1 ],
                checkerboard_colors = [GREEN, GREEN],
                fill_opacity = 0.7
            )   
        ))
        tpp.add_updater( lambda x: x.become( 
            Surface(
                lambda u, v: axs.c2p( np.cos(t.get_value())*u + np.sin(t.get_value())*(sqrt(3)*v/2), np.cos(t.get_value())*(sqrt(3)*v/2) - np.sin(t.get_value())*u , - v/2 ),
                u_range = [ -1, 1 ],
                v_range = [ -1, 1 ],
                checkerboard_colors = [GREEN, GREEN],
                fill_opacity = 0.7
            )   
        ))
        p.add_updater(lambda x: x.become(Sphere(radius = 0.05).set_color(ORANGE).shift(ax.c2p(np.sin(t.get_value())/2,np.cos(t.get_value())/2,sqrt(3)/2 ))))
        self.add(axs,tpp,zero, ax, sig, p, tp)
        self.wait()
        self.play(t.animate.set_value(3*PI/2), rate_func=rate_functions.linear)
        self.wait()

class balance4(ThreeDScene):
    def construct(self):
        ax = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-2/5,6/5,10], x_length = 6, y_length = 6, z_length = 3.2 , tips = False )
        axs = ax.copy()
        ax.shift([(7/4),-(7/4)*sqrt(3),-1])
        axs.shift([-(7/4),(7/4)*sqrt(3),-1])
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0,1/2,sqrt(3)/2 ))
        zero = Sphere(radius = 0.05)
        zero.set_color(ORANGE).shift(axs.c2p(0,0,0))
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( np.cos(v)*np.sin(u), np.sin(v) , np.cos(v)*np.cos(u)  ),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tp = Surface(
            lambda u, v: ax.c2p( u , sqrt(3)*v/2 + 1/2 , sqrt(3)/2 - v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        tpp = Surface(
            lambda u, v: axs.c2p( u , sqrt(3)*v/2 ,- v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        t = ValueTracker(3*PI/2)
        tp.add_updater( lambda x: x.become( 
            Surface(
                lambda u, v: ax.c2p( np.cos(t.get_value())*u + np.sin(t.get_value())*(sqrt(3)*v/2 + 1/2), np.cos(t.get_value())*(sqrt(3)*v/2 + 1/2) - np.sin(t.get_value())*u , sqrt(3)/2 - v/2 ),
                u_range = [ -1, 1 ],
                v_range = [ -1, 1 ],
                checkerboard_colors = [GREEN, GREEN],
                fill_opacity = 0.7
            )   
        ))
        tpp.add_updater( lambda x: x.become( 
            Surface(
                lambda u, v: axs.c2p( np.cos(t.get_value())*u + np.sin(t.get_value())*(sqrt(3)*v/2), np.cos(t.get_value())*(sqrt(3)*v/2) - np.sin(t.get_value())*u , - v/2 ),
                u_range = [ -1, 1 ],
                v_range = [ -1, 1 ],
                checkerboard_colors = [GREEN, GREEN],
                fill_opacity = 0.7
            )   
        ))
        p.add_updater(lambda x: x.become(Sphere(radius = 0.05).set_color(ORANGE).shift(ax.c2p(np.sin(t.get_value())/2,np.cos(t.get_value())/2,sqrt(3)/2 ))))
        self.add(axs,tpp,zero, ax, sig, p, tp)
        self.wait()
        self.play(t.animate.set_value(2*PI), rate_func=rate_functions.linear)
        self.wait()





class tcurves(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        thm = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift(0.4*RIGHT + 0.2*DOWN)
        l = thm.get_left()[0]
        thm1 = Tex(r"Let ", r"$\Sigma $ ", r"$ \subset \mathbb{R}^3$ be a surface and ", r"$p$ ", r"$ \in $ ",r"$ \Sigma$.",r" Then", font_size = 34).next_to(thm, DOWN)
        thm1.shift((l-thm1.get_left()[0])*RIGHT)
        thm2 = Tex(r"$ T_p\Sigma$ ", r" $= \{  $ ",r"$\gamma^{\prime}(0) $ ",r"$ \vert $ ",r"$ \gamma$ ",r"$ : [-1,1] \to $ ",r"$ \Sigma $",r" smooth with ",r"$ \gamma (0) = p$ ",r"$ \}$.", font_size = 34).next_to(thm1, DOWN)
        thm2.shift((l-thm2.get_left()[0])*RIGHT)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(thm1[1], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(thm1[5], DOWN).shift(0.15*UP)
        sli3 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(thm2[6], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(thm1[3], DOWN).shift(0.15*UP)
        pli2 = Line([0,0,0], [0.8,0,0], color = ORANGE).next_to(thm2[8], DOWN).shift(0.15*UP)
        tli = Line([0,0,0], [0.6,0,0], color = GREEN).next_to(thm2[0], DOWN).shift(0.15*UP)
        gli = Line([0,0,0], [0.3,0,0], color = MAROON_D).next_to(thm2[4], DOWN).shift(0.15*UP)
        wli = Line([0,0,0], [0.5,0,0], color = RED).next_to(thm2[2], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-1,3,10], x_length = 7, y_length = 7, z_length = 2.8 ).shift([0,0,-1])
        def normal(s,t):
            return 2**(1-(s**2 + t**2))
        sig = Surface(
            lambda u, v: ax.c2p( u, v , normal(u,v) ),
            u_range = [ -3, 3 ],
            v_range = [ -3, 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0,1/4, normal(0,1/4)))
        tps = Surface(
            lambda u, v: ax.c2p( sqrt(2)*u, v + 1/4 , -(2**(-1/16))*math.log(2)*v + normal(0,1/4)),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        g1 = ParametricFunction(lambda t : ax.c2p( t ,  1/4 -t , normal(t,1/4 -t) ), t_range = [-2,2], color = MAROON_D)
        t1 = Arrow(ax.c2p( 0 ,  1/4  , normal(0,1/4 ) ), ax.c2p( 1 ,  -3/4  , normal(0,1/4) + 2**(-1/16)*math.log(2) ), color = RED)
        g2 = ParametricFunction(lambda t : ax.c2p( t , 1/4 +t , normal( t, 1/4 +t) ), t_range = [-2,2], color = MAROON_D)
        t2 = Arrow(ax.c2p( 0 ,  1/4  , normal(0,1/4 ) ), ax.c2p( 1 ,  5/4  , normal(0,1/4) - 2**(-1/16)*math.log(2) ), color = RED)
        g3 = ParametricFunction(lambda t : ax.c2p( -np.sin(t) , 5/4 - np.cos(t) , normal( 5/4 - np.cos(t) , np.sin(t) ) ), t_range = [-PI/2,PI/2], color = MAROON_D)
        t3 = Arrow(ax.c2p( 0 ,  1/4  , normal(0,1/4 ) ), ax.c2p( -2 , 1/4  , normal(0,1/4)  ), color = RED)
        self.add_fixed_in_frame_mobjects(thm, thm1, thm2, pli, tli, gli, pli2,  wli, sli, sli2, sli3)
        self.remove(thm, thm1, thm2, pli, tli, gli, pli2,  wli, sli, sli3, sli2)
        self.play(Create(thm), Create(thm1), Create(thm2), Create(pli), Create(tli),  Create(gli), Create(pli2), Create(wli), Create(sli), Create(sli2), Create(sli3))
        self.wait()
        self.play(Create(ax), Create(sig), Create(p), Create(tps))
        self.wait()
        self.play(Create(g1))
        self.wait()
        self.play(Create(t1))
        self.wait()
        self.play(FadeOut(g1), FadeOut(t1))
        self.wait()
        self.play(Create(g2))
        self.wait()
        self.play(Create(t2))
        self.wait()
        self.play(FadeOut(g2), FadeOut(t2))
        self.wait()
        self.play(Create(g3))
        self.wait()
        self.play(Create(t3))
        self.wait()
        self.play(FadeOut(g3), FadeOut(t3))
        self.wait()
        



class tptc(ThreeDScene):
    def construct(self):
        tan = Tex(r"Proof", font_size = 34, color = BLUE).to_edge(UL).shift(0.4*RIGHT + 0.5*DOWN)
        l = tan.get_left()[0]
        tan1 = Tex(r" For $\phi:$ ",r"$ U $ ", r"$ \to $ ", r"$\Sigma$ ", r"a chart, ", r"$p$ ", r"$ = \phi ($ ", r"$q$ ", r"$)$,  and ", r"$W$ ", r"$= a $", r"$\frac{\partial \phi}{\partial u} (q) $ ", r"$ + b $ ", r"$\frac{\partial \phi}{\partial v} (q) $,", font_size = 34).next_to(tan, DOWN)
        tan1.shift((l-tan1.get_left()[0])*RIGHT)
        alpeq = Tex(r"$t \mapsto $ ", r"$ q + t(a,b) $ ", font_size = 34).next_to(tan1, DOWN)
        alpeq.shift((l-alpeq.get_left()[0])*RIGHT)
        gameq = Tex(r"$\gamma (t) $ ", r"$ : =  \phi ($ ", r"$ q + t(a,b) $ ", r"$)$ satisfies ", r"$\gamma^{\prime } (0) = W$.", font_size = 34).next_to(tan1, DOWN)
        gameq.shift((l-gameq.get_left()[0])*RIGHT)
        uli = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(tan1[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(tan1[3], DOWN).shift(0.15*UP)
        pla = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(tan1[5], DOWN).shift(0.15*UP)
        qla = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(tan1[7], DOWN).shift(0.15*UP)
        wli = Line([0,0,0], [0.3,0,0], color = RED).next_to(tan1[9], DOWN).shift(0.15*UP)
        e3li = Line([0,0,0], [0.5,0,0], color = YELLOW).next_to(tan1[11], DOWN).shift(0.15*UP)
        e4li = Line([0,0,0], [0.5,0,0], color = YELLOW).next_to(tan1[13], DOWN).shift(0.15*UP)
        alpli = Line([0,0,0], [1.5,0,0], color = MAROON_D).next_to(alpeq, DOWN).shift(0.15*UP)
        gamli = Line([0,0,0], [0.6,0,0], color = MAROON_D).next_to(gameq[0], DOWN).shift(0.15*UP)
        gamwli = Line([0,0,0], [1.2,0,0], color = RED).next_to(gameq[4], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-2/5,6/5,10], x_length = 3, y_length = 3, z_length = 1.6 , tips = False )
        ax.shift([-(1.8),(1.8)*sqrt(3),-1.5])
        axu = ThreeDAxes(x_range = [-2,2,10], y_range = [-2,2,10], z_range = [-2,3,10], x_length = 3, y_length = 3, z_length = 0.1 , tips = False )
        axu.shift([(1.8),-sqrt(3)*(1.8),-1.5])
        phi = CurvedArrow([-1.5,-1,0], [1.5,-1,0], radius = -5)
        phil = Tex(r"$\phi$", font_size = 34).next_to(phi, UP)
        q = Sphere(radius = 0.05)
        q.set_color(ORANGE).shift(axu.c2p(0,PI/6,0))
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0,1/2,sqrt(3)/2 ))
        e1 = Arrow( axu.c2p(-0.2,PI/6,0 ), axu.c2p(1.2,PI/6,0), color = YELLOW )
        e2 = Arrow( axu.c2p(0,PI/6-0.2,0 ), axu.c2p(0,1.2+PI/6,0), color = YELLOW )
        e3 = Arrow( ax.c2p(-0.2,1/2,sqrt(3)/2 ), ax.c2p(1.2,1/2,sqrt(3)/2 ), color = YELLOW )
        e4 = Arrow( ax.c2p(0,1/2 - sqrt(3)/10 ,sqrt(3)/2 + 1/10 ), ax.c2p(0,1/2 + sqrt(3)*(1.2)/2,sqrt(3)/2 - 0.6), color = YELLOW )
        e5 = e4.copy()
        w = Arrow( ax.c2p( 0.1,1/2 - sqrt(3)*(1.2)/20 ,sqrt(3)/2 + 0.06 ), ax.c2p(-1,1/2 + sqrt(3)*(1.2)/2  ,sqrt(3)/2 - 0.6  ), color = RED )
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( np.cos(v)*np.sin(u), np.sin(v) , np.cos(v)*np.cos(u)  ),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        u = Surface(
            lambda u, v: axu.c2p( u , v , 0),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        tps = Surface(
            lambda u, v: ax.c2p( u , sqrt(3)*v/2 + 1/2 , sqrt(3)/2 - v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        alp = ParametricFunction(lambda t : axu.c2p( - t ,  PI/6 + t , 0 ) , t_range = [-1,1], color = MAROON_D)
        alp2 = alp.copy()
        gam = ParametricFunction(lambda t : ax.c2p( np.cos(PI/6 + t)*np.sin(-t), np.sin(PI/6 + t) , np.cos(PI/6 + t)*np.cos(-t) ), t_range = [-1,1], color = MAROON_D)
        self.add_fixed_in_frame_mobjects(tan, tan1, phi, phil, uli, sli, pla, qla, wli, alpli, alpeq, gameq, gamli, gamwli, e3li, e4li)
        self.remove(tan, tan1, phi, phil, uli, sli, pla, qla, wli, alpli, alpeq, gameq, gamli, gamwli, e3li, e4li)
        self.play(Create(uli),  Create(sli), Create(tan), Create(tan1), Create(wli), Create(pla), Create(qla), Create(ax), Create(axu), Create(sig), Create(u), Create(phi), Create(phil), Create(q), Create(p), Create(e1), Create(e2), Create(e3), Create(e4), Create(tps), Create(e3li), Create(e4li))
        self.wait()
        self.play(Create(w))
        self.wait()
        self.play(Create(alp), Create(alpeq), Create(alpli))
        self.wait()
        self.add(alp2)
        self.play(Transform(alp2, gam), Create(gameq[0]), Create(gameq[3]), Create(gameq[4]), Create(gameq[1]), Transform(alpeq[1], gameq[2]), FadeOut(alpeq[0]), Transform(alpli, gamli), Create(e5), Create(gamwli))
        self.wait()
        self.play(FadeOut(tan1), FadeOut(alpeq[1]), FadeOut(gameq), FadeOut(e1), FadeOut(e2), FadeOut(e3), FadeOut(e4), FadeOut(e5), FadeOut(alp), FadeOut(alp2), FadeOut(uli), FadeOut(sli), FadeOut(pla), FadeOut(qla), FadeOut(wli), FadeOut(alpli), FadeOut(gamli), FadeOut(gamwli), FadeOut(w), FadeOut(e3li), FadeOut(e4li))
        self.wait()


class tctp(ThreeDScene):
    def construct(self):
        tan = Tex(r"Proof", font_size = 34, color = BLUE).to_edge(UL).shift(0.4*RIGHT + 0.5*DOWN)
        l = tan.get_left()[0]
        tan1 = Tex(r" For $\phi:$ ",r"$ U $ ", r"$ \to $ ", r"$\Sigma$ ", r"a chart, ", r"$p$ ", r"$ = \phi ($ ", r"$q$ ", r"$)$,  and ", r"$\gamma$", r" $: [-1,1] \to $ ", r"$\Sigma $ ", r" smooth with ", r"$\gamma(0) = p$,", font_size = 34).next_to(tan, DOWN)
        tan1.shift((l-tan1.get_left()[0])*RIGHT)
        alpeq = Tex(r"$\alpha (t) $ ", r"$ : =  \phi ^{-1} (\gamma (t)) $.", r" Then ", r"$\gamma ^{\prime } (0)$ ", r"$  = $ ", r"$\frac{\partial \phi}{\partial u} (q)$ ",r"$ \frac{d \alpha_1}{dt} (0) +$ ",r"$ \frac{\partial \phi}{\partial v} (q)$ ",r"$ \frac{d \alpha_2}{dt} (0) $", font_size = 34).next_to(tan1, DOWN)
        alpeq.shift((l-alpeq.get_left()[0])*RIGHT)
        uli = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(tan1[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(tan1[3], DOWN).shift(0.15*UP)
        pla = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(tan1[5], DOWN).shift(0.15*UP)
        qla = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(tan1[7], DOWN).shift(0.15*UP)
        gamli = Line([0,0,0], [0.3,0,0], color = MAROON_D).next_to(tan1[9], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(tan1[11], DOWN).shift(0.15*UP)
        alpli = Line([0,0,0], [0.6,0,0], color = MAROON_D).next_to(alpeq[0], DOWN).shift(0.15*UP)
        pli2 = Line([0,0,0], [1,0,0], color = ORANGE).next_to(tan1[13], DOWN).shift(0.15*UP)
        gamwli = Line([0,0,0], [0.6,0,0], color = RED).next_to(alpeq[3], DOWN).shift(0.15*UP)
        e3li = Line([0,0,0], [0.5,0,0], color = YELLOW).next_to(alpeq[5], DOWN).shift(0.15*UP)
        e4li = Line([0,0,0], [0.5,0,0], color = YELLOW).next_to(alpeq[7], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-2/5,6/5,10], x_length = 3, y_length = 3, z_length = 1.6 , tips = False )
        ax.shift([-(1.8),(1.8)*sqrt(3),-1.5])
        axu = ThreeDAxes(x_range = [-2,2,10], y_range = [-2,2,10], z_range = [-2,3,10], x_length = 3, y_length = 3, z_length = 0.1 , tips = False )
        axu.shift([(1.8),-sqrt(3)*(1.8),-1.5])
        phi = CurvedArrow([-1.5,-1,0], [1.5,-1,0], radius = -5)
        phil = Tex(r"$\phi$", font_size = 34).next_to(phi, UP)
        q = Sphere(radius = 0.05)
        q.set_color(ORANGE).shift(axu.c2p(0,PI/6,0))
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0,1/2,sqrt(3)/2 ))
        #e1 = Arrow( axu.c2p(-0.2,PI/6,0 ), axu.c2p(1.2,PI/6,0), color = YELLOW )
        #e2 = Arrow( axu.c2p(0,PI/6-0.2,0 ), axu.c2p(0,1.2+PI/6,0), color = YELLOW )
        e3 = Arrow( ax.c2p(-0.2,1/2,sqrt(3)/2 ), ax.c2p(1.2,1/2,sqrt(3)/2 ), color = YELLOW )
        e4 = Arrow( ax.c2p(0,1/2 - sqrt(3)/10 ,sqrt(3)/2 + 1/10 ), ax.c2p(0,1/2 + sqrt(3)*(1.2)/2,sqrt(3)/2 - 0.6), color = YELLOW )
        #e5 = e4.copy()
        w = Arrow( ax.c2p( -PI / (20*sqrt(3)),1/2 - sqrt(3)/20 ,sqrt(3)/2 + 1/20 ), ax.c2p( sqrt(3)*PI / 8 ,1/2 + sqrt(3)*3/8 ,sqrt(3)/2  -3/8 ), color = RED )
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( np.cos(v)*np.sin(u), np.sin(v) , np.cos(v)*np.cos(u)  ),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        u = Surface(
            lambda u, v: axu.c2p( u , v , 0),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        tps = Surface(
            lambda u, v: ax.c2p( u , sqrt(3)*v/2 + 1/2 , sqrt(3)/2 - v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        alp = ParametricFunction(lambda t : axu.c2p(  t**2 - PI**2/36  ,  t , 0 ) , t_range = [-1,1], color = MAROON_D)
        gam = ParametricFunction(lambda t : ax.c2p( np.cos(t)*np.sin(t**2- PI**2/36), np.sin( t) , np.cos( t)*np.cos(t**2- PI**2/36) ), t_range = [-1,1], color = MAROON_D)
        gam2 = gam.copy()
        self.add_fixed_in_frame_mobjects(tan, tan1, phi, phil, uli, sli, pla, qla, alpli, alpeq, gamli, gamwli, sli2, pli2, e3li, e4li)
        self.remove(tan1, phi, phil, uli, sli, pla, qla, alpli, alpeq, gamli, gamwli, sli2, pli2, e3li, e4li)
        self.add(ax, axu, u, sig, p, q, tps, phi, phil)
        self.wait()
        self.play(Create(uli),  Create(sli), Create(tan1), Create(pla), Create(qla), Create(gam), Create(gamli), Create(sli2), Create(pli2))
        self.wait()
        self.add(gam2)
        self.play(Write(alpeq[0]), Write(alpeq[1]), Create(alpli), Transform(gam2, alp))
        self.wait()
        self.play(Write(alpeq[2]), Write(alpeq[3]), Create(gamwli), Create(w), Write(alpeq[4]), Create(alpeq[5]), Create(alpeq[6]), Create(alpeq[7]), Create(alpeq[8]))
        self.wait()
        self.play(Create(e3li), Create(e4li), Create(e3), Create(e4))
        self.wait()


class ddd(ThreeDScene):
    def construct(self):
        ddef = Tex(r"Definition", font_size = 34, color = BLUE).to_edge(UL).shift(0.4*RIGHT + 0.2*DOWN)
        l = ddef.get_left()[0]
        dd1 = Tex(r" For $f:$ ",r"$ \Sigma $ ", r"$ \to  \mathbb{R}$ smooth, ", r"$ W $ ", r"$ \in $ ", r"$ T_p \Sigma$, ", r" the ", r"directional derivative ", r" is defined as", font_size = 34).next_to(ddef, DOWN)
        dd1.shift((l-dd1.get_left()[0])*RIGHT)
        dd1[7].set_color(BLUE)
        dd2 = Tex(r"$D_W f  : = (f \circ $ ", r"$ \gamma $ ", r"$ )^{\prime}(0)$, where ", r"$\gamma :$ ", r"$\to $", r"$ \Sigma$ ", r" satisfies ", r"$\gamma^{\prime} (0) = W$.", font_size = 34).next_to(dd1, DOWN)
        dd2.shift((l-dd2.get_left()[0])*RIGHT)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(dd1[1], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(dd2[5], DOWN).shift(0.15*UP)
        gamli = Line([0,0,0], [0.3,0,0], color = MAROON_D).next_to(dd2[1], DOWN).shift(0.15*UP)
        gamli2 = Line([0,0,0], [0.3,0,0], color = MAROON_D).next_to(dd2[3], DOWN).shift(0.15*UP)
        gamwli = Line([0,0,0], [0.9,0,0], color = RED).next_to(dd2[7], DOWN).shift(0.15*UP)
        wli = Line([0,0,0], [0.3,0,0], color = RED).next_to(dd1[3], DOWN).shift(0.15*UP)
        pp0 = Tex(r"Proposition", font_size = 34, color = BLUE).next_to(dd2, DOWN).shift(0.1*DOWN)
        pp0.shift((l-pp0.get_left()[0])*RIGHT)
        pp1 = Tex(r"$D_Wf$ does not depend on ", r"$\gamma$,", r" and $d_pf : $ ", r"$T_p \Sigma$", r" $  \to \mathbb{R}$ is linear.", font_size = 34).next_to(pp0, DOWN)
        pp1.shift((l-pp1.get_left()[0])*RIGHT)
        gamli3 = Line([0,0,0], [0.3,0,0], color = MAROON_D).next_to(pp1[1], DOWN).shift(0.15*UP)
        tpli = Line([0,0,0], [0.3,0,0], color = GREEN).next_to(dd1[5], DOWN).shift(0.15*UP)
        tpli2 = Line([0,0,0], [0.5,0,0], color = GREEN).next_to(pp1[3], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-2/5,6/5,10], x_length = 5, y_length = 5, z_length = 8/3 , tips = False )
        ax.shift([0,0,-2.5])
        p = Sphere(radius = 0.1)
        p.set_color(ORANGE).shift(ax.c2p(0,1/2,sqrt(3)/2 ))
        p2 = p.copy()
        w = Arrow( ax.c2p( 0 ,1/2 ,sqrt(3)/2 ), ax.c2p( 3*sqrt(3)*PI /32 ,1/2 + sqrt(3)*9/32 ,sqrt(3)/2  -9/32 ), color = RED )
        w2 = w.copy()
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( np.cos(v)*np.sin(u), np.sin(v) , np.cos(v)*np.cos(u)  ),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tps = Surface(
            lambda u, v: ax.c2p( u , sqrt(3)*v/2 + 1/2 , sqrt(3)/2 - v/2 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        gam = ParametricFunction(lambda t : ax.c2p( np.cos(t)*np.sin(t**2- PI**2/36), np.sin( t) , np.cos( t)*np.cos(t**2- PI**2/36) ), t_range = [-1,1], color = MAROON_D)
        self.add_fixed_in_frame_mobjects( ddef, dd1, dd2, sli, sli2, gamli, gamli2, gamwli, wli, pp0, pp1, gamli3, tpli, tpli2)
        self.remove(ddef, dd1, dd2, sli, sli2, gamli, gamli2, gamwli, wli, pp0, pp1, gamli3, tpli, tpli2)
        self.play(Create(ddef), Create(dd1), Create(dd2), Create(sli), Create(sli2), Create(gamli), Create(gamli2), Create(gamwli), Create(wli))
        self.wait()
        self.play(Create(ax), Create(sig), Create(w), Create(p))
        self.wait()
        self.play(Create(gam), Create(p2), Create(w2))
        self.wait()
        self.play(Write(pp0), Write(pp1), Create(gamli3), Create(tpli), Create(tpli2), Create(tps))
        self.wait()


class ddpr(ThreeDScene):
    def construct(self):
        pro = Tex(r"Proof", font_size = 34, color = BLUE).to_edge(UL).shift(0.4*RIGHT + 0.5*DOWN)
        l = pro.get_left()[0]
        pro1 = Tex(r" For a chart $\phi:$ ",r"$ U $ ", r"$ \to $ ", r"$\Sigma$,", r" we can write $f \circ \gamma = (f \circ \phi ) \circ ( \phi ^{-1} \circ \gamma )$,", font_size = 34).next_to(pro, DOWN)
        pro1.shift((l-pro1.get_left()[0])*RIGHT)
        pro2 = Tex(r"If ", r"$\alpha $ ", r"$  =  \phi ^{-1} \circ  \gamma $, and ", r"$q $ ", r"$= \phi ^{-1}( $ ", r"$p $ ", r"$)$, ", r" then ", font_size = 34).next_to(pro1, DOWN)
        pro2.shift((l-pro2.get_left()[0])*RIGHT)
        pro3 = Tex(r"$ (f \circ \gamma )^{\prime }(0)   = \frac{\partial (f \circ \phi)}{ \partial u} (q) \frac{d \alpha_1}{dt} (0)  + \frac{\partial (f \circ \phi)}{ \partial v} (q) \frac{d \alpha_2}{dt} (0)  $ ", r"$=  W_1 \frac{\partial (f \circ \phi)}{ \partial u} (q)  + W_2 \frac{\partial (f \circ \phi)}{ \partial v} (q) $", font_size = 34).next_to(pro2, DOWN)
        pro3.shift((l-pro3.get_left()[0])*RIGHT)
        wb = Tex(r"$ W $", r" $ = \gamma^{\prime } (0) = W_1 $ ", r"$\frac{\partial \phi}{\partial u}$", r" $+ W_2$ ", r"$\frac{\partial \phi}{\partial v}$", font_size = 34)
        wb.shift(2.7*RIGHT + 2.7*DOWN)
        uli = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(pro1[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(pro1[3], DOWN).shift(0.15*UP)
        pla = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(pro2[5], DOWN).shift(0.15*UP)
        qla = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(pro2[3], DOWN).shift(0.15*UP)
        ala = Line([0,0,0], [0.3,0,0], color = MAROON_D).next_to(pro2[1], DOWN).shift(0.15*UP)
        wbla = Line([0,0,0], [0.3,0,0], color = RED).next_to(wb[0], DOWN).shift(0.15*UP)
        b1la = Line([0,0,0], [0.3,0,0], color = YELLOW).next_to(wb[2], DOWN).shift(0.15*UP)
        b2la = Line([0,0,0], [0.3,0,0], color = YELLOW).next_to(wb[4], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-2/5,6/5,10], x_length = 3, y_length = 3, z_length = 1.6 , tips = False )
        ax.shift([-(1.8),(1.8)*sqrt(3),-1.5])
        axu = ThreeDAxes(x_range = [-2,2,10], y_range = [-2,2,10], z_range = [-2,3,10], x_length = 3, y_length = 3, z_length = 0.1 , tips = False )
        axu.shift([(1.8),-sqrt(3)*(1.8),-1.5])
        phi = CurvedArrow([-1.5,-1,0], [1.5,-1,0], radius = -5)
        phil = Tex(r"$\phi$", font_size = 34).next_to(phi, UP)
        q = Sphere(radius = 0.05)
        q.set_color(ORANGE).shift(axu.c2p(0,PI/6,0))
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0,1/2,sqrt(3)/2 ))
        w = Arrow( ax.c2p( -PI / (20*sqrt(3)),1/2 - sqrt(3)/20 ,sqrt(3)/2 + 1/20 ), ax.c2p( sqrt(3)*PI / 8 ,1/2 + sqrt(3)*3/8 ,sqrt(3)/2  -3/8 ), color = RED )
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        sig = Surface(
            lambda u, v: ax.c2p( np.cos(v)*np.sin(u), np.sin(v) , np.cos(v)*np.cos(u)  ),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        u = Surface(
            lambda u, v: axu.c2p( u , v , 0),
            u_range = [ -PI/2, PI/2 ],
            v_range = [ -PI/2, PI/2 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.7
        )
        #tps = Surface(
        #    lambda u, v: ax.c2p( u , sqrt(3)*v/2 + 1/2 , sqrt(3)/2 - v/2 ),
        #    u_range = [ -1, 1 ],
        #    v_range = [ -1, 1 ],
        #    checkerboard_colors = [GREEN, GREEN],
        #    fill_opacity = 0.7
        #)
        e3 = Arrow( ax.c2p(-0.2,1/2,sqrt(3)/2 ), ax.c2p(1.2,1/2,sqrt(3)/2 ), color = YELLOW )
        e4 = Arrow( ax.c2p(0,1/2 - sqrt(3)/10 ,sqrt(3)/2 + 1/10 ), ax.c2p(0,1/2 + sqrt(3)*(1.2)/2,sqrt(3)/2 - 0.6), color = YELLOW )
        alp = ParametricFunction(lambda t : axu.c2p(  t**2 - PI**2/36  ,  t , 0 ) , t_range = [-1,1], color = MAROON_D)
        gam = ParametricFunction(lambda t : ax.c2p( np.cos(t)*np.sin(t**2- PI**2/36), np.sin( t) , np.cos( t)*np.cos(t**2- PI**2/36) ), t_range = [-1,1], color = MAROON_D)
        gam2 = gam.copy()
        self.add_fixed_in_frame_mobjects(pro, pro1, pro2, phi, phil, uli, sli, pla, qla, ala, pro3, wb, wbla, b1la, b2la)
        self.remove(pro, pro1, pro2, phi, phil, uli, sli, pla, qla, ala, pro3, wb, wbla, b1la, b2la)
        self.play(Create(uli),  Create(sli), Create(pro), Create(pro1))
        self.wait()
        self.play(Create(ax), Create(axu), Create(u), Create(sig), Create(p), Create(q), Create(gam), Create(w), Create(phi), Create(phil))
        self.wait()
        self.play(Transform(gam2, alp), Write(pro2), Write(pro3[0]), Create(pla), Create(qla), Create(ala))
        self.wait()
        self.play(Write(pro3[1]), Write(e3), Write(e4), Write(wb), Create(wbla), Create(b1la), Create(b2la))
        self.wait()


class dlin(ThreeDScene):
    def construct(self):
        obs1 = Tex(r"If $f:$ ", r"$\Sigma_1$", r" $\to$ ", r"$\Sigma _2$", r" is smooth and $f($ ", r"$p$", r" $) =$ ", r"$q$", r", then $dpf :$ ", r"$T_p\Sigma _1$", r" $\to $ ", r"$T_{q}\Sigma _2$", r" is linear.", font_size = 34).to_edge(UL).shift(0.4*RIGHT + 1.2*DOWN)
        s1li = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(obs1[1], DOWN).shift(0.15*UP)
        s2li = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(obs1[3], DOWN).shift(0.15*UP)
        pla = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(obs1[5], DOWN).shift(0.15*UP)
        qla = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(obs1[7], DOWN).shift(0.15*UP)
        tpli = Line([0,0,0], [0.6,0,0], color = GREEN).next_to(obs1[9], DOWN).shift(0.15*UP)
        tqli = Line([0,0,0], [0.6,0,0], color = GREEN).next_to(obs1[11], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,3,10], x_length = 5, y_length = 5, z_length = 2.5 , tips = False )
        ax2 = ax.copy()
        ax.shift([(1.8),-(1.8)*sqrt(3),-1])
        ax2.shift([-(1.8),(1.8)*sqrt(3),-1])
        f = CurvedArrow([-1,0,0], [1,0,0], radius = -6)
        fl = Tex(r"$f$", font_size = 34).next_to(f, UP)
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0,1,2**(3/4) ))
        q = Sphere(radius = 0.05)
        q.set_color(ORANGE).shift(ax2.c2p(0,1,1/2 ))
        s1 = Surface(
            lambda u, v: ax.c2p( u , v , 2*2**(- (u**2 + v**2)/4) ),
            u_range = [-4,4],
            v_range = [-4,4],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        s11 = s1.copy()
        s2 = Surface(
            lambda u, v: ax2.c2p( u , v ,  (np.cos(PI*u/2) +  np.cos(PI*v/2))/2 ),
            u_range = [-4,4],
            v_range = [-4,4],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tps = Surface(
            lambda u, v: ax.c2p( u , v +1, 2**(3/4) - v*math.log(2)*2**(-1/4) ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        tpss = tps.copy()
        tqs = Surface(
            lambda u, v: ax2.c2p( u , v+1 , 1/2 -PI*v/4 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        self.add_fixed_in_frame_mobjects(obs1, s1li, s2li, pla, qla, tpli, tqli, f , fl)
        self.remove(obs1, s1li, s2li, pla, qla, tpli, tqli, f , fl)
        self.play(Create(obs1), Create(s1li), Create(s2li), Create(pla), Create(qla), Create(tpli), Create(tqli), Create(f), Create(fl), Create(ax), Create(ax2), Create(s1), Create(s2), Create(p), Create(q))
        self.wait()
        self.add(s11)
        self.play(Transform(s11, s2))
        self.remove(s2)
        self.wait()
        self.play(Create(tps), Create(tqs))
        self.wait()
        self.add(tpss)
        self.play(Transform(tpss, tqs))
        self.remove(tqs)
        self.wait()



class jaco(ThreeDScene):
    def construct(self):
        ti = Tex(r"Jacobian of functions", font_size = 34, color = BLUE).to_edge(UL).shift(0.4*RIGHT + 0.4*DOWN )
        l = ti.get_left()[0]
        jac0 = Tex(r"If $f:$ ", r"$\Sigma$", r" $\to \mathbb{R}^k$ is smooth, its Jacobian at ", r"$p$", r" $ \in \Sigma$ is defined as", font_size = 34).next_to(ti, DOWN)
        jac0.shift((l - jac0.get_left()[0])*RIGHT)
        jac1 = Tex(r"Jac$(f)(p) : = $Jac$(d_pf)$.", font_size = 34).next_to(jac0, DOWN)
        jac1.shift(jac1.get_center()[0]*LEFT)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(jac0[1], DOWN).shift(0.15*UP)
        pla = Line([0,0,0], [0.2,0,0], color = ORANGE).next_to(jac0[3], DOWN).shift(0.15*UP)
        ax = ThreeDAxes(x_range = [-5,5,10], y_range = [-5,5,10], z_range = [-2,3,10], x_length = 5, y_length = 5, z_length = 2.5 , tips = False )
        ax2 = ax.copy()
        ax.shift([(1.8),-(1.8)*sqrt(3),-1])
        ax2.shift([-(1.8),(1.8)*sqrt(3),-1])
        f = CurvedArrow([-1,0,0], [1,0,0], radius = -6)
        fl = Tex(r"$f$", font_size = 34).next_to(f, UP)
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift(ax.c2p(0,1,2**(3/4) ))
        p2 = p.copy()
        q = Sphere(radius = 0.05)
        q.set_color(ORANGE).shift(ax2.c2p(0,1,1/2 ))
        s1 = Surface(
            lambda u, v: ax.c2p( u , v , 2*2**(- (u**2 + v**2)/4) ),
            u_range = [-4,4],
            v_range = [-4,4],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        s11 = s1.copy()
        s2 = Surface(
            lambda u, v: ax2.c2p( u , v ,  (np.cos(PI*u/2) +  np.cos(PI*v/2))/2 ),
            u_range = [-4,4],
            v_range = [-4,4],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.7
        )
        tps = Surface(
            lambda u, v: ax.c2p( u , v +1, 2**(3/4) - v*math.log(2)*2**(-1/4) ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        tpss = tps.copy()
        tqs = Surface(
            lambda u, v: ax2.c2p( u , v+1 , 1/2 -PI*v/4 ),
            u_range = [ -1, 1 ],
            v_range = [ -1, 1 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        self.add_fixed_in_frame_mobjects(jac0, jac1, sli, pla, f , fl, ti)
        self.remove(jac0, jac1, pla, sli, f , fl, ti)
        self.play(Create(ti), Create(jac0), Create(jac1), Create(sli), Create(pla), Create(f), Create(fl), Create(ax), Create(ax2), Create(s1),  Create(p))
        self.wait()
        self.add(s11, p2)
        self.play(Transform(s11, s2), Transform(p2,q))
        self.wait()
        self.play(Create(tps))
        self.wait()
        self.add(tpss)
        self.play(Transform(tpss, tqs))
        self.wait()


class area(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        axes = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-1,3,10], x_length = 4, y_length = 4, z_length = 8/3 )
        dom = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-1,1,2], x_length = 4, y_length = 4, z_length = 0.01, tips = False )
        axes.shift([-sqrt(3), 3, -1.5])
        dom.shift([sqrt(3), -3,-1.5])
        u = Surface(
            lambda u, v: dom.c2p( u,v,0 ),
            u_range = [-2.5, 2.5],
            v_range = [-2.5,2.5],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        )
        sig = Surface(
            lambda u, v: axes.c2p( u,v, 1 + np.sin(2*u)/3 + np.sin(2*v)/3 ),
            u_range = [-2.5, 2.5],
            v_range = [-2.5,2.5],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7
        )
        f = CurvedArrow([-1,-0.5,0], [1,-0.5,0], radius = -4.5)
        fl = Tex(r"$\phi$", font_size = 34).next_to(f, UP).shift(0.15*DOWN)
        cap = Tex(r"Area", font_size = 36, color = BLUE).to_edge(UL).shift(RIGHT + (0.4)* DOWN)
        l = cap.get_left()[0]
        cap1 = Tex(r"If a chart $\phi :$ ", r"$ U$", r" $ \to$ ", r"$ \Sigma$", r" covers ", r"$\Sigma$", r", the area is defined as:", font_size = 36).next_to(cap, DOWN)
        cap1.shift((l - cap1.get_left()[0])*RIGHT)
        cap2 = Tex(r"Area$(\Sigma ) = \int_U$ Jac$(\phi )$", r" $ = \int_U \vert \frac{\partial \phi}{ \partial u} \times \frac{\partial \phi }{\partial v} \vert dudv $", font_size = 36).next_to(cap1, DOWN)
        cap2.shift((cap2.get_center()[0])*LEFT)
        uli = Line([0,0,0], [0.3,0,0], color = GREEN).next_to(cap1[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(cap1[3], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(cap1[5], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(cap, cap1, cap2, f, fl, uli, sli, sli2)
        self.remove(cap, cap1, cap2, f, fl, uli, sli, sli2)
        self.play(Write(cap), Write(cap1), Write(cap2[0]), Create(axes), Create(dom), Create(sig), Create(u), Create(f), Create(fl), Create(uli), Create(sli), Create(sli2))  
        self.wait()
        self.play(Write(cap2[1]))
        self.wait()


class int0(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        axes = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-1,3,10], x_length = 4, y_length = 4, z_length = 8/3 )
        dom = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-1,1,2], x_length = 4, y_length = 4, z_length = 0.01, tips = False )
        axes.shift([-sqrt(3), 3, -1.5])
        dom.shift([sqrt(3), -3,-1.5])
        u = Surface(
            lambda u, v: dom.c2p( u,v,0 ),
            u_range = [-2.5, 2.5],
            v_range = [-2.5,2.5],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.7
        ) 
        rd = Surface(
            lambda u, v: dom.c2p( u,v,0 ),
            u_range = [-2,-1],
            v_range = [1,2],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7,
            resolution = 5
        )
        rdd = rd.copy()
        r = Surface(
            lambda u, v: axes.c2p( u,v, 1 + np.sin(2*u)/3 + np.sin(2*v)/3 ),
            u_range = [-2,-1],
            v_range = [1,2],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7,
            resolution = 5
        )
        sig = Surface(
            lambda u, v: axes.c2p( u,v, 1 + np.sin(2*u)/3 + np.sin(2*v)/3 ),
            u_range = [-2.5, 2.5],
            v_range = [-2.5,2.5],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7
        )
        f = CurvedArrow([-1,-0.5,0], [1,-0.5,0], radius = -4.5)
        fl = Tex(r"$\phi$", font_size = 34).next_to(f, UP).shift(0.15*DOWN)
        cap = Tex(r"Integrals on surfaces", font_size = 36, color = BLUE).to_edge(UL).shift(RIGHT + (0.4)* DOWN)
        l = cap.get_left()[0]
        cap1 = Tex(r"If a chart $\phi :$ ", r"$ U$", r" $ \to$ ", r"$ \Sigma$", r" covers ", r"$ R $", r" $\subset $ ", r"$\Sigma$", r", and $f: $ ", r"$\Sigma$", r" $  \to \mathbb{R}$ is continuous,", font_size = 36).next_to(cap, DOWN)
        cap1.shift((l - cap1.get_left()[0])*RIGHT)
        cap2 = Tex(r"the integral of $f$ over $R$ is defined as:", font_size = 36).next_to(cap1, DOWN)
        cap2.shift((l - cap2.get_left()[0])*RIGHT)
        formula = Tex(r"$\int_R f : = \int_{\phi^{-1}(R)} (f \circ \phi  ) \vert \frac{\partial \phi}{ \partial u} \times \frac{\partial \phi }{\partial v} \vert dudv $.", font_size = 36).next_to(cap2, DOWN)
        formula.shift((formula.get_center()[0])*LEFT)
        uli = Line([0,0,0], [0.3,0,0], color = GREEN).next_to(cap1[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(cap1[3], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(cap1[7], DOWN).shift(0.15*UP)
        sli3 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(cap1[9], DOWN).shift(0.15*UP)
        rli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(cap1[5], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(cap, cap1, cap2, f, fl, uli, sli, sli2, rli, sli3, formula)
        self.remove(cap1, cap2, f, fl, uli, sli, sli2, rli, sli3, formula)
        self.add( f, fl, dom, axes, u, sig)
        self.wait()
        self.play(Write(cap1), Create(uli), Create(sli), Create(sli2), Create(rli), Create(sli3), Create(r), Create(rd))  
        self.wait()
        self.add(rdd)
        self.play(Transform(rdd,r))
        self.remove(r)
        self.wait()
        self.play(Write(cap2), Write(formula))
        self.wait()

class int1(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        axes = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-1,3,10], x_length = 5, y_length = 5, z_length = 10/3 )
        axes.shift([0,0, -1.5])
        r = Surface(
            lambda u, v: axes.c2p( u,v, 2**( -( u**2 + v**2 )/8 ) ),
            u_range = [0,1],
            v_range = [0,1],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7,
            resolution = 5
        )
        sig1 = Surface(
            lambda u, v: axes.c2p( u,v, 2**( -( u**2 + v**2 )/4 ) ),
            u_range = [0,2.5],
            v_range = [-2.5,2.5],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7,
            resolution = 5
        )
        cap = Tex(r"Integrals on surfaces", font_size = 36, color = BLUE).to_edge(UL).shift(RIGHT + (0.4)* DOWN)
        l = cap.get_left()[0]
        cap1 = Tex(r"If ", r"$ R $", r" $\subset $ ", r"$\Sigma$", r" is not covered by a single chart, and $f: $ ", r"$\Sigma$", r" $ \to \mathbb{R}$ is continuous,", font_size = 36).next_to(cap, DOWN)
        cap1.shift((l - cap1.get_left()[0])*RIGHT)
        cap2 = Tex(r"the integral of $f$ over ", r"$R$", r" is defined as:", font_size = 36).next_to(cap1, DOWN)
        cap2.shift((l - cap2.get_left()[0])*RIGHT)
        formula = Tex(r"$\int_R f : = \sum_i \int_{R_i} f $", font_size = 36).next_to(cap2, DOWN)
        formula.shift((formula.get_center()[0])*LEFT)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(cap1[3], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(cap1[5], DOWN).shift(0.15*UP)
        rli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(cap1[1], DOWN).shift(0.15*UP)
        rli2 = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(cap2[1], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(cap, cap1, cap2, sli, sli2, rli2, formula, rli)
        self.remove(cap1, cap2, sli, sli2, rli, rli2, formula, rli)
        self.wait()
        self.play(Write(cap1), Create(sli), Create(sli2), Create(rli), Create(r))  
        self.wait()
        self.play(Write(cap2), Write(formula), Create(rli2))
        self.wait()



class int2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi= 68 * DEGREES, theta= 30*DEGREES)
        axes = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-1,3,10], x_length = 5, y_length = 5, z_length = 10/3 )
        axes.shift([0,0,-1])
        sig1 = Surface(
            lambda u, v: axes.c2p(((1-2*np.sin(PI*u/4)**2)/2 + 2*np.sin(PI*u/4)**2 *np.sin(v))* np.sin(v) , u, ((1-2*np.sin(PI*u/4)**2)/2 + 2*np.sin(PI*u/4)**2 *np.sin(v))* np.cos(v)),
            u_range = [0,1],
            v_range = [0,PI],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7,
            resolution = 20
        )
        sig2 = Surface(
            lambda u, v: axes.c2p(-((1-2*np.sin(PI*u/4)**2)/2 + 2*np.sin(PI*u/4)**2 *np.sin(v))* np.sin(v) , u, ((1-2*np.sin(PI*u/4)**2)/2 + 2*np.sin(PI*u/4)**2 *np.sin(v))* np.cos(v)),
            u_range = [0,1],
            v_range = [0,PI],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7,
            resolution = 20
        )
        sig3 = Surface(
            lambda u, v: axes.c2p(((1-2*np.sin(PI*u/4)**2)/2 + 2*np.sin(PI*u/4)**2 *np.sin(v))* np.sin(v) , -u, ((1-2*np.sin(PI*u/4)**2)/2 + 2*np.sin(PI*u/4)**2 *np.sin(v))* np.cos(v)),
            u_range = [0,1],
            v_range = [0,PI],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7,
            resolution = 20
        )
        sig4 = Surface(
            lambda u, v: axes.c2p(-((1-2*np.sin(PI*u/4)**2)/2 + 2*np.sin(PI*u/4)**2 *np.sin(v))* np.sin(v) , - u, ((1-2*np.sin(PI*u/4)**2)/2 + 2*np.sin(PI*u/4)**2 *np.sin(v))* np.cos(v)),
            u_range = [0,1],
            v_range = [0,PI],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7,
            resolution = 20
        )
        s1 = Surface(
            lambda u, v: axes.c2p( np.sin(v)**2 + u /2  , u + 1 , np.sin(v) * np.cos(v)),
            u_range = [0,1],
            v_range = [0,PI],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7,
            resolution = 20
        )
        s2 = Surface(
            lambda u, v: axes.c2p( - np.sin(v)**2 - u /2 , u + 1 , np.sin(v) * np.cos(v)),
            u_range = [0,1],
            v_range = [0,PI],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7,
            resolution = 20
        )
        s3 = Surface(
            lambda u, v: axes.c2p( np.sin(v)**2 + u /2  , - u - 1 , np.sin(v) * np.cos(v)),
            u_range = [0,1],
            v_range = [0,PI],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7,
            resolution = 20
        )
        s4 = Surface(
            lambda u, v: axes.c2p( - np.sin(v)**2 - u /2 , - u - 1 , np.sin(v) * np.cos(v)),
            u_range = [0,1],
            v_range = [0,PI],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7,
            resolution = 20
        )
        t1 = Surface(
            lambda u, v: axes.c2p( (1 + ( np.cos(v) )/2)*np.cos(u) , 2 +(1 + ( np.cos(v) )/2)*np.sin(u) , np.sin(v)/2  ),
            u_range = [0,PI/2],
            v_range = [0,TAU],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7,
            resolution = 20
        )
        t2 = Surface(
            lambda u, v: axes.c2p( (1 + ( np.cos(v) )/2)*np.cos(u) , - 2 - (1 + ( np.cos(v) )/2)*np.sin(u) , np.sin(v)/2  ),
            u_range = [0,PI/2],
            v_range = [0,TAU],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7,
            resolution = 20
        )
        t3 = Surface(
            lambda u, v: axes.c2p( (1 + ( np.cos(v) )/2)*np.cos(u) , 2 +(1 + ( np.cos(v) )/2)*np.sin(u) , np.sin(v)/2  ),
            u_range = [PI/2,PI],
            v_range = [0,TAU],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.7,
            resolution = 20
        )
        t4 = Surface(
            lambda u, v: axes.c2p( (1 + ( np.cos(v) )/2)*np.cos(u) , - 2 - (1 + ( np.cos(v) )/2)*np.sin(u) , np.sin(v)/2  ),
            u_range = [PI/2,PI],
            v_range = [0,TAU],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7,
            resolution = 20
        )
        cap = Tex(r"Integrals on surfaces", font_size = 36, color = BLUE).to_edge(UL).shift(RIGHT + (0.4)* DOWN)
        l = cap.get_left()[0]
        cap1 = Tex(r"If ", r"$ R $", r" $\subset $ ", r"$\Sigma$", r" is not covered by a single chart, and $f: $ ", r"$\Sigma$", r" $ \to \mathbb{R}$ is continuous,", font_size = 36).next_to(cap, DOWN)
        cap1.shift((l - cap1.get_left()[0])*RIGHT)
        cap2 = Tex(r"the integral of $f$ over ", r"$R$", r" is defined as:", font_size = 36).next_to(cap1, DOWN)
        cap2.shift((l - cap2.get_left()[0])*RIGHT)
        formula = Tex(r"$\int_R f : = \sum_i \int_{R_i} f $", font_size = 36).next_to(cap2, DOWN)
        formula.shift((formula.get_center()[0])*LEFT)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(cap1[3], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(cap1[5], DOWN).shift(0.15*UP)
        rli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(cap1[1], DOWN).shift(0.15*UP)
        rli2 = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(cap2[1], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(cap, cap1, cap2, sli, sli2, rli2, formula, rli)
        self.remove(cap1, cap2, sli, sli2, rli, rli2, formula)
        self.wait()
        self.play(Create(cap1), Create(sli), Create(sli2), Create(rli), Create(sig1), Create(sig2), Create(sig3), Create(sig4), Create(s1), Create(s2), Create(s3), Create(s4), Create(t1), Create(t2), Create(t3), Create(t4))
        self.wait()
        self.play(Create(cap2), Create(formula), Create(rli2))
        self.wait()



class totest(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        axes = ThreeDAxes(x_range = [-3,3,10], y_range = [-3,3,10], z_range = [-1,3,10], x_length = 5, y_length = 5, z_length = 10/3 )
        t1 = Surface(
            lambda u, v: axes.c2p( (1 + ( np.cos(v) )/2)*np.cos(u) , 1.5 +(1 + ( np.cos(v) )/2)*np.sin(u) , np.sin(v)/2  ),
            u_range = [0,PI],
            v_range = [0,TAU],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.7,
            resolution = 10
        )
        self.wait()
        self.play(Create(t1))
        self.wait()

class btest(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        ax = ThreeDAxes(x_range = [-1.5,1.5,10], y_range = [-1.5,1.5,10], z_range = [-2/5,6/5,10], x_length = 6, y_length = 6, z_length = 3.2 , tips = False )
        p = Dot()
        self.play(Create(p))
        self.wait()
        t = ValueTracker(0)
        p.add_updater( lambda x: x.become( Dot([0,0, t.get_value()]) ))
        self.play(t.animate.set_value(3), rate_func=rate_functions.ease_in_sine)
        self.wait()



class logtest(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)
        ax = ThreeDAxes(x_range = [-5,5,1], y_range = [-1.5,1.5,10], z_range = [-2/5,6/5,10], x_length = 6, y_length = 6, z_length = 3.2 , tips = False )
        p1 = Dot( ax.c2p(math.log(2),0,0), color = YELLOW)
        self.play(Create(p1), Create(ax))
        self.wait()


class dete(Scene):
    def construct(self):
        ti = Tex(r"Jacobian", font_size = 34, color = BLUE).to_edge(UL).shift(0.6*RIGHT + 0.5*DOWN)
        tl = Tex(r"$T =  \begin{pmatrix} 1 & 1/2 \\  1/2 & 1 \end{pmatrix} $", font_size = 34).shift(1.9*UP)
        det = Tex(r"Jac$ (T)  = $ ", r"$ \dfrac{Area(T(A))}{Area(A)}$", r" $ = $ ", r"$\dfrac{Area(T(B))}{Area(B)}$", r" $ = \vert $det$(T) \vert =  \dfrac{3}{4} $", font_size = 34).shift(2*DOWN)
        ula = Line([0,0,0], [1.5,0,0], color = PURPLE_B).next_to(det[1], DOWN).shift(0.15*UP)
        ulb = Line([0,0,0], [1.5,0,0], color = ORANGE).next_to(det[3], DOWN).shift(0.15*UP)
        t = CurvedArrow([-1.5,0.9,0], [2,0.7,0], radius = -6)
        ax = NumberPlane()
        ax1 = ax.copy()
        a = Circle(color = PURPLE_B, fill_opacity = 0.5, radius = 0.3).shift((0.5)*LEFT+(0.5)*UP)
        a1 = a.copy()
        a.shift(4*LEFT+(0.5)*UP)
        a1.apply_function( lambda p: np.array([p[0] + p[1]/2 , p[0]/2 + p[1], 0, ] ) )
        a1.shift(4*RIGHT + (0.7)* UP)
        a0 = a.copy()
        b = Triangle(color=ORANGE, fill_opacity=0.5)
        b.scale(0.35)
        b.shift((0.5)*RIGHT+(0.2)*DOWN)
        b1 = b.copy()
        b.shift(4*LEFT+(0.5)*UP)
        b1.apply_function( lambda p: np.array([p[0] + p[1]/2 , p[0]/2 + p[1], 0, ] ) )
        b1.shift(4*RIGHT +(0.7)* UP)
        b0 = b.copy()
        ax1.apply_function( lambda p: np.array([p[0] + p[1]/2 , p[0]/2 + p[1], 0, ] ) )
        ax.scale(0.25)
        ax.shift(4*LEFT+(0.5)*UP)
        ax1.scale(0.25)
        ax1.shift(4*RIGHT+(0.7)*UP)
        ax0 = ax.copy()
        self.play(Create(ax), Create(t), Create(tl), Create(ti))
        self.wait()
        self.add(ax0)
        self.play(Transform(ax0, ax1))
        self.wait()
        self.play(Create(a), Create(b))
        self.wait()
        self.add(a0, b0)
        self.play(Transform(a0, a1), Transform(b0,b1))
        self.wait()
        self.play(Write(det), Create(ula), Create(ulb))
        self.wait()
        