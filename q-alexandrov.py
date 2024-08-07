from asyncio import threads
from this import d
from tkinter import E
from manim import *
from numpy import sqrt
import math



class ti(Scene):
    def construct(self):
        t1 = Text("Comparison geometry", font_size=60)
        self.play(Write(t1))
        self.wait()        



class gbf3(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/4)
        plab = Tex(r"$K > 0$", font_size = 30).shift(4.5*LEFT + 0.9*UP )
        zlab = Tex(r"$K = 0$", font_size = 30).shift( 0.9 * UP )
        nlab = Tex(r"$K < 0$", font_size = 30).shift(4.5*RIGHT + 0.9*UP )
        boxp = SurroundingRectangle(plab, buff = .1, color = RED)
        boxz = SurroundingRectangle(zlab, buff = .1, color = PURPLE_B)
        boxn = SurroundingRectangle(nlab, buff = .1, color = BLUE)
        r = 1.2
        h = -1.1
        ht = 3.2
        m = 1.6
        mn = 1.3
        op = 0.2
        q = 0.18
        def f(u,v):
            return [r*np.cos(u)*np.sin(v) + ht , r*np.sin(u)*np.sin(v) - ht , r*np.cos(v) + h - 0.15  ]
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
            v_range = [ 0 , PI ],
            checkerboard_colors = [RED, RED],
            fill_opacity = op*0.75, 
            resolution = [ 32, 32 ]
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
        self.add_fixed_in_frame_mobjects(plab, zlab, nlab, boxp, boxz, boxn)
        self.remove(plab, zlab, nlab, boxp, boxz, boxn)
        self.wait()
        self.play(
            Create(plab), Create(zlab), Create(nlab), 
            Create(boxp), Create(boxz), Create(boxn), 
            Create(sp), Create(deltap), Create(gpa), Create(gpb), Create(gpc), 
            Create(sz), Create(deltaz), Create(gza), Create(gzb), Create(gzc), 
            Create(sn), Create(deltan), Create(gna), Create(gnb), Create(gnc)
        )
        self.wait()





class ctd(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=78 * DEGREES, theta=  PI/4)
        de = Tex(r"Definition", font_size = 34, color = BLUE).to_edge(UL).shift(0.5*RIGHT + 0.5*DOWN)
        de1 = Tex(r"A ", r"triangle", r" $[x,y,z]$ in a surface is a collection of three points $x$, $y$, $z$ $\in \Sigma$,", font_size = 34).next_to(de, DOWN)
        de2 = Tex(r"with a choice of minimizing geodesics ", r"$[x,y]$", r", ", r"$[x,z]$", r", ", r"$[y,z]$", r".", font_size = 34).next_to(de1, DOWN)
        de3 = Tex(r"A ", r"comparison triangle", r" is a triangle $[\overline{x}, \overline{y}, \overline{z} ]$ in $\mathbb{R}^2$ with", font_size = 34).next_to(de2, DOWN)
        de4 = Tex(r"$d(\overline{x},\overline{y}) = d(x,y)$ ", r"$d(\overline{x},\overline{z}) = d(x,z)$ ", r"$d(\overline{y},\overline{z}) = d(y, z)$ ", font_size = 34).next_to(de3, DOWN)
        l = de.get_left()[0]
        de1.shift((l - de1.get_left()[0])*RIGHT)
        de2.shift((l - de2.get_left()[0])*RIGHT)
        de3.shift((l - de3.get_left()[0])*RIGHT)
        de4.shift(  de4.get_center()[0]  *LEFT)
        de1[1].set_color(BLUE)
        de3[1].set_color(BLUE)
        de4[0].shift(0.5 * LEFT)
        de4[2].shift(0.5 * RIGHT)
        xyli = Line( [0,0,0] , [0.6,0,0], color = ORANGE).next_to(de2[1], DOWN).shift(0.15*UP)
        xzli = Line( [0,0,0] , [0.6,0,0], color = PURPLE).next_to(de2[3], DOWN).shift(0.15*UP)
        yzli = Line( [0,0,0] , [0.6,0,0], color = GREEN ).next_to(de2[5], DOWN).shift(0.15*UP)
        xy2li = Line( [0,0,0] , [2,0,0], color = ORANGE ).next_to(de4[0], DOWN).shift(0.15*UP)
        xz2li = Line( [0,0,0] , [2,0,0], color = PURPLE ).next_to(de4[1], DOWN).shift(0.15*UP)
        yz2li = Line( [0,0,0] , [2,0,0], color = GREEN  ).next_to(de4[2], DOWN).shift(0.15*UP)
        r = 1.5
        m = 2.2
        a = 2
        h = - 1.8
        s = Surface(
            lambda u, v: [ r * np.cos(u)*np.sin(v) + m , r * np.sin(u) * np.sin(v) - m , r * np.cos(v) + h ],
            u_range = [ 0 , TAU ],
            v_range = [ 0 ,  PI ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 32, 24]
        )
        x  = Sphere(radius = 0.05).set_color(BLUE).shift([   r + m , - m ,  h ])
        y  = Sphere(radius = 0.05).set_color(BLUE).shift([   m , r - m ,  h   ])
        z  = Sphere(radius = 0.05).set_color(BLUE).shift([   m ,  sqrt(3) * r /2 - m , r /2  + h  ])
        xy = ParametricFunction(lambda t : [  r * np.cos(t)   + m  ,  r * np.sin(t)   - m ,   h  ] , t_range = [0, PI/2 ], color = ORANGE)
        xz = ParametricFunction(lambda t : [  r * np.cos(t)   + m  ,  r * np.sin(t) *sqrt(3) / 2  - m ,  r * np.sin(t) / 2  + h ] , t_range = [0, PI/2 ], color = PURPLE)
        yz = ParametricFunction(lambda t : [  m  ,  r * np.cos(t)  - m ,  r * np.sin(t)  + h  ] , t_range = [ 0 , PI / 6 ], color = GREEN )
        r2 = Surface(
            lambda u, v: [ u - m , v + m , h ],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 16, 16]
        )
        xb = Sphere(radius = 0.05).set_color(BLUE).shift([   1.5 - m  ,        m , h ])
        yb = Sphere(radius = 0.05).set_color(BLUE).shift([  -1.5 - m  ,  0.5 + m , h ])
        zb = Sphere(radius = 0.05).set_color(BLUE).shift([  -1.5 - m  , -0.8 + m , h ])
        xbyb = Line( xb.get_center() , yb.get_center() , color = ORANGE)
        xbzb = Line( xb.get_center() , zb.get_center() , color = PURPLE)
        ybzb = Line( yb.get_center() , zb.get_center() , color = GREEN )
        self.add_fixed_in_frame_mobjects(de, de1, de2, de3, de4, xyli, xzli, yzli, xy2li, xz2li, yz2li)
        self.remove(de, de1, de2, de3, de4, xyli, xzli, yzli, xy2li, xz2li, yz2li)
        self.play(Create(de), Create(de1), Create(de2), Create(s), Create(x), Create(y), Create(z), Create(xy), Create(xz), Create(yz), Create(xyli), Create(xzli), Create(yzli))
        self.wait()
        self.play(Create(de3), Create(de4), Create(xy2li), Create(xz2li), Create(yz2li), Create(r2), Create(xb), Create(yb), Create(zb), Create(xbyb), Create(xbzb), Create(ybzb))
        self.wait()


class cad(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=78 * DEGREES, theta=  PI/4)
        de = Tex(r"Definition", font_size = 34, color = BLUE).to_edge(UL).shift(0.5*RIGHT + 0.5*DOWN)
        de1 = Tex(r"A ", r"triangle", r" $[x,y,z]$ in a surface is a collection of three points $x$, $y$, $z$ $\in \Sigma$,", font_size = 34).next_to(de, DOWN)
        de2 = Tex(r"with a choice of minimizing geodesics ", r"$[x,y]$", r", ", r"$[x,z]$", r", ", r"$[y,z]$", r".", font_size = 34).next_to(de1, DOWN)
        de3 = Tex(r"A ", r"comparison triangle", r" is a triangle $[\overline{x}, \overline{y}, \overline{z} ]$ in $\mathbb{R}^2$ with", font_size = 34).next_to(de2, DOWN)
        de4 = Tex(r"$d(\overline{x},\overline{y}) = d(x,y)$ ", r"$d(\overline{x},\overline{z}) = d(x,z)$ ", r"$d(\overline{y},\overline{z}) = d(y, z)$ ", font_size = 34).next_to(de3, DOWN)
        de5 = Tex(r"The ", r"comparison angle ", r"$\tilde{\measuredangle}(x _y ^z ) $", r" is defined to be ", r"$\measuredangle ( \overline{x}_{\overline{y} }^{\overline{z}} )$", r" $ \in [0, \pi ]$.", font_size = 34).shift(2*UP)
        de6 = Tex(r"In general, one has ", r"$\measuredangle (x_y^z)$", r" $ \neq $ ", r"$ \tilde{\measuredangle}(x_y^z)$.", font_size = 34).next_to(de5, DOWN)
        l = de.get_left()[0]
        de1.shift((l - de1.get_left()[0])*RIGHT )
        de2.shift((l - de2.get_left()[0])*RIGHT )
        de3.shift((l - de3.get_left()[0])*RIGHT )
        de4.shift(   de4.get_center()[0] * LEFT )
        de6.shift(   de6.get_center()[0] * LEFT )
        de1[1].set_color(BLUE)
        de3[1].set_color(BLUE)
        de5[1].set_color(BLUE)
        de4[0].shift(0.5 * LEFT)
        de4[2].shift(0.5 * RIGHT)
        xyli = Line( [0,0,0] , [0.6,0,0], color = ORANGE).next_to(de2[1], DOWN).shift(0.15*UP)
        xzli = Line( [0,0,0] , [0.6,0,0], color = PURPLE).next_to(de2[3], DOWN).shift(0.15*UP)
        yzli = Line( [0,0,0] , [0.6,0,0], color = GREEN ).next_to(de2[5], DOWN).shift(0.15*UP)
        xy2li = Line( [0,0,0] , [2,0,0], color = ORANGE ).next_to(de4[0], DOWN).shift(0.15*UP)
        xz2li = Line( [0,0,0] , [2,0,0], color = PURPLE ).next_to(de4[1], DOWN).shift(0.15*UP)
        yz2li = Line( [0,0,0] , [2,0,0], color = GREEN  ).next_to(de4[2], DOWN).shift(0.15*UP)
        cali  = Line( [0,0,0] , [1,0,0], color = YELLOW ).next_to(de5[2], DOWN).shift(0.15*UP)
        cadli = Line( [0,0,0] , [1,0,0], color = YELLOW ).next_to(de5[4], DOWN).shift(0.15*UP)
        ca2li = Line( [0,0,0] , [1,0,0], color = YELLOW ).next_to(de6[3], DOWN).shift(0.15*UP)
        aali  = Line( [0,0,0] , [1,0,0], color = TEAL   ).next_to(de6[1], DOWN).shift(0.15*UP)
        r = 1.5
        m = 2.2
        a = 2
        h = - 1.8
        rs = 0.25
        th = 0.3
        s = Surface(
            lambda u, v: [ r * np.cos(u)*np.sin(v) + m , r * np.sin(u) * np.sin(v) - m , r * np.cos(v) + h ],
            u_range = [ 0 , TAU ],
            v_range = [ 0 ,  PI ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 32, 24]
        )
        x  = Sphere(radius = 0.05).set_color(BLUE).shift([   r + m , - m ,  h ])
        y  = Sphere(radius = 0.05).set_color(BLUE).shift([   m , r - m ,  h   ])
        z  = Sphere(radius = 0.05).set_color(BLUE).shift([   m ,  sqrt(3) * r /2 - m , r /2  + h  ])
        xy = ParametricFunction(lambda t : [  r * np.cos(t)   + m  ,  r * np.sin(t)   - m ,   h  ] , t_range = [0, PI/2 ], color = ORANGE)
        xz = ParametricFunction(lambda t : [  r * np.cos(t)   + m  ,  r * np.sin(t) *sqrt(3) / 2  - m ,  r * np.sin(t) / 2  + h ] , t_range = [0, PI/2 ], color = PURPLE)
        yz = ParametricFunction(lambda t : [  m  ,  r * np.cos(t)  - m ,  r * np.sin(t)  + h  ] , t_range = [ 0 , PI / 6 ], color = GREEN )
        r2 = Surface(
            lambda u, v: [ u - m , v + m , h ],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 16, 16]
        )
        xb = Sphere(radius = 0.05).set_color(BLUE).shift([   1.5 - m  ,        m , h ])
        yb = Sphere(radius = 0.05).set_color(BLUE).shift([  -1.5 - m  ,  0.5 + m , h ])
        zb = Sphere(radius = 0.05).set_color(BLUE).shift([  -1.5 - m  , -0.8 + m , h ])
        xbyb = Line( xb.get_center() , yb.get_center() , color = ORANGE)
        xbzb = Line( xb.get_center() , zb.get_center() , color = PURPLE)
        ybzb = Line( yb.get_center() , zb.get_center() , color = GREEN )
        ca = ParametricFunction(lambda t : [  1.5 - m  - ( 3 * np.cos(t) + 0.5 * np.sin(t) ) * rs , m  + ( 0.5 * np.cos(t) - 3 * np.sin(t)) * rs , h  ] , t_range = [0, PI / 7 ], color = YELLOW)
        aa = ParametricFunction(lambda t : [  r * np.cos(th) + m  , r * np.sin(th) * np.cos(t) - m , r * np.sin(th) * np.sin(t) + h  ] , t_range = [0, PI / 6 ], color = TEAL)
        old = VGroup(de, de1, de2, de3, de4, xyli, xzli, yzli, xy2li, xz2li, yz2li)
        self.add_fixed_in_frame_mobjects(de, de1, de2, de3, de4, xyli, xzli, yzli, xy2li, xz2li, yz2li, de5, de6, cali, cadli, ca2li, aali)
        self.remove(de5, de6, cali, cadli, ca2li, aali)
        self.add(s, x, y, z, xy, xz, yz, r2, xb, yb, zb, xbyb, xbzb, ybzb)
        self.wait()
        self.play(FadeOut(old), Create(de5), Create(cali), Create(cadli), Create(ca))
        self.wait()
        self.play(Create(de6), Create(ca2li), Create(aali), Create(aa))
        self.wait()






class monot(Scene):
    def construct(self):
        le = Tex(r"Monotonicity Lemma", font_size = 34, color = BLUE).to_edge(UL).shift(0.5*RIGHT + 0.5*DOWN)
        le1 = Tex(r"If $a,b,c, a^{\prime}, b^{\prime}, c^{\prime} \in \mathbb{R}^3$ are such that", font_size = 34).next_to(le, DOWN)
        le2 = Tex(r"$d(a ^{\prime} , b^{\prime} ) = d (a,b) $, ", r" $d(a ^{\prime} , c^{\prime}) = d(a,c)$,", font_size = 34).next_to(le1, DOWN)
        le3 = Tex(r"then ", font_size = 34).next_to(le2, DOWN)
        le4 = Tex(r"$\measuredangle (a _b ^c)$", r" $ \leq$ ", r"$ \measuredangle (a^{\prime c^{\prime}}_{b^{\prime}})$", r" if and only if ", r"$d(b, c)$", r" $ \leq$ ", r"$ d(b^{\prime}, c^{\prime })$", font_size = 34).next_to(le2, DOWN)
        l = le.get_left()[0]
        le1.shift((l - le1.get_left()[0])*RIGHT )
        le2.shift( le2.get_center()[0]   *LEFT  )
        le2[0].shift(0.5 * LEFT )
        le2[1].shift(0.5 * RIGHT)
        le3.shift((l - le3.get_left()[0])*RIGHT )
        le4.shift( le4.get_center()[0]   *LEFT  )
        abli = Line( [0,0,0] , [2,0,0], color = ORANGE).next_to(le2[0], DOWN).shift(0.15*UP)
        acli = Line( [0,0,0] , [2,0,0], color = PURPLE).next_to(le2[1], DOWN).shift(0.15*UP)
        aali = Line( [0,0,0] , [0.8,0,0], color = YELLOW).next_to(le4[0], DOWN).shift(0.15*UP)
        apli = Line( [0,0,0] , [0.8,0,0], color = TEAL  ).next_to(le4[2], DOWN).shift(0.15*UP)
        bcli = Line( [0,0,0] , [0.8,0,0], color = GREEN ).next_to(le4[4], DOWN).shift(0.15*UP)
        bcpp = Line( [0,0,0] , [0.8,0,0], color = PINK  ).next_to(le4[6], DOWN).shift(0.15*UP)
        r = 1.5
        w = 2.3
        h = - 2.2
        th = - PI / 8
        a = 1.5
        c2 = a * r
        a  = Dot([0,0,0]  , color = BLUE)
        b  = Dot([r, 0 ,0], color = BLUE)
        c  = Dot([0, c2 ,0], color = BLUE)
        ap = Dot([0, 0 ,0], color = BLUE)
        bp = Dot([r*np.cos(th), r*np.sin(th), 0], color = BLUE)
        cp = Dot([ c2 * np.cos(PI/2 - th), c2 * np.sin(PI/2 - th), 0], color = BLUE)
        ab = Line(a.get_center(), b.get_center(), color = ORANGE)
        apbp = Line(ap.get_center(), bp.get_center(), color = ORANGE)
        ac = Line(a.get_center(), c.get_center(), color = PURPLE)
        apcp = Line(ap.get_center(), cp.get_center(), color = PURPLE)
        bc = Line(b.get_center(), c.get_center(), color = GREEN)
        bpcp = Line(bp.get_center(), cp.get_center(), color = PINK)
        caa = Angle(ab, ac, color = YELLOW, radius = 0.2)
        aap = Angle(apbp, apcp, color = TEAL, radius = 0.2)
        al = Tex(r"$a$", font_size = 34).next_to(a, DOWN).shift(0.2*LEFT  + 0.2*UP )
        bl = Tex(r"$b$", font_size = 34).next_to(b, DOWN).shift(0.2*RIGHT + 0.2*UP )
        cl = Tex(r"$c$", font_size = 34).next_to(c, DOWN).shift(0.2*UP + 0.2 * LEFT)
        apl = Tex(r"$a^{\prime}$", font_size = 34).next_to(ap, DOWN).shift(0.2*LEFT + 0.2*UP)
        bpl = Tex(r"$b^{\prime}$", font_size = 34).next_to(bp, DOWN).shift(0.2*RIGHT + 0.2 *UP)
        cpl = Tex(r"$c^{\prime}$", font_size = 34).next_to(cp, DOWN).shift(0.2*UP + 0.2 * LEFT)
        abc = VGroup(a, b, c, ab, ac, bc, caa, al, bl, cl)
        abcp = VGroup(ap, bp, cp, apbp, apcp, bpcp, aap, apl, bpl, cpl)
        abc.shift(w*LEFT + 1.1 * h * UP)
        abcp.shift(w*RIGHT + h *UP)
        self.play(Create(abc), Create(abcp), Create(le), Create(le1), Create(le2), Create(le3), Create(abli), Create(acli), Create(aali), Create(apli), Create(bcli), Create(bcpp), Create(le4))
        self.wait()




class chc(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/4)
        thm  = Tex(r"Theorem", font_size = 34, color = BLUE).to_edge(UL).shift(0.5*RIGHT + 0.2*DOWN)
        th1 = Tex(r"Let $\Sigma$ be a complete surface, and $[ x y z ]$ a triangle in $\Sigma$.", font_size = 34).next_to(thm, DOWN)
        th2 = Tex(r"$\bullet$ If $\Sigma$ has non-positive curvature and is simply connected, then", font_size = 34).next_to(th1, DOWN)
        th3 = Tex(r"$\measuredangle (x_y^z)$", r" $\leq $ ", r"$\tilde{\measuredangle}(x_y^z)$", font_size = 34).next_to(th2, DOWN)
        th4 = Tex(r"$ \bullet$ If $\Sigma$ has non-negative curvature, then", font_size = 34).next_to(th3, DOWN)
        th5 = Tex(r"$\measuredangle (x_y^z)$", r" $\geq $ ", r"$\tilde{\measuredangle}(x_y^z)$", font_size = 34).next_to(th4, DOWN)
        l = thm.get_left()[0]
        th1.shift((l - th1.get_left()[0])*RIGHT)
        th2.shift((l - th2.get_left()[0])*RIGHT)
        th4.shift((l - th4.get_left()[0])*RIGHT)
        th3.shift(   th3.get_center()[0] *LEFT )
        th5.shift(   th5.get_center()[0] *LEFT )
        aali  = Line( [0,0,0] , [1,0,0], color = TEAL   ).next_to(th3[0], DOWN).shift(0.15*UP)
        cali  = Line( [0,0,0] , [1,0,0], color = YELLOW ).next_to(th3[2], DOWN).shift(0.15*UP)
        aa2li = Line( [0,0,0] , [1,0,0], color = TEAL   ).next_to(th5[0], DOWN).shift(0.15*UP)
        ca2li = Line( [0,0,0] , [1,0,0], color = YELLOW ).next_to(th5[2], DOWN).shift(0.15*UP)
        r = 1.5
        m = 2.2
        a = 2
        h = - 2
        rs = 0.25
        th = 0.3
        s = Surface(
            lambda u, v: [ r * np.cos(u)*np.sin(v) + m , r * np.sin(u) * np.sin(v) - m , r * np.cos(v) + h ],
            u_range = [ 0 , TAU ],
            v_range = [ 0 ,  PI ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 32, 24]
        )
        x  = Sphere(radius = 0.05).set_color(BLUE).shift([   r + m , - m ,  h ])
        y  = Sphere(radius = 0.05).set_color(BLUE).shift([   m , r - m ,  h   ])
        z  = Sphere(radius = 0.05).set_color(BLUE).shift([   m ,  sqrt(3) * r /2 - m , r /2  + h  ])
        xy = ParametricFunction(lambda t : [  r * np.cos(t)   + m  ,  r * np.sin(t)   - m ,   h  ] , t_range = [0, PI/2 ], color = ORANGE)
        xz = ParametricFunction(lambda t : [  r * np.cos(t)   + m  ,  r * np.sin(t) *sqrt(3) / 2  - m ,  r * np.sin(t) / 2  + h ] , t_range = [0, PI/2 ], color = PURPLE)
        yz = ParametricFunction(lambda t : [  m  ,  r * np.cos(t)  - m ,  r * np.sin(t)  + h  ] , t_range = [ 0 , PI / 6 ], color = GREEN )
        r2 = Surface(
            lambda u, v: [ u - m , v + m , h ],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 16, 16]
        )
        xb = Sphere(radius = 0.05).set_color(BLUE).shift([   1.5 - m  ,        m , h ])
        yb = Sphere(radius = 0.05).set_color(BLUE).shift([  -1.5 - m  ,  0.5 + m , h ])
        zb = Sphere(radius = 0.05).set_color(BLUE).shift([  -1.5 - m  , -0.8 + m , h ])
        xbyb = Line( xb.get_center() , yb.get_center() , color = ORANGE)
        xbzb = Line( xb.get_center() , zb.get_center() , color = PURPLE)
        ybzb = Line( yb.get_center() , zb.get_center() , color = GREEN )
        ca = ParametricFunction(lambda t : [  1.5 - m  - ( 3 * np.cos(t) + 0.5 * np.sin(t) ) * rs , m  + ( 0.5 * np.cos(t) - 3 * np.sin(t)) * rs , h  ] , t_range = [0, PI / 7 ], color = YELLOW)
        aa = ParametricFunction(lambda t : [  r * np.cos(th) + m  , r * np.sin(th) * np.cos(t) - m , r * np.sin(th) * np.sin(t) + h  ] , t_range = [0, PI / 6 ], color = TEAL)
        zz = 10
        hh = -1
        q = 5
        sn = Surface(
            lambda u, v: [ u + m , v - m , ( u ** 2 - v ** 2 ) / zz + hh ],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 32, 24]
        )
        xn  = Sphere(radius = 0.05).set_color(BLUE).shift([        m ,     - m ,            hh ])
        yn  = Sphere(radius = 0.05).set_color(BLUE).shift([- 1.5 + m ,     - m ,  2.25/zz + hh ])
        zn  = Sphere(radius = 0.05).set_color(BLUE).shift([        m , 1.5 - m , -2.25/zz + hh  ])
        xyn = ParametricFunction(lambda t : [                                     - t + m ,                                - m ,                                                             t**2 / zz + hh  ] , t_range = [0, 1.5 ], color = ORANGE)
        xzn = ParametricFunction(lambda t : [                                           m ,                              t - m ,                                                            -t**2 / zz + hh  ] , t_range = [0, 1.5 ], color = PURPLE)
        yzn = ParametricFunction(lambda t : [ - 1.5 + t + ( 9/16 - ( t - 3/4)**2 )/ q + m , t - ( 9/16 - ( t - 3/4)**2 ) / q - m , ( 2.25 - 3 * t + ( 4 * t - 3 ) * ( 9/16 - ( t - 3/4)**2 ) / q )/ zz + hh  ] , t_range = [0, 1.5 ], color = GREEN )
        r22 = Surface(
            lambda u, v: [ u - m , v + m , hh ],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 16, 16]
        )
        leng = 1.7
        thet = 0.3
        xbn = Sphere(radius = 0.05).set_color(BLUE).shift([                       - m ,                         m , hh ])
        ybn = Sphere(radius = 0.05).set_color(BLUE).shift([ - leng * np.cos(thet) - m , - leng * np.sin(thet) + m , hh ])
        zbn = Sphere(radius = 0.05).set_color(BLUE).shift([   leng * np.sin(thet) - m ,   leng * np.cos(thet) + m , hh ])
        rad = 0.25
        can = ParametricFunction(lambda t : [ - m + rad * np.cos(t) ,   m  + rad * np.sin(t), hh  ] , t_range = [ PI/ 2 - thet , PI + thet ], color = YELLOW)
        aan = ParametricFunction(lambda t : [   m + rad * np.cos(t) , - m  + rad * np.sin(t), hh  ] , t_range = [ PI/ 2        , PI        ], color = TEAL  )
        dotsn = VGroup(xbn, ybn, zbn, can)
        dotsn.shift([0.5,-0.5,0])
        xbybn = Line( xbn.get_center() , ybn.get_center() , color = ORANGE)
        xbzbn = Line( xbn.get_center() , zbn.get_center() , color = PURPLE)
        ybzbn = Line( ybn.get_center() , zbn.get_center() , color = GREEN )
        nfig = VGroup(sn, xn, yn, zn, xyn, xzn, yzn, r22, xbn, ybn, zbn, xbybn, xbzbn, ybzbn, can, aan)
        pfig = VGroup(s, x, y, z, xy, xz, yz, r2, xb, yb, zb, xbyb, xbzb, ybzb, ca, aa)
        self.add_fixed_in_frame_mobjects(thm, th1, th2, th3, th4, th5, cali, aali, ca2li, aa2li)
        self.remove(thm, th1, th2, th3, th4, th5, cali, aali, ca2li, aa2li)
        knt = VGroup(thm, th1, th2, th3, cali, aali)
        kpt = VGroup(th4, th5, ca2li, aa2li)
        self.play(Create(knt), Create(nfig))
        self.wait()
        self.play(FadeOut(nfig))
        self.wait()
        self.play(Create(kpt), Create(pfig))
        self.wait()




class pi1(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/4)
        bad = Tex(r"$\measuredangle (x_y^z)$", r" $ = \pi > \pi / 3 =$ ", r"$\tilde{\measuredangle } (x_y^z)$", font_size = 34).shift(2*UP)
        mali = Line([0,0,0], [1,0,0], color = TEAL  ).next_to(bad[0], DOWN).shift(0.15*UP)
        tmli = Line([0,0,0], [1,0,0], color = YELLOW).next_to(bad[2], DOWN).shift(0.15*UP)
        r = 1
        m = 2.2
        a = 2
        rr = r * 2 * PI / sqrt(27)
        rad = 0.25
        h = -1.2
        s = Surface(
            lambda u, v: [ r * np.cos(u)*np.cosh(v/2) + m , r * np.sin(u) * np.cosh(v/2) - m , np.sinh(v/1.5) + h ],
            u_range = [ 0 , TAU ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 32, 24]
        )
        x  = Sphere(radius = 0.05).set_color(BLUE).shift([    r + m ,                  - m , h ])
        y  = Sphere(radius = 0.05).set_color(BLUE).shift([ -r/2 + m ,  sqrt(3) * r / 2 - m , h ])
        z  = Sphere(radius = 0.05).set_color(BLUE).shift([ -r/2 + m , -sqrt(3) * r / 2 - m , h ])
        xy = ParametricFunction(lambda t : [ r * np.cos(t) + m  ,  r * np.sin(t)   - m ,          h          ] , t_range = [         0 , 2 * PI / 3 ], color = ORANGE)
        xz = ParametricFunction(lambda t : [ r * np.cos(t) + m  ,  r * np.sin(t)   - m ,          h          ] , t_range = [4 * PI / 3 ,        TAU ], color = PURPLE)
        yz = ParametricFunction(lambda t : [ r * np.cos(t) + m  ,  r * np.sin(t)   - m ,          h          ] , t_range = [2 * PI / 3 , 4 * PI / 3 ], color = GREEN )
        aa = ParametricFunction(lambda t : [ r             + m  ,  rad * np.cos(t) - m , rad * np.sin(t) + h ] , t_range = [         0 ,         PI ], color = TEAL  )
        r2 = Surface(
            lambda u, v: [ u - m , v + m , h ],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 16, 16]
        )
        xb = Sphere(radius = 0.05).set_color(BLUE).shift([  rr   - m  ,                 m , h ])
        yb = Sphere(radius = 0.05).set_color(BLUE).shift([ -rr/2 - m  ,  sqrt(3)*rr/2 + m , h ])
        zb = Sphere(radius = 0.05).set_color(BLUE).shift([ -rr/2 - m  , -sqrt(3)*rr/2 + m , h ])
        xbyb = Line( xb.get_center() , yb.get_center() , color = ORANGE)
        xbzb = Line( xb.get_center() , zb.get_center() , color = PURPLE)
        ybzb = Line( yb.get_center() , zb.get_center() , color = GREEN )
        ca = ParametricFunction(lambda t : [ rr - m - rad * np.cos(t) , m + rad * np.sin(t), h ] , t_range = [- PI / 6, PI / 6 ], color = YELLOW)
        fig = VGroup(s, x, y, z, xy, xz, yz, aa, r2, xb, yb, zb, xbyb, xbzb, ybzb, ca)
        self.add_fixed_in_frame_mobjects(bad, mali, tmli)
        self.remove(bad, mali, tmli)
        self.play(Create(bad), Create(mali), Create(tmli), Create(fig))
        self.wait()



class knegalt(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/4)
        ex = Tex(r"Cartan--Hadamard Theorem", font_size = 34, color = BLUE).shift(  UP + 3.3 *LEFT )
        ex1 = Tex(r"If $\Sigma \subset \mathbb{R}^3$ is a complete surface of non-negative Gauss curvature,", font_size = 34).next_to(ex, DOWN)
        ex2 = Tex(r"then for all $p \in \Sigma$, the map $\text{exp}_p : T_p \Sigma \to \Sigma $ is non-singular,", font_size = 34).next_to(ex1, DOWN)
        ex3 = Tex(r"non-contracting, and a covering map. ", font_size = 34).next_to(ex2, DOWN)
        thm  = Tex(r"Theorem", font_size = 34, color = BLUE).to_edge(UL).shift(0.5*RIGHT + 0.2 *DOWN)
        th1 = Tex(r"Let $\Sigma$ be a complete surface, and $[ x y z ]$ a triangle in $\Sigma$.", font_size = 34).next_to(thm, DOWN)
        th2 = Tex(r"If $\Sigma$ has non-positive curvature and is simply connected, then", font_size = 34).next_to(th1, DOWN)
        th3 = Tex(r"$\measuredangle (x_y^z)$", r" $\leq $ ", r"$\tilde{\measuredangle}(x_y^z)$", font_size = 34).next_to(th2, DOWN)
        pr  = Tex(r"Proof", font_size = 34, color = BLUE).next_to(th3, DOWN)
        pr1 = Tex(r"Since $\Sigma $ is simply connected, $\exp_x : T_x \Sigma \to \Sigma$ is a homeomorphism.", font_size = 34).next_to(pr, DOWN)
        pr2 = Tex(r"Let $v ,w \in T_x \Sigma$ with ", r"$\exp_x(v) = y$", r", ", r"$\exp_x (w) = z$", r", ", r"$\gamma$", r" $ : [0,1] \to \Sigma $",  font_size = 34).next_to(pr1, DOWN)
        pr3 = Tex(r"a minimizing curve between $y$ and $z$, and ", r"$\eta$", r" $: [0,1] \to T_x \Sigma $ its lift.", font_size = 34).next_to(pr2, DOWN)
        l = ex.get_left()[0]
        ex1.shift((l - ex1.get_left()[0])*RIGHT)
        ex2.shift((l - ex2.get_left()[0])*RIGHT)
        ex3.shift((l - ex3.get_left()[0])*RIGHT)
        thm.shift((l - thm.get_left()[0])*RIGHT)
        th1.shift((l - th1.get_left()[0])*RIGHT)
        th2.shift((l - th2.get_left()[0])*RIGHT)
        th3.shift(   th3.get_center()[0] * LEFT)
        pr.shift((l - pr.get_left()[0])*RIGHT)
        pr1.shift((l - pr1.get_left()[0])*RIGHT)
        pr2.shift((l - pr2.get_left()[0])*RIGHT)
        pr3.shift((l - pr3.get_left()[0])*RIGHT)
        vyli = Line([0,0,0], [1.8,0,0], color = BLUE ).next_to(pr2[1], DOWN).shift(0.15*UP)
        wzli = Line([0,0,0], [1.8,0,0], color = BLUE ).next_to(pr2[3], DOWN).shift(0.15*UP)
        gli  = Line([0,0,0], [0.4,0,0], color = PINK ).next_to(pr2[5], DOWN).shift(0.15*UP)
        eli  = Line([0,0,0], [0.4,0,0], color = GREEN).next_to(pr3[1], DOWN).shift(0.15*UP)
        tool = VGroup(ex , ex1, ex2, ex3)
        thms = VGroup(thm, th1, th2, th3)
        m = 2.4
        a = 1.9
        h = - 2.5
        zz = 10
        q = 5
        s = Surface(
            lambda u, v: [ u - m , v + m , ( u ** 2 - v ** 2 ) / zz + h ],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 32, 24]
        )
        x  = Sphere(radius = 0.05).set_color(GOLD).shift([      - m ,     + m ,            h ])
        y  = Sphere(radius = 0.05).set_color(BLUE).shift([- 1.5 - m ,     + m ,  2.25/zz + h ])
        z  = Sphere(radius = 0.05).set_color(BLUE).shift([      - m , 1.5 + m , -2.25/zz + h ])
        xy = ParametricFunction(lambda t : [                                     - t - m ,                                    m ,                                                           t**2 / zz + h  ] , t_range = [0, 1.5 ], color = ORANGE)
        xz = ParametricFunction(lambda t : [                                         - m ,                                t + m ,                                                          -t**2 / zz + h  ] , t_range = [0, 1.5 ], color = PURPLE)
        yz = ParametricFunction(lambda t : [ - 1.5 + t + ( 9/16 - ( t - 3/4)**2 )/ q - m , t - ( 9/16 - ( t - 3/4)**2 ) / q + m , ( 2.25 - 3 * t + ( 4 * t - 3 ) * ( 9/16 - ( t - 3/4)**2 ) / q )/ zz + h  ] , t_range = [0, 1.5 ], color = PINK  )
        txs = Surface(
            lambda u, v: [ u + m , v - m , h ],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 16, 16]
        )
        eta = ParametricFunction(lambda t : [ - 1.7 + t + ( 289/400 - ( t - 0.85)**2 )/ q + m , t - ( 289/400 - ( t - 0.85)**2 ) / q - m , h ], t_range = [0, 1.7], color = GREEN)
        zer = Sphere(radius = 0.05).set_color(GOLD).shift([         m ,     - m , h ])
        v   = Sphere(radius = 0.05).set_color(BLUE).shift([ - 1.7 + m ,     - m , h ])
        w   = Sphere(radius = 0.05).set_color(BLUE).shift([         m , 1.7 - m , h ])
        ar = Arrow([-0.9,-2.3,0], [0.9,-2.3,0], max_tip_length_to_length_ratio = 0.1 , stroke_width = 3)
        arla = Tex(r"$\exp_x$", font_size = 34).next_to(ar, UP).shift(0.15*DOWN)
        fig = VGroup(s, x, y, z, txs, zer, ar, arla)
        self.add_fixed_in_frame_mobjects(tool, thms, vyli, wzli, gli, eli, pr, pr1, pr2, pr3, ar, arla)
        self.remove(tool, thms, vyli, wzli, gli, eli, pr, pr1, pr2, pr3, ar, arla)
        self.play(Create(tool))
        self.wait()
        self.play(FadeOut(tool))
        self.wait()
        self.play(Create(thms))
        self.wait()
        self.play(Create(pr), Create(pr1), Create(fig))
        self.wait()
        pr22 = VGroup(pr2[0], pr2[1], pr2[2], pr2[3], pr2[4], vyli, wzli)
        self.play(Create(pr22), Create(v), Create(w))
        self.wait()
        self.play(Create(pr2[5]), Create(pr2[6]), Create(pr3), Create(eli), Create(gli), Create(yz), Create(eta))
        self.wait()
        old = VGroup(thm, th1, th2, th3, pr, pr1)
        new = VGroup(pr2, pr3, gli, eli, vyli, wzli)
        self.play(FadeOut(old), new.animate.shift( (thm.get_center()[1] - pr2.get_center()[1] - 0.5 )*UP))
        self.wait()


class kneg2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/4)
        ex = Tex(r"Corollary", font_size = 34, color = BLUE).shift( 5 *LEFT )
        thm  = Tex(r"Theorem", font_size = 34, color = BLUE).to_edge(UL).shift(0.5*RIGHT + 0.2 *DOWN)
        pr2 = Tex(r"Let $v ,w \in T_x \Sigma$ with ", r"$\exp_x(v) = y$", r", ", r"$\exp_x (w) = z$", r", ", r"$\gamma$", r" $ : [0,1] \to \Sigma $",  font_size = 34).next_to(thm, RIGHT).shift(0.5*DOWN)
        pr3 = Tex(r"a minimizing curve between $y$ and $z$, and ", r"$\eta$", r" $: [0,1] \to T_x \Sigma $ its lift.", font_size = 34).next_to(pr2, DOWN)
        pr4 = Tex(r"Since $\exp_x$ is non-contracting,", font_size = 34).next_to(pr3, DOWN)
        pr5 = Tex(r"$d (v,w) \leq $ ", r"length$(\eta )$", r" $ \leq $ ", r"length$(\gamma)$", r" $ = d(y,z)$.", font_size = 34).next_to(pr4, DOWN)
        pr6 = Tex(r"Also, ", r"$\vert v \vert = d(x,y)$", r", ", r"$\vert w  \vert = d(x,z)$", r" and ", r"$\measuredangle (0_v^w) = \measuredangle (x_y^z)$", r",", font_size = 34).next_to(pr5, DOWN)
        pr7 = Tex(r"so by monotonicity, $\measuredangle (x_y^z) \leq \tilde{\measuredangle}(x_y^z)$", font_size = 34).next_to(pr6, DOWN).shift(0.2*DOWN)
        l = ex.get_left()[0]
        pr2.shift((l - pr2.get_left()[0])*RIGHT)
        pr3.shift((l - pr3.get_left()[0])*RIGHT)
        pr4.shift((l - pr4.get_left()[0])*RIGHT)
        pr5.shift(   pr5.get_center()[0] *LEFT )
        pr6.shift((l - pr6.get_left()[0])*RIGHT)
        pr7.shift(   pr7.get_center()[0] *LEFT )
        vyli = Line([0,0,0], [1.8,0,0], color = BLUE  ).next_to(pr2[1], DOWN).shift(0.15*UP)
        wzli = Line([0,0,0], [1.8,0,0], color = BLUE  ).next_to(pr2[3], DOWN).shift(0.15*UP)
        gli  = Line([0,0,0], [0.4,0,0], color = PINK  ).next_to(pr2[5], DOWN).shift(0.15*UP)
        eli  = Line([0,0,0], [0.4,0,0], color = GREEN ).next_to(pr3[1], DOWN).shift(0.15*UP)
        g2li = Line([0,0,0], [1,0,0],   color = PINK  ).next_to(pr5[3], DOWN).shift(0.15*UP)
        e2li = Line([0,0,0], [1,0,0],   color = GREEN ).next_to(pr5[1], DOWN).shift(0.15*UP)
        vyll = Line([0,0,0], [1.5,0,0], color = ORANGE).next_to(pr6[1], DOWN).shift(0.15*UP)
        wzll = Line([0,0,0], [1.5,0,0], color = PURPLE).next_to(pr6[3], DOWN).shift(0.15*UP)
        call = Line([0,0,0], [1.7,0,0], color = YELLOW).next_to(pr6[5], DOWN).shift(0.15*UP)
        m = 2.4
        a = 1.9
        h = - 2.5
        zz = 10
        q = 5
        s = Surface(
            lambda u, v: [ u - m , v + m , ( u ** 2 - v ** 2 ) / zz + h ],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 32, 24]
        )
        x  = Sphere(radius = 0.05).set_color(GOLD).shift([      - m ,     + m ,            h ])
        y  = Sphere(radius = 0.05).set_color(BLUE).shift([- 1.5 - m ,     + m ,  2.25/zz + h ])
        z  = Sphere(radius = 0.05).set_color(BLUE).shift([      - m , 1.5 + m , -2.25/zz + h ])
        xy = ParametricFunction(lambda t : [                                     - t - m ,                                    m ,                                                           t**2 / zz + h  ] , t_range = [0, 1.5 ], color = ORANGE)
        xz = ParametricFunction(lambda t : [                                         - m ,                                t + m ,                                                          -t**2 / zz + h  ] , t_range = [0, 1.5 ], color = PURPLE)
        yz = ParametricFunction(lambda t : [ - 1.5 + t + ( 9/16 - ( t - 3/4)**2 )/ q - m , t - ( 9/16 - ( t - 3/4)**2 ) / q + m , ( 2.25 - 3 * t + ( 4 * t - 3 ) * ( 9/16 - ( t - 3/4)**2 ) / q )/ zz + h  ] , t_range = [0, 1.5 ], color = PINK  )
        txs = Surface(
            lambda u, v: [ u + m , v - m , h ],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 16, 16]
        )
        eta = ParametricFunction(lambda t : [ - 1.7 + t + ( 289/400 - ( t - 0.85)**2 )/ q + m , t - ( 289/400 - ( t - 0.85)**2 ) / q - m , h ], t_range = [0, 1.7], color = GREEN)
        zer = Sphere(radius = 0.05).set_color(GOLD).shift([         m ,     - m , h ])
        v   = Sphere(radius = 0.05).set_color(BLUE).shift([ - 1.7 + m ,     - m , h ])
        w   = Sphere(radius = 0.05).set_color(BLUE).shift([         m , 1.7 - m , h ])
        v0 = Line(zer.get_center(), v.get_center(), color = ORANGE)
        w0 = Line(zer.get_center(), w.get_center(), color = PURPLE)
        ar = Arrow([-0.9,-2.3,0], [0.9,-2.3,0], max_tip_length_to_length_ratio = 0.1 , stroke_width = 3)
        arla = Tex(r"$\exp_x$", font_size = 34).next_to(ar, UP).shift(0.15*DOWN)
        rad = 0.3
        angt = ParametricFunction(lambda t : [  m  + rad * np.cos(t) , -m + rad * np.sin(t) , h] , t_range = [PI/2 , PI], color = YELLOW )
        anga = ParametricFunction(lambda t : [ -m  + rad * np.cos(t) ,  m + rad * np.sin(t) , h] , t_range = [PI/2 , PI], color = YELLOW )
        fig = VGroup(s, x, y, z, txs, zer, v, w, yz, eta)
        self.add_fixed_in_frame_mobjects(vyli, wzli, gli, eli, pr2, pr3, ar, arla, pr4, pr5, pr6, pr7, g2li, e2li, vyll, wzll, call)
        self.remove(pr4, pr5, pr6, pr7, g2li, e2li, vyll, wzll, call)
        self.add(fig)
        self.wait()
        self.play(Write(pr4), Write(pr5), Create(g2li), Create(e2li))
        self.wait()
        self.play(Write(pr6), Create(vyll), Create(wzll), Create(call), Create(xy), Create(xz), Create(v0), Create(w0), Create(angt), Create(anga))
        self.wait()
        self.play(Create(pr7))
        self.wait()


class kpl(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/4)
        lem = Tex(r"Lemma", font_size = 34, color = BLUE).shift(5*LEFT + 2 * UP)
        le1 = Tex(r"If $K \geq 0$ in a neighborhood of $p$, then there is $r > 0$ such that", font_size = 34).next_to(lem, DOWN)
        le2 = Tex(r" if $x,y,z \in B(p, r)$, then ",  font_size = 34).next_to(le1, DOWN)
        le3 = Tex(r"$\measuredangle (x_y^z)$", r" $ \geq$ ", r"$\tilde{\measuredangle }(x_y^z)$", r".", font_size = 34).next_to(le2, DOWN)
        pr  = Tex(r"Exercise", font_size = 34, color = BLUE).next_to(le3, DOWN)
        pr1 = Tex(r"There is $r>0$ such that for all $x \in B(p, r)$, the map ", font_size = 34).next_to(pr, DOWN)
        pr2 = Tex(r"$\exp_x : T_x \Sigma \to \Sigma $ is non-singular and is injective in $B(0, 2r)$.", font_size = 34).next_to(pr1, DOWN)
        po  = Tex(r"Proof of Lemma", font_size = 34, color = BLUE).next_to(lem, DOWN).shift(1.5*UP)
        po1 = Tex(r"Pick $v,w \in T_x \Sigma$ with ", r"$\exp_x(v) = y$", r", ", r"$\exp_x(w) = z$", r", ", font_size = 34).next_to(po, DOWN)
        po2 = Tex(r"$\eta$", r" $:[0,1] \to T_x \Sigma$ the segment from $v$ to $w$, and ", r"$\gamma$", r" $: = \exp_x (\eta )$.", font_size = 34).next_to(po1, DOWN)
        po3 = Tex(r"Since $\exp_x$ is non-expanding,", font_size = 34).next_to(po2, DOWN)
        po4 = Tex(r"$d(u,v) =$ ", r"length$(\eta)$", r" $\geq $ ", r"length$(\gamma)$", r" $\geq d(y,z)$", font_size = 34).next_to(po3, DOWN)
        po5 = Tex(r"Also, ", r"$\vert v \vert = d(x,y)$", r", ", r"$\vert w \vert = d(x,z)$", r", ", r"$\measuredangle (0_y^z) = \measuredangle (x_y^z)$", r",", font_size = 34).next_to(po4, DOWN)
        po6 = Tex(r"so by monotonicity, $\measuredangle (x_y^z) \geq \tilde{\measuredangle}(x_y^z)$", font_size = 34).next_to(po5, DOWN)
        l = lem.get_left()[0]
        le1.shift((l - le1.get_left()[0])*RIGHT)
        le2.shift((l - le2.get_left()[0])*RIGHT)
        le3.shift(   le3.get_center()[0] *LEFT)
        pr.shift((l - pr.get_left()[0])*RIGHT)
        pr1.shift((l - pr1.get_left()[0])*RIGHT)
        pr2.shift((l - pr2.get_left()[0])*RIGHT)
        po.shift((l - po.get_left()[0])*RIGHT)
        po1.shift((l - po1.get_left()[0])*RIGHT)
        po2.shift((l - po2.get_left()[0])*RIGHT)
        po3.shift((l - po3.get_left()[0])*RIGHT)
        po4.shift(   po4.get_center()[0]*  LEFT)
        po5.shift((l - po5.get_left()[0])*RIGHT)
        po6.shift(   po6.get_center()[0]*  LEFT)
        toli = Line([0,0,0], [2  ,0,0], color = BLUE  ).next_to(le3, DOWN).shift(0.15*UP)
        vyli = Line([0,0,0], [1.6,0,0], color = BLUE  ).next_to(po1[1], DOWN).shift(0.15*UP)
        wzli = Line([0,0,0], [1.6,0,0], color = BLUE  ).next_to(po1[3], DOWN).shift(0.15*UP)
        etli = Line([0,0,0], [0.3,0,0], color = PINK  ).next_to(po2[0], DOWN).shift(0.15*UP)
        gli  = Line([0,0,0], [0.3,0,0], color = GREEN ).next_to(po2[2], DOWN).shift(0.15*UP)
        et2l = Line([0,0,0], [1.2,0,0], color = PINK  ).next_to(po4[1], DOWN).shift(0.15*UP)
        g2l  = Line([0,0,0], [1.2,0,0], color = GREEN ).next_to(po4[3], DOWN).shift(0.15*UP)
        vy2  = Line([0,0,0], [1.4,0,0], color = ORANGE).next_to(po5[1], DOWN).shift(0.15*UP)
        wz2  = Line([0,0,0], [1.4,0,0], color = PURPLE).next_to(po5[3], DOWN).shift(0.15*UP)
        aa2  = Line([0,0,0], [1.7,0,0], color = YELLOW).next_to(po5[5], DOWN).shift(0.15*UP)
        old = VGroup(lem, le1, le2, le3, toli, pr, pr1, pr2)
        #################################################################################
        r = 1.5
        m = 2.5
        h = -2.9
        hh = -2.2
        le = np.sin(PI/4)
        a = 1.9
        s = Surface(
            lambda u, v: [ r * np.sin(v) * np.cos(u) - m , r * np.sin(v) * np.sin(u) + m , r * np.cos(v) + h],
            u_range = [ 0, TAU ],
            v_range = [ 0, PI / 2 ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 32, 12]
        )
        x  = Sphere(radius = 0.05).set_color(GOLD).shift([        - m ,        + m ,                    r + h ])
        y  = Sphere(radius = 0.05).set_color(BLUE).shift([ le * r - m ,        + m , r * sqrt(1 - le **2) + h ])
        z  = Sphere(radius = 0.05).set_color(BLUE).shift([        - m , le * r + m , r * sqrt(1 - le **2) + h ])
        xy = ParametricFunction(lambda t: [    np.sin(t) * r - m ,                 m ,                                 np.cos(t)*r + h ], t_range = [0,PI/4], color = ORANGE)
        xz = ParametricFunction(lambda t: [                  - m , np.sin(t) * r + m ,                                 np.cos(t)*r + h ], t_range = [0,PI/4], color = PURPLE)
        yz = ParametricFunction(lambda t: [ le * r * (1 - t) - m ,    le * r * t + m , r * sqrt(1 - (le * (1-t))**2 - (le * t)**2) + h ], t_range = [0,  1 ],  color = GREEN )
        tps = Surface(
            lambda u, v: [ u + m, v - m , hh],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 16, 16]
        )
        xb  = Sphere(radius = 0.05).set_color(GOLD).shift([        m,        -m , hh ])
        yb  = Sphere(radius = 0.05).set_color(BLUE).shift([ 1.45 + m ,      - m , hh ])
        zb  = Sphere(radius = 0.05).set_color(BLUE).shift([        m , 1.45 - m , hh ])
        v0  = Line(xb.get_center(), yb.get_center(), color = ORANGE)
        w0  = Line(xb.get_center(), zb.get_center(), color = PURPLE)
        eta = Line(yb.get_center(), zb.get_center(), color = PINK  )
        ar = Arrow([-0.3,-2.1,0], [1.5,-2.1,0], max_tip_length_to_length_ratio = 0.1 , stroke_width = 3)
        arla = Tex(r"$\exp_x$", font_size = 34).next_to(ar, UP).shift(0.15*DOWN)
        fig = VGroup(s, x, y, z, tps, xb, yb, zb, eta, yz, ar, arla)
        ##################################################################################
        self.add_fixed_in_frame_mobjects(lem, le1, le2, le3, toli, pr, pr1, pr2, po, po1, po2, po3, po4, po5, po6, vyli, wzli, etli, gli, et2l, g2l, ar, arla, vy2, wz2, aa2)
        self.remove(lem, le1, le2, le3, toli, pr, pr1, pr2, po, po1, po2, po3, po4, po5, po6, vyli, wzli, etli, gli, et2l, g2l, ar, arla, vy2, wz2, aa2)
        self.play(Create(lem), Create(le1), Create(le2), Create(le3), Create(toli))
        self.wait()
        self.play(Create(pr), Create(pr1), Create(pr2))
        self.wait()
        self.play(Create(po), Create(po1), Create(po2), FadeOut(old), Create(vyli), Create(wzli), Create(etli), Create(gli), Create(fig))
        self.wait()
        self.play(Create(po3), Create(po4), Create(et2l), Create(g2l))
        self.wait()
        self.play(Create(po5), Create(vy2), Create(wz2), Create(aa2), Create(xy), Create(xz), Create(v0), Create(w0))
        self.wait()
        self.play(Create(po6))
        self.wait()



class topo(Scene):
    def construct(self):
        topo0 = Tex(r"Theorem (Toponogov)", color = BLUE, font_size = 34).to_edge(UL).shift(RIGHT + 1.5*DOWN)
        topo1 = Tex(r"Let $\Sigma$ be a complete surface such that for each $p \in \Sigma$,", font_size = 34).next_to(topo0, DOWN)
        topo2 = Tex(r"there is $r > 0 $ such that for any $a, b, c \in B(p, r)$,", font_size = 34).next_to(topo1, DOWN)
        topo3 = Tex(r"$\measuredangle (a_b^c) \geq \tilde{\measuredangle}(a_b^c)$.", font_size = 34).next_to(topo2, DOWN)
        topo4 = Tex(r"Then for each $x,y,z \in \Sigma $, one has", font_size = 34).next_to(topo3, DOWN)
        topo5 = Tex(r"$\measuredangle (x_y^z) \geq \tilde{\measuredangle}(x_y^z)$.", font_size = 34).next_to(topo4, DOWN)
        l = topo0.get_left()[0]
        topo1.shift((l - topo1.get_left()[0])*RIGHT)
        topo2.shift((l - topo2.get_left()[0])*RIGHT)
        topo3.shift( topo3.get_center()[0]*LEFT)
        topo4.shift((l - topo4.get_left()[0])*RIGHT)
        topo5.shift( topo5.get_center()[0]*LEFT)
        all = VGroup(topo0, topo1, topo2, topo3, topo4, topo5)
        self.play(Create(all))
        self.wait()





class cal(Scene):
    def construct(self):
        lem = Tex(r"Lemma", color = BLUE, font_size = 34).to_edge(UL).shift(0.5*RIGHT + 0.2*DOWN)
        le1 = Tex(r"For ", r"$u$", r", ", r"$w$", r" $\in \mathbb{R}^2$ with $\vert u \vert = 1$, the function $f: \mathbb{R} \to \mathbb{R}$", font_size = 34).next_to(lem, DOWN)
        le2 = Tex(r"given by $f(t) = \vert w + t u \vert$ is convex.", font_size = 34).next_to(le1, DOWN)
        pr = Tex(r"Proof", font_size = 34, color = BLUE).next_to(le2, DOWN)
        pr1 = Tex(r"$f (t) $", r" $= \sqrt{ w \cdot w + 2t u \cdot w + t^2  } $", font_size = 34).next_to(pr, DOWN).shift(0.1*UP)
        pr2 = Tex(r"$f^{\prime}(t) $", r" $ = \frac{1}{ \sqrt{ w \cdot w + 2t u \cdot w + t^2  } } \left(  u \cdot w + t   \right)   $", font_size = 34).next_to(pr1, DOWN)
        pr3 = Tex(r"$f^{\prime \prime }(t) $ ", r"$ =  \frac{1}{ \sqrt{ w \cdot w + 2t u \cdot w + t^2  } } - \frac{1}{ \sqrt{  w \cdot w + 2t u \cdot w + t^2 }^3 } \left(   u \cdot w + t   \right) ^2   $", font_size = 34).next_to(pr2, DOWN)
        pr4 = Tex(r"$ =  \frac{1}{ \sqrt{ w \cdot w + 2t u \cdot w + t^2  } ^3 } \left[  w \cdot w + 2t u \cdot w + t^2  - \left(   u \cdot w + t   \right) ^2 \right]  $", font_size = 34).next_to(pr3, DOWN)
        pr5 = Tex(r"$w \cdot w + $ ", r"$2tu \cdot w$", r" $ + $ ", r"$t^2$", r" $ - (u \cdot w)^2 -$ ", r"$ 2 t u \cdot w$", r" $  -$ ", r"$ t^2$", r" $ = $ ", r"$w \cdot w - (u \cdot w)^2 $ ", r" $\geq 0$", font_size = 34).next_to(pr4, DOWN)
        l = lem.get_left()[0]
        le1.shift((l - le1.get_left()[0])*RIGHT)
        le2.shift((l - le2.get_left()[0])*RIGHT)
        pr.shift((l - pr.get_left()[0])*RIGHT)
        pr1.shift(pr1.get_center()[0]*LEFT)
        pr2.shift(pr2.get_center()[0]*LEFT)
        pr3.shift(pr3.get_center()[0]*LEFT)
        pr4.shift((pr3[1].get_left()[0] - pr4[0].get_left()[0])*RIGHT)
        pr5.shift(pr5.get_center()[0]*LEFT + 0.2 *DOWN)
        wli  = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(le1[3], DOWN).shift(0.15*UP)
        uli  = Line([0,0,0], [0.3,0,0], color = YELLOW).next_to(le1[1], DOWN).shift(0.15*UP)
        can1 = Line([0,0,0], [1  , 0.2, 0], color = RED).next_to(pr5[1], DOWN).shift(0.4*UP)
        can2 = Line([0,0,0], [0.5, 0.2, 0], color = RED).next_to(pr5[3], DOWN).shift(0.5*UP)
        can3 = Line([0,0,0], [1  , 0.2, 0], color = RED).next_to(pr5[5], DOWN).shift(0.4*UP)
        can4 = Line([0,0,0], [0.5, 0.2, 0], color = RED).next_to(pr5[7], DOWN).shift(0.5*UP)
        pr55 = VGroup(pr5[0], pr5[1], pr5[2], pr5[3], pr5[4], pr5[5], pr5[6], pr5[7], pr5[8])
        ##########################################
        w = Dot([-1,-1,0], color = ORANGE)
        u = Arrow([1,-2,0], [3,-1,0], color = YELLOW)
        g = ParametricFunction(lambda t : [ -1 + 2 * t , t -1 , 0 ] , t_range = [-1,2], color = BLUE )
        ##########################################
        self.play(Create(lem), Create(le1), Create(le2))
        self.wait()
        self.play(Create(uli), Create(wli), Create(u), Create(w), Create(g))
        self.wait()
        self.play(Create(pr), Create(pr1), FadeOut(g), FadeOut(u), FadeOut(w))
        self.wait()
        self.play(Create(pr2))
        self.wait()
        self.play(Create(pr3))
        self.wait()
        self.play(Create(pr4))
        self.wait()
        self.play(Create(pr55))
        self.wait()
        self.play(Create(can1), Create(can3), Create(can2), Create(can4))
        self.wait()
        self.play(Create(pr5[9]))
        self.wait()
        self.play(Create(pr5[10]))
        self.wait()
        proof = VGroup(pr, pr1, pr2, pr3, pr4, pr5, can1, can2, can3, can4, uli, wli)
        self.play(FadeOut(proof))
        self.wait()



class anex(Scene):
    def construct(self):
        lem = Tex(r"Lemma", color = BLUE, font_size = 34).to_edge(UL).shift(0.5*RIGHT + 0.2*DOWN)
        le1 = Tex(r"For ", r"$u$", r", ", r"$w$", r" $\in \mathbb{R}^2$ with $\vert u \vert = 1$, the function $f: \mathbb{R} \to \mathbb{R}$", font_size = 34).next_to(lem, DOWN)
        le2 = Tex(r"given by $f(t) = \vert w + t u \vert$ is convex.", font_size = 34).next_to(le1, DOWN)
        co = Tex(r"Corollary", font_size = 34, color = BLUE).next_to(le2, DOWN)
        co1 = Tex(r"Let ", r"$a$", r", ", r"$b$", r", ", r"$c$", r", ", r"$d$", r" $ \in \mathbb{R}^2$ be distinct points with $d \in [a,c]$. Then", font_size = 34).next_to(co, DOWN)
        co2 = Tex(r"$\frac{\vert ab \vert + \vert ac \vert - \vert bc \vert }{\vert ac \vert} \leq \frac{\vert ab \vert + \vert ad \vert - \vert bd \vert }{\vert ad \vert} $", font_size = 34).next_to(co1, DOWN)
        co3 = co2.copy().shift((lem.get_center()[1]-co.get_center()[1])*UP)
        pr  = Tex(r"Proof", font_size = 34, color = BLUE).next_to(co3, DOWN)
        pr1 = Tex(r"Set ", r"$w$", r" $ = a$, ", r"$u$", r" $ = (c - a) / \vert c - a \vert $, and ", r"$\gamma (t)$", r" $ = w + t u$.", font_size = 34).next_to(pr, DOWN)
        pr2 = Tex(r"Then $h : \mathbb{R} \to \mathbb{R}$ given by", font_size = 34).next_to(pr1, DOWN)
        pr3 = Tex(r"$h(t) : = \vert ab \vert + t - \vert b \gamma (t) \vert $", font_size = 34).next_to(pr2, DOWN)
        pr4 = Tex(r"is concave and $h ( 0 ) = 0$. ", r"Then", font_size = 34).next_to(pr3, DOWN)
        pr5 = Tex(r"$\dfrac{h (\vert ac \vert )}{\vert ac \vert } \leq \dfrac{h (\vert ad \vert)}{\vert ad \vert}$", font_size = 34).next_to(pr4, DOWN)
        l = lem.get_left()[0]
        le1.shift((l - le1.get_left()[0])*RIGHT)
        le2.shift((l - le2.get_left()[0])*RIGHT)
        co.shift( (l -  co.get_left()[0])*RIGHT)
        co1.shift((l - co1.get_left()[0])*RIGHT)
        co2.shift(   co2.get_center()[0] *LEFT )
        pr.shift( (l -  pr.get_left()[0])*RIGHT)
        pr1.shift( (l -  pr1.get_left()[0])*RIGHT)
        pr2.shift( (l -  pr2.get_left()[0])*RIGHT)
        pr3.shift(     pr3.get_center()[0] *LEFT )
        pr4.shift( (l -  pr4.get_left()[0])*RIGHT)
        pr5.shift(     pr5.get_center()[0] *LEFT )
        ali  = Line([0,0,0], [0.25,0,0], color = ORANGE).next_to(co1[1], DOWN).shift(0.15*UP)
        bli  = Line([0,0,0], [0.25,0,0], color = PURPLE).next_to(co1[3], DOWN).shift(0.15*UP)
        cli  = Line([0,0,0], [0.25,0,0], color = GREEN ).next_to(co1[5], DOWN).shift(0.15*UP)
        dli  = Line([0,0,0], [0.25,0,0], color = RED   ).next_to(co1[7], DOWN).shift(0.15*UP)
        wli  = Line([0,0,0], [0.3 ,0,0], color = ORANGE).next_to(pr1[1], DOWN).shift(0.15*UP)
        uli  = Line([0,0,0], [0.3 ,0,0], color = YELLOW).next_to(pr1[3], DOWN).shift(0.15*UP)
        gli  = Line([0,0,0], [0.4 ,0,0], color = BLUE  ).next_to(pr1[5], DOWN).shift(0.15*UP)
        lemma = VGroup(lem, le1, le2)
        corollary = VGroup(co, co1, co2, ali, bli, cli, dli)
        ##########################################
        a = Dot([-2, -2.5 ,0], color = ORANGE)
        b = Dot([ 1, -0.7 ,0], color = PURPLE)
        c = Dot([ 2,  -3  ,0], color = GREEN )
        d = Dot([ 0,-2.75 ,0], color = RED   )
        ab  = Line(a.get_center(), b.get_center())
        ac  = Line(a.get_center(), c.get_center())
        bc  = Line(b.get_center(), c.get_center())
        u = Arrow([-1,-3,0], [1,-3.25,0], color = YELLOW)
        g = ParametricFunction(lambda t : [ -2 + 2 * t , -2.5 - 0.25 * t , 0 ] , t_range = [-0.5,2.5], color = BLUE )
        fig = VGroup(ab, bc, a,b,c,d,u,g)
        ##########################################
        self.add(lemma)
        self.wait()
        self.play(Create(co), Create(co1), Create(co2), Create(ab), Create(ac), Create(bc), Create(a), Create(b), Create(c), Create(d))
        self.wait()
        self.play(FadeOut(lemma), corollary.animate.shift((lem.get_center()[1]-co.get_center()[1] - 0.5)*UP))
        self.wait()
        self.play(Create(pr), Create(pr1), Create(u), Create(g), Create(gli), Create(wli), Create(uli))
        self.remove(ac)
        self.wait()
        self.play(Create(pr2), Create(pr3), Create(pr4[0]), fig.animate.shift(3.7*RIGHT + 0.2 *DOWN))
        self.wait()
        self.play(Create(pr4[1]), Create(pr5))
        self.wait()


class msd(ThreeDScene):
    def construct(self):
        packs = TexTemplate()
        packs.add_to_preamble(r"\usepackage{graphicx}")
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/4)
        de = Tex(r"Definition", font_size = 34, color = BLUE).shift(5.2*LEFT + 2.8 * UP)
        de1 = Tex(r"Let $x,y,z \in \Sigma$. The ", r"model side ",  r" is defined as", font_size = 34).next_to(de, DOWN)
        de2 = Tex(r"$ \tilde{\rotatebox[origin=c]{90}{$\prec$}} [ x_y^z ]  = d ( y^{\prime} , z^{\prime} ) $,",  tex_template = packs, font_size = 34).next_to(de1, DOWN)
        de22 = Tex(r"where $[x^{\prime} y^{\prime}z^{\prime}] \in \mathbb{R}^2$ is a triangle with", font_size = 34).next_to(de2, DOWN)
        de3 = Tex(r"$ d ( x^{\prime} , y^{\prime} ) = d (x,y) $, ", r"$ d ( x^{\prime} , z^{\prime} ) = d (x,z) $, ",  r"$ \measuredangle (x^{\prime z ^{\prime}}_{y ^{\prime}}) = \measuredangle (x_y^z)$", font_size = 34).next_to(de22, DOWN)
        prop = Tex(r"Proposition", font_size = 34, color = BLUE).next_to(de3, DOWN).shift(0.8*DOWN)
        pro1 = Tex(r"$\measuredangle (x_y^z) \geq \tilde{\measuredangle} (x_y^z)$ ", r" if and only if ", r" $ d (y , z ) \leq \tilde{\rotatebox[origin=c]{90}{$\prec$}} [ x_y^z ] $", tex_template = packs, font_size = 34).next_to(prop, DOWN)
        l = de.get_left()[0]
        de1.shift((l - de1.get_left()[0])*RIGHT)
        de2.shift( de2.get_center()[0]*LEFT )
        de22.shift((l - de22.get_left()[0])*RIGHT)
        de3.shift( de3.get_center()[0]*LEFT )
        prop.shift((l - prop.get_left()[0])*RIGHT)
        pro1.shift( pro1.get_center()[0]*LEFT )
        de3[0].shift(0.3*LEFT )
        de3[2].shift(0.3*RIGHT)
        de1[1].set_color(BLUE )
        pro1[0].shift(0.3*LEFT )
        pro1[2].shift(0.3*RIGHT)
        mali = Line([0,0,0], [2  ,0,0], color = GREEN ).next_to(de2, DOWN).shift(0.15*UP)
        xyli = Line([0,0,0], [2.2,0,0], color = ORANGE).next_to(de3[0], DOWN).shift(0.15*UP)
        xzli = Line([0,0,0], [2.2,0,0], color = PURPLE).next_to(de3[1], DOWN).shift(0.15*UP)
        thli = Line([0,0,0], [2  ,0,0], color = YELLOW).next_to(de3[2], DOWN).shift(0.15*UP)
        p111 = Line([0,0,0], [1.5,0,0], color = YELLOW).next_to(pro1[0], DOWN).shift(0.15*UP)
        p222 = Line([0,0,0], [1.5,0,0], color = GREEN ).next_to(pro1[2], DOWN).shift(0.15*UP)
        propo = VGroup(prop, pro1, p111, p222)
        #################################################################################
        r = 1.5
        m = 2.4
        h = -2.2
        hh = -1.6
        le = np.sin(PI/4)
        a = 1.9
        q = 1
        rad = 0.3
        s = Surface(
            lambda u, v: [ r * np.sin(v) * np.cos(u) - m , r * np.sin(v) * np.sin(u) + m , r * np.cos(v) + h],
            u_range = [ 0, TAU ],
            v_range = [ 0, PI / 2 ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 32, 12]
        )
        x  = Sphere(radius = 0.05).set_color(GOLD).shift([        - m ,        + m ,                    r + h ])
        y  = Sphere(radius = 0.05).set_color(BLUE).shift([ le * r - m ,        + m , r * sqrt(1 - le **2) + h ])
        z  = Sphere(radius = 0.05).set_color(BLUE).shift([        - m , le * r + m , r * sqrt(1 - le **2) + h ])
        xy = ParametricFunction(lambda t: [    np.sin(t) * r - m ,                 m ,                                 np.cos(t)*r + h ], t_range = [0,PI/4], color = ORANGE)
        xz = ParametricFunction(lambda t: [                  - m , np.sin(t) * r + m ,                                 np.cos(t)*r + h ], t_range = [0,PI/4], color = PURPLE)
        tps = Surface(
            lambda u, v: [ u + m, v - m , hh],
            u_range = [ - a , a ],
            v_range = [ - a , a ],
            checkerboard_colors = [RED, RED],
            fill_opacity = 0.15, 
            stroke_width = 0.2,
            resolution = [ 16, 16]
        )
        xb  = Sphere(radius = 0.05).set_color(GOLD).shift([        m - q,       -m - q, hh ])
        yb  = Sphere(radius = 0.05).set_color(BLUE).shift([ 1.55 + m - q,      - m - q, hh ])
        zb  = Sphere(radius = 0.05).set_color(BLUE).shift([        m - q, 1.55 - m - q, hh ])
        v0  = Line(xb.get_center(), yb.get_center(), color = ORANGE)
        w0  = Line(xb.get_center(), zb.get_center(), color = PURPLE)
        ms  = Line(yb.get_center(), zb.get_center(), color = GREEN)
        th1 = ParametricFunction(lambda t: [ -m + rad * np.cos(t) ,   m + rad * np.sin(t) , r + h  ], t_range = [0, PI / 2], color = YELLOW)
        th2 = ParametricFunction(lambda t: [  m + rad * np.cos(t) - q , - m + rad * np.sin(t) - q ,    hh  ], t_range = [0, PI / 2], color = YELLOW)
        fig0 = VGroup(s, x, y, z, tps, xy, xz, th1)
        fig1 = VGroup(xb, yb, zb, v0, w0, th2)
        ##################################################################################
        self.add_fixed_in_frame_mobjects(de, de1, de2, de3, mali, xyli, xzli, thli, de22, propo)
        self.remove(de, de1, de2, de3, mali, xyli, xzli, thli, de22, propo)
        self.play(Create(de), Create(de1), Create(de2), Create(de3), Create(mali), Create(xyli), Create(xzli), Create(thli), Create(de22), Create(fig0))
        self.wait()
        self.play(Create(fig1))
        self.wait()
        self.play(Create(ms))
        self.wait()
        self.play(FadeOut(fig1), FadeOut(fig0), FadeOut(ms), Create(propo))
        self.wait()



class chc2(Scene):
    def construct(self):
        packs = TexTemplate()
        packs.add_to_preamble(r"\usepackage{graphicx}")
        thm  = Tex(r"Theorem", font_size = 34, color = BLUE).to_edge(UL).shift(0.5*RIGHT + 1.5*DOWN)
        th1 = Tex(r"Let $\Sigma$ be a complete surface, and $[ x y z ]$ a triangle in $\Sigma$.", font_size = 34).next_to(thm, DOWN)
        th2 = Tex(r"$\bullet$ If $\Sigma$ has non-positive curvature and is simply connected, then", font_size = 34).next_to(th1, DOWN)
        th3 = Tex(r"$\tilde{\rotatebox{90}{$\prec$}} (x_y^z)$", r" $\leq $ ", r"$  d(y,z ) $", tex_template = packs, font_size = 34).next_to(th2, DOWN)
        th4 = Tex(r"$ \bullet$ If $\Sigma$ has non-negative curvature, then", font_size = 34).next_to(th3, DOWN)
        th5 = Tex(r"$\tilde{\rotatebox{90}{$\prec$}} (x_y^z)$", r" $\geq $ ", r"$d (y,z)$", tex_template = packs, font_size = 34).next_to(th4, DOWN)
        l = thm.get_left()[0]
        th1.shift((l - th1.get_left()[0])*RIGHT)
        th2.shift((l - th2.get_left()[0])*RIGHT)
        th4.shift((l - th4.get_left()[0])*RIGHT)
        th3.shift(   th3.get_center()[0] *LEFT )
        th5.shift(   th5.get_center()[0] *LEFT )
        aali  = Line( [0,0,0] , [0.8,0,0], color = TEAL   ).next_to(th3[2], DOWN).shift(0.15*UP)
        cali  = Line( [0,0,0] , [0.8,0,0], color = YELLOW ).next_to(th3[0], DOWN).shift(0.15*UP)
        aa2li = Line( [0,0,0] , [0.8,0,0], color = TEAL   ).next_to(th5[2], DOWN).shift(0.15*UP)
        ca2li = Line( [0,0,0] , [0.8,0,0], color = YELLOW ).next_to(th5[0], DOWN).shift(0.15*UP)
        all = VGroup(thm, th1, th2, th3, th4, th5, aali, cali, aa2li, ca2li)
        self.play(Create(all))
        self.wait(2)
        


class anex2(Scene):
    def construct(self):
        packs = TexTemplate()
        packs.add_to_preamble(r"\usepackage{graphicx}")
        exe = Tex(r"Exercise", color = BLUE, font_size = 34).shift(1.2*UP + 5*LEFT)
        ex1 = Tex(r"For $x, y,z, w \in \Sigma$ distinct points with $w \in [x,z]$, one has", font_size = 34).next_to(exe, DOWN)
        ex2 = Tex(r"$ \dfrac{d(x,y) + d(x,z) - \tilde{\rotatebox{90}{$\prec$}}[x_y^z ] }{d(x,z)} \leq \dfrac{d(x,y) + d(x,w) - \tilde{\rotatebox[origin=c]{90}{$\prec$}} [x_y^w]}{d(x,w)}   $", tex_template = packs, font_size = 34).next_to(ex1, DOWN)
        l = exe.get_left()[0]
        ex1.shift((l - ex1.get_left()[0])*RIGHT)
        ex2.shift(- ex2.get_center()[0]*RIGHT + 0.5*DOWN)
        exer = VGroup(exe, ex1, ex2)
        self.play(Create(exer))
        self.wait()



class indl(Scene):
    def construct(self):
        packs = TexTemplate()
        packs.add_to_preamble(r"\usepackage{graphicx}")
        il = Tex(r"Induction Lemma", color = BLUE, font_size = 34).to_edge(UL).shift(0.3*DOWN + 0.3 * RIGHT)
        il1 = Tex(r"Let $\Sigma$ be a complete surface, $x, y,z \in \Sigma$, and set $\ell : = d(x,y) + d(x,z)$.", font_size = 34).next_to(il, DOWN)
        il2 = Tex(r"Assume for all $a,b,c \in B(x, 10 \ell )$ with $d(a,b) + d (a,c) \leq 2 \ell / 3$ one has", tex_template = packs, font_size = 34).next_to(il1, DOWN)
        il3 = Tex(r"$d(b,c) \leq \tilde{ \rotatebox{90}{$\prec$}}[a_b^c]$.", tex_template = packs, font_size = 34).next_to(il2, DOWN)
        il4 = Tex(r"Then", tex_template = packs, font_size = 34).next_to(il3, DOWN)
        il5 = Tex(r"$d(y,z) \leq \tilde{\rotatebox{90}{$\prec$}}[x_y^z]$.", tex_template = packs, font_size = 34).next_to(il4, DOWN)
        pr  = Tex(r"Proof", font_size = 34, color = BLUE).next_to(il5, DOWN)
        pr1 = Tex(r"If $d(x,z) \geq d(x,y)$, ", r"take $x_1 \in [x,z]$", font_size = 34).next_to(pr, DOWN)
        pr2 = Tex(r"with ", r"$d(x,y) +  2 d(x, x_1) = 2 \ell / 3 $", font_size = 34).next_to(pr1, DOWN)
        l = il.get_left()[0]
        il1.shift((l - il1.get_left()[0])*RIGHT)
        il2.shift((l - il2.get_left()[0])*RIGHT)
        il3.shift(il3.get_center()[0]*LEFT)
        il4.shift((l - il4.get_left()[0])*RIGHT)
        il5.shift(il5.get_center()[0]*LEFT)
        pr.shift( (l -  pr.get_left()[0])*RIGHT)
        pr1.shift((l - pr1.get_left()[0])*RIGHT)
        pr2.shift((l - pr2.get_left()[0])*RIGHT)
        exer = VGroup(il, il1, il2, il3, il4, il5)
        surv = VGroup(pr1, pr2)
        old = VGroup(il, il1, il2, il3, il4, il5, pr)
        ########################################################################
        x = Dot([0, 0,0], color = BLUE)
        y = Dot([1, 1,0], color = BLUE)
        z = Dot([2,-1,0], color = BLUE)
        xx = x.copy()
        x1 = Dot([0.5,-0.25,0], color = BLUE)
        xl = Tex(r"$x$", font_size = 34).next_to(x, DOWN).shift(0.2*UP  +0.2*LEFT )
        yl = Tex(r"$y$", font_size = 34).next_to(y,   UP).shift(0.2*DOWN+0.2*RIGHT)
        zl = Tex(r"$z$", font_size = 34).next_to(z, DOWN).shift(0.2*UP  +0.2*RIGHT)
        x1l= Tex(r"$x_1$",font_size= 34).next_to(x1,DOWN).shift(0.1*UP  +0.1*LEFT )
        xy = Line(x.get_center(), y.get_center(), stroke_width= 2)
        xz = Line(x.get_center(), z.get_center(), stroke_width= 2)
        yz = Line(y.get_center(), z.get_center(), stroke_width= 2)
        fig = VGroup(xl, yl, zl, xy, xz, yz, x, y, z, xx, x1, x1l)
        fig.shift(2*RIGHT + DOWN)
        fig0 = VGroup(xy, xz, yz, xl, yl, zl, x,y,z)
        ########################################################################
        self.play(Create(exer))
        self.wait()
        self.play(Create(pr), Create(pr1[0]), Create(fig0))
        self.wait()
        self.play(Create(pr1[1]), Create(pr2))
        self.wait()
        self.play(Transform(xx,x1), run_time = 2)
        self.play(Create(x1l))
        self.wait()
        self.play(pr1.animate.shift( (il.get_center()[1] - pr1.get_center()[1] - 0.7 )*UP ), pr2.animate.shift( (il.get_center()[1] - pr2.get_center()[1] - 0.7 )*UP + (pr1.get_right() - pr2.get_left() + 0.15 ) *RIGHT )   , FadeOut(old))
        self.wait()





class ind1(Scene):
    def construct(self):
        packs = TexTemplate()
        packs.add_to_preamble(r"\usepackage{graphicx}")
        pr1 = Tex(r"If $d(x,z) \geq d(x,y)$, ", r"take $x_1 \in [x,z]$", font_size = 34).to_edge(UL).shift( DOWN + 0.3 * RIGHT)
        pr2 = Tex(r"with ", r"$d(x,y) +  2 d(x, x_1) = 2 \ell / 3 $", font_size = 34).next_to(pr1, DOWN)
        pr2.shift( (pr1.get_center()[1] - pr2.get_center()[1])*UP + (pr1.get_right() - pr2.get_left() + 0.15 ) *RIGHT  )
        pr3 = Tex(r"$r_0 : = d(x,y)+ d(x,z)$,", r"$r_1 : = d(x_1, y) + d(x_1, z)$", font_size = 34).next_to(pr1, DOWN)
        pr4 = Tex(r"$s_0 : =  \tilde{\rotatebox{90}{$\prec$}} [x_y^z] $,", r"$s_1 : =  \tilde{\rotatebox{90}{$\prec$}} [x_{1y}^{\, \text{ }  z}]$", tex_template = packs, font_size = 34).next_to(pr1, DOWN)
        pr5 = Tex(r"$r_1$", r" $ \leq $ ", r"$ r_0$", font_size = 34).next_to(pr2, DOWN)
        pr6 = Tex(r"Claim", r": $s_1 \leq s_0$", font_size = 34).next_to(pr4, DOWN)
        pr7 = Tex(r"$d(\overline{x}, \overline{y})= d(x,y)$, ", r"$d(\overline{x}_1, \overline{y})= d(x_1,y)$, ",r"$d(\overline{x}, \overline{x}_1)= d(x,x_1)$,", font_size = 34).next_to(pr6, DOWN)
        #pr8 = Tex(r"$z^{\prime}$ in the extension of $[\overline{x},\overline{x}_1]$ and $d(z^{\prime}, \overline{x})=d(z,x)$", font_size = 34).next_to(pr7, DOWN)
        #pr9 = Tex(r"and $d(z^{\prime}, \overline{x})=d(z,x)$", font_size = 34).next_to(pr8, DOWN)
        pr8 = Tex(r"$d(x, y) + d(x, x_1) $", r" $\leq 2 \ell / 3$", font_size = 34).next_to(pr7, DOWN)
        pr9= Tex(r"$d(x_1, y) + d(x_1, x)$", r" $ \leq 2 \ell / 3$", font_size = 34).next_to(pr8, DOWN)
        pr10= Tex(r"$\measuredangle (x_y^{x_1})$", r" $\geq$ ", r"$ \measuredangle (\overline{x} _{\overline{y}}^{\overline{x}_1})$", font_size = 34).next_to(pr9 , DOWN)
        pr11= Tex(r"$\measuredangle (x_{1y}^{ \,\text{ } x})$", r" $\geq$ ", r"$ \measuredangle (\overline{x} _{1\overline{y}}^{\,\text{ } \overline{x}})$", font_size = 34).next_to(pr10, DOWN)
        l = pr1.get_left()[0]
        pr3.shift( ( pr3.get_center()[0] + 2.2)*LEFT)
        pr3[0].shift(0.3*LEFT )
        pr3[1].shift(0.3*RIGHT)
        pr4.shift( ( pr4.get_center()[0] + 2.2)*LEFT)
        pr4[0].shift(0.3*LEFT )
        pr4[1].shift(0.3*RIGHT)
        pr5.shift(2*RIGHT + 0.5 * DOWN )
        pr6.shift( ( pr4.get_center()[0] -  pr6.get_center()[0])*RIGHT)
        pr7.shift( ( l -  pr7.get_left()[0])*RIGHT)
        pr8.shift( ( l -  pr8.get_left()[0])*RIGHT)
        pr9.shift( ( l -  pr9.get_left()[0])*RIGHT)
        pr8[1].shift((pr9[1].get_left()[0] - pr8[1].get_left()[0])*RIGHT)
        r1li = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(pr5[0], DOWN).shift(0.15*UP)
        r0li = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(pr5[2], DOWN).shift(0.15*UP)
        clli = Line([0,0,0], [1,0,0]).next_to(pr6[0], DOWN).shift(0.15*UP)
        cc1li = Line([0,0,0], [1.8,0,0], color = PINK).next_to(pr7[0], DOWN).shift(0.15*UP)
        cc2li = Line([0,0,0], [1.8,0,0], color = ORANGE).next_to(pr7[1], DOWN).shift(0.15*UP)
        cc3li = Line([0,0,0], [1.8,0,0], color = RED).next_to(pr7[2], DOWN).shift(0.15*UP)
        goldli = Line([0,0,0], [3.3,0,0], color = GOLD).next_to(pr2[1], DOWN).shift(0.15*UP)
        gol1li = Line([0,0,0], [2.1,0,0], color = GOLD).next_to(pr9[0], DOWN).shift(0.15*UP)
        gol2li = Line([0,0,0], [2.1,0,0], color = GOLD).next_to(pr8[0], DOWN).shift(0.15*UP)
        aa1li = Line([0,0,0], [1,0,0], color = PURPLE).next_to(pr10[0], DOWN).shift(0.15*UP)
        aa2li = Line([0,0,0], [1,0,0], color = PURPLE).next_to(pr11[0], DOWN).shift(0.15*UP)
        ca1li = Line([0,0,0], [1,0,0], color = ORANGE).next_to(pr10[2], DOWN).shift(0.15*UP)
        ca2li = Line([0,0,0], [1,0,0], color = ORANGE).next_to(pr11[2], DOWN).shift(0.15*UP)
        cla = VGroup(pr6, clli)
        ########################################################################
        #def direction(a,b):
        #    return [ b.get_center()[0] - a.get_center()[0] , b.get_center()[1] - a.get_center()[1] , 0 ]
        bu = 0.1
        x = Dot([0, 0,0], color = BLUE)
        y = Dot([1, 1,0], color = BLUE)
        z = Dot([2,-1,0], color = BLUE)
        xx = x.copy()
        x1 = Dot([0.5,-0.25,0], color = BLUE)
        xl = Tex(r"$x$", font_size = 34).next_to(x, DOWN).shift(0.2*UP  +0.2*LEFT )
        yl = Tex(r"$y$", font_size = 34).next_to(y,   UP).shift(0.2*DOWN+0.2*RIGHT)
        zl = Tex(r"$z$", font_size = 34).next_to(z, DOWN).shift(0.2*UP  +0.2*RIGHT)
        x1l= Tex(r"$x_1$",font_size= 34).next_to(x1,DOWN).shift(0.1*UP  +0.1*LEFT )
        xx1 = Line(x.get_center(), x1.get_center(), stroke_width= 2)
        xy = Line(x.get_center(), y.get_center(), stroke_width= 2)
        xz = Line(x.get_center(), z.get_center(), stroke_width= 2)
        yz = Line(y.get_center(), z.get_center(), stroke_width= 2)
        x1y = Line(x1.get_center(), y.get_center(), stroke_width= 2)
        r11= Line(x1.get_center(), z.get_center(), color = ORANGE, stroke_width= 4)
        r12= Line(x1.get_center(), y.get_center(), color = ORANGE, stroke_width= 4)
        r21= Line( x.get_center(), z.get_center(), color = PURPLE, stroke_width= 3)
        r22= Line( x.get_center(), y.get_center(), color = PURPLE, stroke_width= 3)
        r1 = VGroup(r11, r12)
        r2 = VGroup(r21, r22)
        fig = VGroup(xl, yl, zl, xy, xz, yz, x, y, z, xx, x1, x1l)
        fig0 = VGroup(xl, yl, zl, xy, xz, yz, x, y, z, xx, x1, x1l, r11, r12, r21, r22, x1y, xx1)
        fig0.shift(2*RIGHT + DOWN)
        dummy = fig.copy().shift(UP + 1.3 * RIGHT)
        sbox = SurroundingRectangle(dummy, buff = 0.2, stroke_width = 3)
        sl = Tex(r"$\Sigma$", font_size = 34).next_to(z, RIGHT).shift( 1.6 * UP + 1.7 * RIGHT)
        sli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(sl, DOWN).shift(0.15*UP)
        r2l = Tex(r"$\mathbb{R}^2$", font_size = 34).next_to(z, RIGHT).shift( DOWN + 1.5 * LEFT )
        r2li = Line([0,0,0], [0.4,0,0], color = GREEN).next_to(r2l, DOWN).shift(0.15*UP)
        xb  = Dot([-0.5,0.25, 0], color = BLUE)
        x1b = Dot([0,0,0], color = BLUE)
        yb = Dot([0.5,1.2,0], color = BLUE)
        zb = Dot([1.5,-0.75,0], color = BLUE)
        xyb = Line(xb.get_center(),  yb.get_center(), stroke_width= 2)
        xx1b= Line(xb.get_center(), x1b.get_center(), stroke_width= 2)
        yx1b= Line(yb.get_center(), x1b.get_center(), stroke_width= 2)
        zx1b= Line(zb.get_center(), x1b.get_center(), stroke_width= 2)
        yzb = Line(yb.get_center(),  zb.get_center(), stroke_width= 2)
        ybl = Tex(r"$\overline{y}$", font_size = 34).next_to(yb,   UP).shift(0.2*DOWN+0.2*RIGHT)
        zbl = Tex(r"$z^{\prime}$", font_size = 34).next_to(zb, DOWN).shift(0.2*UP  +0.2*RIGHT)
        x1bl= Tex(r"$\overline{x}_1$",font_size= 34).next_to(x1b, DOWN).shift(0.1*UP  +0.1*LEFT )
        xbl= Tex(r"$\overline{x}$",font_size= 34).next_to(xb, DOWN).shift(0.1*UP  +0.1*LEFT )
        model0 = VGroup(xb, x1b, yb, zb, xyb, xx1b, yx1b, zx1b, yzb, ybl, xbl, x1bl, zbl)
        model0.shift(0.6* RIGHT + 2.3 *DOWN)
        model = VGroup(xb, x1b, yb, xyb, xx1b, yx1b, ybl, xbl, x1bl)
        rbox = SurroundingRectangle(model0, buff = 0.2, stroke_width = 3, color = GREEN)
        obser = VGroup(pr5, r1li, r0li)
        ########################################################################
        self.add(fig, pr1, pr2)
        self.wait()
        self.play(Create(pr3))
        self.wait()
        self.play(Create(pr5), Create(r1), Create(r1li), Create(r0li))
        self.wait()
        self.play(Create(r2))
        self.wait()
        self.add(x1y)
        self.play(Create(pr4), FadeOut(pr5), FadeOut(r1li), FadeOut(r0li), FadeOut(pr3), FadeOut(r1), FadeOut(r2))
        self.wait()
        self.play(Create(pr6), Create(clli))
        self.wait()
        fig.add(x1y, xx1)
        self.play(fig.animate.shift(UP + 1.3 * RIGHT), Create(model), Create(rbox), Create(sbox), Create(sl), Create(r2l), Create(sli), Create(r2li), Create(pr7))
        tri1 = Polygon(x.get_center(), x1.get_center(), y.get_center(), fill_opacity = 0.3 , fill_color = GOLD , stroke_width = 0)
        tri2 = Polygon(xb.get_center(), x1b.get_center(), yb.get_center(), fill_opacity = 0.3 , fill_color = GOLD , stroke_width = 0 )
        self.play(Create(tri1), Create(tri2))
        self.wait()
        self.play(Create(pr9), Create(pr8), Create(goldli), Create(gol1li), Create(gol2li))
        self.wait()
        rrr = 0.2
        yxx1 = Angle(Line(x,x1), Line(x,y) , color = PURPLE, radius = 1.3 * rrr)
        xx1y = Angle(Line(x1,y), Line(x1,x), color = PURPLE, radius = rrr)
        yxx1b = Angle(Line(xb, x1b), Line(xb , yb), color = ORANGE, radius = 1.25 * rrr)
        xx1yb = Angle(Line(x1b, yb), Line(x1b,xb ), color = ORANGE, radius = 0.9 * rrr)
        allang = VGroup(yxx1, xx1y, yxx1b, xx1yb)
        self.play(Create(allang), Create(pr10), Create(pr11), Create(aa1li), Create(aa2li), Create(ca1li), Create(ca2li))
        self.wait()
        remo = VGroup(pr1, pr2, goldli, tri1, tri2, gol1li, gol2li, pr8, pr9)
        small1 = VGroup(pr10, aa1li, ca1li)
        small2 = VGroup(pr11, aa2li, ca2li)
        self.play( 
            pr4.animate.shift(LEFT+ 0.2 * LEFT + 0.1 *DOWN), 
            FadeOut(remo), 
            cla.animate.shift((pr4.get_center()[1]- pr6.get_center()[1] - 0.1 )*UP + (pr4.get_right()[0]- pr6.get_left()[0] - 0.3)*RIGHT),
            small1.animate.shift((pr6.get_center()[1] - pr10.get_center()[1]-0.2)*UP),
            small2.animate.shift((pr6.get_center()[1] - pr11.get_center()[1]-0.2)*UP + (pr10.get_right()[0] - pr11.get_left()[0] + 1 )*RIGHT),
            pr7.animate.shift((pr1.get_center()[1]-pr7.get_center()[1])*UP + 0.2 *RIGHT)
            )
        self.wait()








class ind2(Scene):
    def construct(self):
        packs = TexTemplate()
        packs.add_to_preamble(r"\usepackage{graphicx}")
        pr1 = Tex(r"If $d(x,z) \geq d(x,y)$, ", r"take $x_1 \in [x,z]$", font_size = 34).to_edge(UL).shift( DOWN + 0.3 * RIGHT)
        pr4 = Tex(r"$s_0 : =  \tilde{\rotatebox{90}{$\prec$}} [x_y^z] $,", r"$s_1 : =  \tilde{\rotatebox{90}{$\prec$}} [x_{1y}^{\, \text{ }  z}]$", tex_template = packs, font_size = 34).next_to(pr1, DOWN)
        pr6 = Tex(r"Claim", r": $s_1 \leq s_0$", font_size = 34).next_to(pr4, DOWN)
        pr7 = Tex(r"$d(\overline{x}, \overline{y})= d(x,y)$, ", r"$d(\overline{x}_1, \overline{y})= d(x_1,y)$, ",r"$d(\overline{x}, \overline{x}_1)= d(x,x_1)$,", font_size = 34).next_to(pr6, DOWN)
        pr8 = Tex(r"$z^{\prime}$ in the extension of $[\overline{x},\overline{x}_1]$ and $d(z^{\prime}, \overline{x})=d(z,x)$", font_size = 34).next_to(pr7, DOWN)
        pr10= Tex(r"$\measuredangle (x_y^{x_1})$", r" $\geq$ ", r"$ \measuredangle (\overline{x} _{\overline{y}}^{\overline{x}_1})$", font_size = 34).next_to(pr7 , DOWN)
        pr11= Tex(r"$\measuredangle (x_{1y}^{ \,\text{ } x})$", r" $\geq$ ", r"$ \measuredangle (\overline{x} _{1\overline{y}}^{\,\text{ } \overline{x}})$", font_size = 34).next_to(pr10, DOWN)
        pr12= Tex(r"$\measuredangle (x_{1y}^{\, \text{ } z})  $ ", r"$=  \pi -  $ ", r"$\measuredangle  ( x_{1y}^{\, \text{ } x} )$", font_size = 34).next_to(pr8, DOWN)
        pr13= Tex(r"$\leq  $ ", r"$ \pi -$ ", r"$\measuredangle (  \overline{x}_{1 \overline{y}}^{\, \text{ } \overline{x}} ) $", font_size = 34).next_to(pr12, DOWN)
        pr14= Tex(r"$ = $ ", r" $ \measuredangle (  \overline{x}_{1 \overline{y}}^{\, \text{ } \overline{z}} ) $", font_size = 34).next_to(pr13, DOWN)
        pr15= Tex(r"$ \Rightarrow $ $ \text{ }$ ", r"$\tilde{\rotatebox{90}{$\prec$}} [x_{1y}^{\, \text{ }  z}]$", r" $ \leq $ ", r"$ d (\overline{y}, z^{\prime} ) $", tex_template = packs, font_size = 34).next_to(pr12, DOWN)
        pr16= Tex(r"$ \measuredangle ( x_y^z ) = \measuredangle (x_y^{x_1})$", r" $ \geq $ ", r"$ \measuredangle (\overline{x}_{\overline{y}}^{\overline{x}_1} ) = \measuredangle (\overline{x}_{\overline{y}}^{z^{\prime}} ) $", tex_template = packs, font_size = 34).next_to(pr15, DOWN)
        pr17= Tex(r"$ \Rightarrow $ $ \text{ }$ ", r"$\tilde{\rotatebox{90}{$\prec$}} [x_{y}^{ z}]$", r" $ \geq $ ", r"$ d (\overline{y}, z^{\prime} ) $", tex_template = packs, font_size = 34).next_to(pr16, DOWN)
        l = pr1.get_left()[0]
        pr4.shift( ( pr4.get_center()[0] + 2.2)*LEFT)
        pr4[0].shift(0.3*LEFT )
        pr4[1].shift(0.3*RIGHT)
        pr6.shift( ( pr4.get_center()[0] -  pr6.get_center()[0])*RIGHT)
        pr7.shift( ( l -  pr7.get_left()[0])*RIGHT)
        pr8.shift( ( l -  pr8.get_left()[0])*RIGHT + 0.1 * UP)
        pr12.shift( ( l -  pr12.get_left()[0] + 1.2 )*RIGHT)
        pr13.shift( (  pr12[1].get_left()[0] - pr13.get_left()[0] )*RIGHT)
        pr14.shift( (  pr12[1].get_left()[0] - pr14.get_left()[0] )*RIGHT)
        pr15.shift( (  pr14.get_center()[0] - pr15.get_center()[0] - 0.6 )*RIGHT + 0.1 *DOWN)
        pr16.shift( (  pr15.get_center()[0] - pr16.get_center()[0] + 0.3 )*RIGHT + 0.2 *DOWN)
        pr17.shift( (  pr15.get_center()[0] - pr17.get_center()[0] )*RIGHT + 0.3 *DOWN)
        clli = Line([0,0,0], [1,0,0]).next_to(pr6[0], DOWN).shift(0.15*UP)
        aa1li = Line([0,0,0], [1,0,0], color = PURPLE).next_to(pr10[0], DOWN).shift(0.15*UP)
        aa2li = Line([0,0,0], [1,0,0], color = PURPLE).next_to(pr11[0], DOWN).shift(0.15*UP)
        ca1li = Line([0,0,0], [1,0,0], color = ORANGE).next_to(pr10[2], DOWN).shift(0.15*UP)
        ca2li = Line([0,0,0], [1,0,0], color = ORANGE).next_to(pr11[2], DOWN).shift(0.15*UP)
        aa3li = Line([0,0,0], [1,0,0], color = RED).next_to(pr12[2], DOWN).shift(0.15*UP)
        ca3li = Line([0,0,0], [1,0,0], color = RED).next_to(pr13[2], DOWN).shift(0.15*UP)
        aa4li = Line([0,0,0], [1.8,0,0], color = PURPLE).next_to(pr16[0], DOWN).shift(0.15*UP)
        ca4li = Line([0,0,0], [1.8,0,0], color = ORANGE).next_to(pr16[2], DOWN).shift(0.15*UP)
        s1li2 = Line([0,0,0], [0.8,0,0], color = BLUE_E).next_to(pr15[1], DOWN).shift(0.15*UP)
        yzli  = Line([0,0,0], [1,0,0], color = RED).next_to(pr15[3], DOWN).shift(0.15*UP)
        s0li2 = Line([0,0,0], [0.8,0,0], color = TEAL).next_to(pr17[1], DOWN).shift(0.15*UP)
        yz2li = Line([0,0,0], [1,0,0],   color = RED ).next_to(pr17[3], DOWN).shift(0.15*UP)
        cla = VGroup(pr6, clli)
        oldtext = VGroup(pr4, pr6, pr7, pr10, pr11, clli, aa1li, aa2li, ca1li, ca2li)
        small1 = VGroup(pr10, aa1li, ca1li)
        small2 = VGroup(pr11, aa2li, ca2li)
        ########################################################################
        #def direction(a,b):
        #    return [ b.get_center()[0] - a.get_center()[0] , b.get_center()[1] - a.get_center()[1] , 0 ]
        bu = 0.1
        x = Dot([0, 0,0], color = BLUE)
        y = Dot([1, 1,0], color = BLUE)
        z = Dot([2,-1,0], color = BLUE)
        xx = x.copy()
        x1 = Dot([0.5,-0.25,0], color = BLUE)
        xl = Tex(r"$x$", font_size = 34).next_to(x, DOWN).shift(0.2*UP  +0.2*LEFT )
        yl = Tex(r"$y$", font_size = 34).next_to(y,   UP).shift(0.2*DOWN+0.2*RIGHT)
        zl = Tex(r"$z$", font_size = 34).next_to(z, DOWN).shift(0.2*UP  +0.2*RIGHT)
        x1l= Tex(r"$x_1$",font_size= 34).next_to(x1,DOWN).shift(0.1*UP  +0.1*LEFT )
        xx1 = Line(x.get_center(), x1.get_center(), stroke_width= 2)
        xy = Line(x.get_center(), y.get_center(), stroke_width= 2)
        xz = Line(x.get_center(), z.get_center(), stroke_width= 2)
        yz = Line(y.get_center(), z.get_center(), stroke_width= 2)
        x1y = Line(x1.get_center(), y.get_center(), stroke_width= 2)
        fig0 = VGroup(xl, yl, zl, xy, xz, yz, x, y, z, xx, x1, x1l, x1y, xx1)
        fig0.shift(3.3*RIGHT)
        sbox = SurroundingRectangle(fig0, buff = 0.2, stroke_width = 3)
        sl = Tex(r"$\Sigma$", font_size = 34).next_to(z, RIGHT).shift( 0.6 * UP + 0.4 * RIGHT)
        sli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(sl, DOWN).shift(0.15*UP)
        r2l = Tex(r"$\mathbb{R}^2$", font_size = 34).next_to(z, RIGHT).shift( 2 * DOWN + 2.8 * LEFT )
        r2li = Line([0,0,0], [0.4,0,0], color = GREEN).next_to(r2l, DOWN).shift(0.15*UP)
        xb  = Dot([-0.5,0.25, 0], color = BLUE)
        x1b = Dot([0,0,0], color = BLUE)
        yb = Dot([0.5,1.2,0], color = BLUE)
        zb = Dot([1.5,-0.75,0], color = BLUE)
        xyb = Line(xb.get_center(),  yb.get_center(), stroke_width= 2)
        xx1b= Line(xb.get_center(), x1b.get_center(), stroke_width= 2)
        yx1b= Line(yb.get_center(), x1b.get_center(), stroke_width= 2)
        zx1b= Line(zb.get_center(), x1b.get_center(), stroke_width= 2)
        yzb = Line(yb.get_center(),  zb.get_center(), stroke_width= 2)
        ybl = Tex(r"$\overline{y}$", font_size = 34).next_to(yb,   UP).shift(0.2*DOWN+0.2*RIGHT)
        zbl = Tex(r"$z^{\prime}$", font_size = 34).next_to(zb, DOWN).shift(0.2*UP  +0.2*RIGHT)
        x1bl= Tex(r"$\overline{x}_1$",font_size= 34).next_to(x1b, DOWN).shift(0.1*UP  +0.1*LEFT )
        xbl= Tex(r"$\overline{x}$",font_size= 34).next_to(xb, DOWN).shift(0.1*UP  +0.1*LEFT )
        model0 = VGroup(xb, x1b, yb, zb, xyb, xx1b, yx1b, zx1b, yzb, ybl, xbl, x1bl, zbl)
        model0.shift(0.6* RIGHT + 2.3 *DOWN)
        model = VGroup(xb, x1b, yb, xyb, xx1b, yx1b, ybl, xbl, x1bl)
        rbox = SurroundingRectangle(model0, buff = 0.2, stroke_width = 3, color = GREEN)
        rrr = 0.2
        yxx1 = Angle(Line(x,x1), Line(x,y) , color = PURPLE, radius = 1.3 * rrr)
        xx1y = Angle(Line(x1,y), Line(x1,x), color = PURPLE, radius = rrr)
        yxx1b = Angle(Line(xb, x1b), Line(xb , yb), color = ORANGE, radius = 1.25 * rrr)
        xx1yb = Angle(Line(x1b, yb), Line(x1b,xb ), color = ORANGE, radius = 0.9 * rrr)
        allang = VGroup(yxx1, xx1y, yxx1b, xx1yb)
        ########################################################################
        self.add(allang, fig0, model, rbox, sbox, sl, r2l, sli, r2li, oldtext)
        self.wait()
        self.play( 
            pr4.animate.shift(LEFT+ 0.2 * LEFT + 0.1 *DOWN), 
            cla.animate.shift((pr4.get_center()[1]- pr6.get_center()[1] - 0.1 )*UP + (pr4.get_right()[0]- pr6.get_left()[0] - 0.3)*RIGHT),
            small1.animate.shift((pr6.get_center()[1] - pr10.get_center()[1]-0.2)*UP),
            small2.animate.shift((pr6.get_center()[1] - pr11.get_center()[1]-0.2)*UP + (pr10.get_right()[0] - pr11.get_left()[0] + 1 )*RIGHT),
            pr7.animate.shift((pr1.get_center()[1]-pr7.get_center()[1])*UP + 0.2 *RIGHT)
            )
        self.wait()
        self.play(Create(pr8), Create(zb), Create(zx1b), Create(zbl), FadeOut(allang), FadeOut(aa1li), FadeOut(aa2li), FadeOut(ca1li), FadeOut(ca2li))
        self.wait()
        aa2li.set_color(RED)
        ca2li.set_color(RED)
        self.play(Create(pr12))
        self.wait()
        self.play(Create(pr13), Create(aa2li), Create(ca2li), Create(aa3li), Create(ca3li))
        self.wait()
        self.play(Create(pr14))
        self.wait()
        trash = VGroup(pr14[0] , pr13[1], pr13[2], pr12[1], pr12[2] , ca3li, aa3li, aa2li, ca2li)
        self.play(pr14[1].animate.shift((pr12.get_center()[1]-pr14.get_center()[1])*UP), pr13[0].animate.shift((pr12.get_center()[1]-pr13.get_center()[1])*UP), FadeOut(trash))
        self.wait()
        s1li  = Line([0,0,0], [1.2,0,0], color = BLUE_E).next_to(pr4[1], DOWN).shift(0.15*UP)
        s0li  = Line([0,0,0], [1.2,0,0], color = TEAL  ).next_to(pr4[0], DOWN).shift(0.15*UP)
        self.play(Create(pr15), Create(s1li), Create(s1li2), Create(yzli), yzb.animate.set_color(RED))
        self.wait()
        self.play(Create(pr16), Create(aa4li), Create(ca4li), Create(aa1li), Create(ca1li), Create(yxx1), Create(yxx1b))
        self.wait()
        self.play(Create(pr17), Create(s0li), Create(s0li2), Create(yz2li))
        self.wait()
        all = VGroup(pr7, sl, sli, r2l, r2li, pr16, pr14[1], pr12[0], pr13[0], pr8, pr10, pr11, aa1li, ca1li, aa4li, ca4li, yxx1, yxx1b, zb, zbl, yzb, zx1b, sbox, rbox)
        self.play(FadeOut(all), FadeOut(fig0), FadeOut(model))
        self.wait()


class ind3(Scene):
    def construct(self):
        packs = TexTemplate()
        packs.add_to_preamble(r"\usepackage{graphicx}")
        pr1 = Tex(r"Construct $x_1, x_2, x_3, \ldots$", font_size = 34).to_edge(UL).shift( 0.7*DOWN + 0.3 * RIGHT)
        pr2 = Tex(r"If $d(x_n, y) \leq d(x_n, z)$, then $x_{n+1} \in [x_n, z]$", font_size = 34).next_to(pr1, DOWN)
        pr3 = Tex(r"with $d(x_n , y) + 2 d(x_n, x_{n+1}) = 2 \ell / 3$.", font_size = 34).next_to(pr2, DOWN)
        pr4 = Tex(r"If $d(x_n, z) < d(x_n, y)$, then $x_{n+1} \in [x_n, y]$", font_size = 34).next_to(pr3, DOWN)
        pr5 = Tex(r"with $d(x_n , z) + 2 d(x_n, x_{n+1}) = 2 \ell / 3$.", font_size = 34).next_to(pr4, DOWN)
        pr6 = Tex(r"$r_n := d(x_n, y) + d(x_n, z)$ ", r" $s_n : = \tilde{\rotatebox{90}{$\prec$}}[x_{ny}^{\, \text{ } z}]$", tex_template = packs,  font_size = 34).next_to(pr5, DOWN)
        pr7 = Tex(r"Case 1", r": $r_n \leq 2 \ell / 3$ for some $n$", font_size = 34).next_to(pr6, DOWN)
        pr8 = Tex(r"$\Rightarrow$ ", r"$s_n$", r" $ \geq d(y,z)$", font_size = 34).next_to(pr7, DOWN)
        pr9 = Tex(r"$s_0 \geq s_1 \geq \ldots \geq s_n \geq d(y,z)$", font_size = 34).next_to(pr8, DOWN)
        l = pr1.get_left()[0]
        pr1.shift((l - pr1.get_left()[0])*RIGHT)
        pr2.shift((l - pr2.get_left()[0])*RIGHT)
        pr3.shift((l - pr3.get_left()[0])*RIGHT)
        pr4.shift((l - pr4.get_left()[0])*RIGHT)
        pr5.shift((l - pr5.get_left()[0])*RIGHT)
        pr6.shift((l - pr6.get_left()[0])*RIGHT)
        pr6[1].shift(0.1*RIGHT)
        pr8.shift((l - pr8.get_left()[0] + 1.5)*RIGHT)
        pr9.shift((pr8.get_center()[0] - pr9.get_center()[0])*RIGHT)
        pr7.shift((pr8.get_center()[0] - pr7.get_center()[0])*RIGHT)
        cali = Line([0,0,0], [0.7,0,0]).next_to(pr7[0], DOWN).shift(0.15*UP)
        rnli = Line([0,0,0], [3  ,0,0], color = ORANGE).next_to(pr6[0], DOWN).shift(0.15*UP)
        snli = Line([0,0,0], [1.5,0,0], color = GREEN).next_to(pr6[1], DOWN).shift(0.15*UP)
        sn2li= Line([0,0,0], [0.4,0,0], color = GREEN).next_to(pr8[1], DOWN).shift(0.15*UP)
        text = VGroup(pr1, pr2, pr3, pr4, pr5, pr7, pr8, pr9, cali, sn2li)
        survi = VGroup(pr6, rnli, snli) 
        #####################################################################
        y  = Dot([0,0,0],       color = BLUE)
        x  = Dot([0,3,0],       color = BLUE)
        z  = Dot([4,-1,0],      color = BLUE)
        x1 = Dot([5/4 , 7/4, 0], color = BLUE)
        x2 = Dot([5/2 , 1/2, 0], color = BLUE)
        x3 = Dot([5/3 , 1/3, 0], color = BLUE)
        x4 = Dot([17/6,-1/3, 0], color = BLUE)
        dots = VGroup(x, y, z, x1, x2, x3, x4)
        dots.shift(1.5 * DOWN + 1.5*RIGHT)
        xl  = Tex(r"$x$"  , font_size = 34).next_to(x , LEFT ).shift(0.2*UP   + 0.2*RIGHT)
        yl  = Tex(r"$y$"  , font_size = 34).next_to(y , LEFT ).shift(0.2*DOWN + 0.2*RIGHT)
        zl  = Tex(r"$z$"  , font_size = 34).next_to(z , RIGHT).shift(0.2*DOWN + 0.2*LEFT )
        x1l = Tex(r"$x_1$", font_size = 34).next_to(x1, UP   ).shift(0.2*DOWN + 0.2*RIGHT)
        x2l = Tex(r"$x_2$", font_size = 34).next_to(x2, UP   ).shift(0.2*DOWN + 0.2*RIGHT)
        x3l = Tex(r"$x_3$", font_size = 34).next_to(x3, UP   ).shift(0.2*DOWN + 0.2*LEFT )
        x4l = Tex(r"$x_4$", font_size = 34).next_to(x4, UP   ).shift(0.2*DOWN + 0.1*RIGHT)
        w  = x.copy()
        w1 = x1.copy()
        w2 = x2.copy()
        w3 = x3.copy()
        xy  = Line(x.get_center() , y.get_center())
        xz  = Line(x.get_center() , z.get_center())
        yz  = Line(y.get_center() , z.get_center())
        x1y = Line(x1.get_center(), y.get_center())
        x2y = Line(x2.get_center(), y.get_center())
        x3z = Line(x3.get_center(), z.get_center())
        x4y = Line(x4.get_center(), y.get_center())
        original = VGroup(x,y,z, xy, xz, yz, xl, yl, zl)
        self.play(Create(original), Create(pr1))
        self.wait()
        self.play(Transform(w, x1))
        self.play(Create(x1y), Create(x1l))
        self.play(Create(pr2), Create(pr3), Create(pr4), Create(pr5))
        self.play(Transform(w1, x2))
        self.play(Create(x2y), Create(x2l))
        self.play(Transform(w2, x3))
        self.play(Create(x3z), Create(x3l))
        self.play(Transform(w3, x4))
        self.play(Create(x4y), Create(x4l))
        self.wait()
        self.play(Create(pr6), Create(rnli), Create(snli))
        self.wait()
        oldf = VGroup(x,y,z,xy,xz,yz, xl, yl, zl, x1y, x2y, x3z, x4y, x1l, x2l, x3l, x4l, x1, x2, x3, x4, w, w1, w2, w3)
        ###########################################################
        yn  = Dot([0,0,0], color = BLUE)
        xn  = Dot([1,2,0], color = BLUE)
        zn  = Dot([3,0,0], color = BLUE)
        xnyn = Line(xn, yn, color = ORANGE)
        xnzn = Line(xn, zn, color = ORANGE)
        ynzn = Line(zn, yn)
        xnl  = Tex(r"$x_n$", font_size = 34).next_to(xn , LEFT ).shift(0.2*UP   + 0.2*RIGHT)
        ynl  = Tex(r"$y$"  , font_size = 34).next_to(yn , LEFT ).shift(0.2*DOWN + 0.2*RIGHT)
        znl  = Tex(r"$z$"  , font_size = 34).next_to(zn , RIGHT).shift(0.2*DOWN + 0.2*LEFT )
        newf = VGroup(xn, yn, zn, xnyn, xnzn, ynzn, xnl, ynl, znl)
        newf.shift(2*RIGHT + 1.5 *DOWN)
        self.play(FadeOut(oldf), Create(newf), Create(pr7), Create(cali))
        self.wait()
        self.play(Create(pr8), Create(sn2li))
        self.wait()
        self.play(Create(pr9))
        self.wait()
        self.play(FadeOut(text), survi.animate.shift((pr1.get_center()[1]-pr6.get_center()[1] - 0.1)*UP))
        self.wait()


class ind4(Scene):
    def construct(self):
        packs = TexTemplate()
        packs.add_to_preamble(r"\usepackage{graphicx}")
        pr6 = Tex(r"$r_n := d(x_n, y) + d(x_n, z)$ ", r" $s_n : = \tilde{\rotatebox{90}{$\prec$}}[x_{ny}^{\, \text{ } z}]$", tex_template = packs,  font_size = 34).to_edge(UL).shift( 0.5*DOWN + 0.3 * RIGHT)
        pr7 = Tex(r"Case 2", r": $r_n > 2 \ell / 3$ for all $n$", font_size = 34).next_to(pr6, DOWN)
        #pr8 = Tex(r"$\lim_{n \to \infty } r_n = r_{\infty} $ ", r" $ \geq d(y,z) $", font_size = 34).next_to(pr7, DOWN)
        pr8 = Tex(r"Claim", r": $r_n - s_n  \leq 18 (r_n - r_{n + 1}) $", font_size = 34).next_to(pr7, DOWN)
        pr9 = Tex(r"Assume $d(x_n, y) \leq d(x_n , z)$, so $x_{n+1} \in [x_n, z]$", font_size = 34).next_to(pr8, DOWN)
        pr10 = Tex(r"$ d (x_n, x_{n + 1}) \geq \ell / 18 $, ", r" $ \tilde{\rotatebox{90}{$\prec$}}[x_{n\,\text{ } y}^{\, \text{ } \, x_{n+1}}] \geq d(y, x_{n+1}) $", tex_template = packs, font_size = 34).next_to(pr9, DOWN)
        pr11 = Tex(r"$ \dfrac{ d(x_n, y) + d(x_n, z) - \tilde{\rotatebox{90}{$\prec$}}[x_{ny}^{\, \text{ }  z}]  }{d(x_n, z)} $", r" $ \leq $ ", r"$ \dfrac{d(x_n, y) + d(x_n, x_{n+1}) -  \tilde{\rotatebox{90}{$\prec$}}[x_{n\,\text{ } y}^{\, \text{ } \, x_{n+1}} ] }{d(x_n, x_{n+1})} $", tex_template = packs, font_size = 34).next_to(pr10, DOWN)
        pr11.shift(0.2*DOWN)
        pr12 = Tex(r"$ \frac{ 1}{ \ell } [ r_n  - s_n ] $", r" $ \frac{18}{\ell} [ d(x_n, y) + d(x_n, x_{n+1}) - d(y,x_{n+1}) ] $", tex_template = packs, font_size = 34).next_to(pr11, DOWN)
        pr14 = Tex(r"$  r_n  - s_n  $ " ,r"$ \leq 18 [ d(x_n, y) + d(x_n, x_{n+1}) - d(y,x_{n+1})  + $ ", r"$d(x_{n+1}, z) - d(x_{n+1}, z)$", r" $] $ ",  font_size = 34).next_to(pr11, DOWN)
        pr15 = Tex(r"$= 18 [ $ ", r"$r_n$", r" $ - $ ", r"$r_{n+1}$", r" $]$", font_size = 34).next_to(pr14, DOWN)
        l = pr6.get_left()[0]
        pr6[1].shift(0.1*RIGHT)
        pr7.shift(( l - pr7.get_left()[0] + 1 )*RIGHT)
        pr7[0].shift(0.1*LEFT)
        pr8.shift(( pr7.get_center()[0] - pr8.get_center()[0])*RIGHT)
        pr9.shift(( l - pr9.get_left()[0])*RIGHT)
        pr10.shift(( l - pr10.get_left()[0])*RIGHT)
        pr11.shift( pr11.get_center()[0]*LEFT)
        pr14.shift( pr14.get_center()[0]*LEFT)
        pr15.shift( ( pr14[1].get_left()[0] - pr15.get_left()[0]) *RIGHT)
        pr10[1].shift(0.2*RIGHT)
        pr12.shift(0.2*DOWN)
        pr12[0].shift((pr11[0].get_center()[0] - pr12[0].get_center()[0])*RIGHT + 0.1 *DOWN)
        pr12[1].shift((pr11[2].get_center()[0] - pr12[1].get_center()[0])*RIGHT + 0.1 *DOWN)
        pr13 = Tex(r"$ \leq $", font_size = 34).next_to(pr12[0], RIGHT).shift(0.3*RIGHT)
        pr112 = Tex(r"$\tilde{\rotatebox{90}{$\leq$}}$", tex_template = packs, font_size = 34).next_to(pr12[0],UP).shift(0.18*DOWN)
        pr1122 = Tex(r"$\tilde{\rotatebox{90}{$\geq$}}$", tex_template = packs, font_size = 34).next_to(pr12[1],UP).shift(0.18*DOWN)
        rnli = Line([0,0,0], [3  ,0,0], color = ORANGE).next_to(pr6[0], DOWN).shift(0.15*UP)
        snli = Line([0,0,0], [1.5,0,0], color = GREEN).next_to(pr6[1], DOWN).shift(0.15*UP)
        cali = Line([0,0,0], [0.7,0,0]).next_to(pr7[0], DOWN).shift(0.15*UP)
        clli = Line([0,0,0], [0.7,0,0]).next_to(pr8[0], DOWN).shift(0.15*UP)
        rn2li = Line([0,0,0], [2.6,0,0], color = ORANGE).shift(0.95*DOWN + 4.85*LEFT)
        sn2li = Line([0,0,0], [0.6,0,0], color = GREEN ).shift(DOWN + 1.5* LEFT) 
        rn3li = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(pr12[0], DOWN).shift(0.15*UP + 0.3*LEFT )
        sn3li = Line([0,0,0], [0.3,0,0], color = GREEN ).next_to(pr12[0], DOWN).shift(0.15*UP + 0.5*RIGHT)
        xzli1 = Line([0,0,0], [1  ,0,0], color = YELLOW).next_to(pr11[0], DOWN).shift(0.2*UP)
        xzli2 = Line([0,0,0], [0.3,0,0], color = YELLOW).next_to(pr12[0], LEFT).shift(0.45*RIGHT + 0.3*DOWN) 
        csli1 = Line([0,0,0], [1.4,0,0], color = PURPLE).next_to(pr11[2], RIGHT).shift(1.7*LEFT + 0.05*UP) 
        cs2i1 = Line([0,0,0], [3  ,0,0], color = PURPLE).next_to(pr10[1], DOWN).shift(0.15*UP) 
        cs3i1 = Line([0,0,0], [1.4,0,0], color = PURPLE).next_to(pr12[1], RIGHT).shift(1.7*LEFT + 0.3*DOWN) 
        dxx1  = Line([0,0,0], [2.2,0,0], color = RED).next_to(pr10[0], DOWN).shift(0.15*UP) 
        dxx2  = Line([0,0,0], [1.3,0,0],  color = RED).next_to(pr11[2], DOWN).shift(0.2*UP)
        dxx3 = Line([0,0,0], [0.3,0,0], color = RED).next_to(pr12[1], LEFT).shift(0.52*RIGHT + 0.3*DOWN) 
        finalcan   = Line([0,0,0], [3,0,0],  color = PINK).next_to(pr14[2], DOWN).shift(0.15*UP)
        finalcana  = Line([0,0,0], [3,0,0] ,  color = BLUE_E).next_to(pr14, DOWN).shift(0.15*UP+2.1*LEFT)
        finalcanb  = Line([0,0,0], [1.4,0,0],  color = TEAL).next_to(pr14, DOWN).shift(0.15*UP+0.9*RIGHT)
        finalcanc  = Line([0,0,0], [1.4,0,0],  color = BLUE_E).next_to(pr14, DOWN).shift(0.15*UP+ 2.8* RIGHT)
        finalcand  = Line([0,0,0], [1.4,0,0],  color = TEAL).next_to(pr14, DOWN).shift(0.15*UP+5*RIGHT)
        finalcane  = Line([0,0,0], [0.3,0,0],  color = BLUE_E).next_to(pr15[1], DOWN).shift(0.15*UP)
        finalcanf  = Line([0,0,0], [0.3,0,0],  color = TEAL).next_to(pr15[3], DOWN).shift(0.15*UP)
        annoy = VGroup(pr11, pr112, pr1122, xzli1, csli1, rn2li, sn2li, dxx2)
        important = VGroup(pr12, pr13, rn3li, sn3li, xzli2, cs3i1, dxx3)
        old = VGroup(pr6, rnli, snli) 
        proofc = VGroup(pr12, pr13, pr14, pr15, pr9, pr10, rn3li, sn3li, xzli2, cs3i1, dxx3, finalcana, finalcanb, finalcanc, finalcand, finalcane, finalcanf, cs2i1, dxx1)
        #####################################################################
        yn  = Dot([0,0,0],       color = BLUE)
        xn  = Dot([1,2,0],       color = BLUE)
        zn  = Dot([3,0,0],       color = BLUE)
        xnn = Dot([5/3, 4/3, 0], color = BLUE)
        xnyn = Line(xn, yn, color = ORANGE)
        xnzn = Line(xn, zn, color = ORANGE)
        ynzn = Line(zn, yn)
        xnny = Line(xnn, yn)
        xnl  = Tex(r"$x_n$", font_size = 34).next_to(xn , LEFT ).shift(0.2*UP   + 0.2*RIGHT)
        ynl  = Tex(r"$y$"  , font_size = 34).next_to(yn , LEFT ).shift(0.2*DOWN + 0.2*RIGHT)
        znl  = Tex(r"$z$"  , font_size = 34).next_to(zn , RIGHT).shift(0.2*DOWN + 0.2*LEFT )
        xnnl = Tex(r"$x_{n+1}$", font_size = 34).next_to(xnn , RIGHT).shift(0.2*UP + 0.2*LEFT )
        newf = VGroup(xn, yn, zn, xnyn, xnzn, ynzn, xnl, ynl, znl, xnn, xnnl, xnny)
        newf.shift(2*RIGHT + 1.5 *DOWN)
        newf0 = VGroup(xn, yn, zn, xnyn, xnzn, ynzn, xnl, ynl, znl)
        ###########################################################
        self.add(old, newf0)
        self.wait()
        self.play(Create(pr7), Create(cali))
        self.wait()
        self.play(Create(pr8), Create(clli))
        self.wait()
        self.play(Create(pr9), Create(xnn), Create(xnnl), xnyn.animate.set_color(WHITE), xnzn.animate.set_color(WHITE) , Create(xnny))
        self.wait()
        self.play(Create(pr10[0]))
        self.wait()
        self.play(Create(pr10[1]))
        self.wait()
        self.play(Create(pr11), newf.animate.shift(2.05*UP) )
        self.wait()
        self.play(Create(pr12[0]), Create(pr112), Create(rn2li), Create(sn2li), Create(xzli1), Create(xzli2), Create(rn3li), Create(sn3li))
        self.wait()
        compa = VGroup(csli1, cs2i1, cs3i1, dxx1, dxx2, dxx3)
        self.play(Create(compa), Create(pr12[1]), Create(pr1122))
        self.wait()
        self.play(Create(pr13))
        self.wait()
        self.play(FadeOut(annoy), important.animate.shift((1.4)*UP))
        self.wait()
        self.play(Create(pr14), Create(finalcan))
        self.wait()
        self.play(FadeOut(finalcan))
        self.wait()
        finallines = VGroup(finalcana, finalcanb, finalcanc, finalcand, finalcane, finalcanf)
        self.play(Create(pr15), Create(finallines))
        self.wait()
        self.play(FadeOut(proofc))
        self.wait()




class ind5(Scene):
    def construct(self):
        packs = TexTemplate()
        packs.add_to_preamble(r"\usepackage{graphicx}")
        pr6  = Tex(r"$r_n := d(x_n, y) + d(x_n, z)$ ", r" $s_n : = \tilde{\rotatebox{90}{$\prec$}}[x_{ny}^{\, \text{ } z}]$", tex_template = packs,  font_size = 34).to_edge(UL).shift( 0.5*DOWN + 0.3 * RIGHT)
        pr7  = Tex(r"Case 2", r": $r_n > 2 \ell / 3$ for all $n$", font_size = 34).next_to(pr6, DOWN)
        pr8  = Tex(r"Claim", r": $r_n - s_n  \leq 18 (r_n - r_{n + 1}) $", font_size = 34).next_to(pr7, DOWN)
        pr9  = Tex(r"For each $n$, one has $s_n \leq r_n$, so ", font_size = 34).next_to(pr8, DOWN)
        pr99 = Tex(r"$\vert r_n - s_n \vert \to \infty$ as $n \to \infty$", font_size = 34).next_to(pr9, DOWN)
        pr10  = Tex(r"$s_0 \geq \lim_{n \to \infty } s_n = \lim_{n \to \infty } r_n  \geq d(y,z) $", font_size = 34).next_to(pr99, DOWN)
        l = pr6.get_left()[0]
        pr6[1].shift(0.1*RIGHT)
        pr7.shift(( l - pr7.get_left()[0] + 1 )*RIGHT)
        pr7[0].shift(0.1*LEFT)
        pr8.shift(( pr7.get_center()[0] - pr8.get_center()[0])*RIGHT)
        pr9.shift( ( pr7.get_center()[0] -  pr9.get_center()[0])*RIGHT)
        pr99.shift( ( pr7.get_center()[0] -  pr99.get_center()[0])*RIGHT)
        pr10.shift(pr10.get_center()[0]*LEFT + DOWN)
        rnli = Line([0,0,0], [3  ,0,0], color = ORANGE).next_to(pr6[0], DOWN).shift(0.15*UP)
        snli = Line([0,0,0], [1.5,0,0], color = GREEN).next_to(pr6[1], DOWN).shift(0.15*UP)
        cali = Line([0,0,0], [0.7,0,0]).next_to(pr7[0], DOWN).shift(0.15*UP)
        clli = Line([0,0,0], [0.7,0,0]).next_to(pr8[0], DOWN).shift(0.15*UP)
        box = SurroundingRectangle(pr10, buff = 0.2)
        still = VGroup(pr6, pr7, pr8, rnli, snli, cali, clli)
        ##########################################################################
        yn  = Dot([0,0,0],       color = BLUE)
        xn  = Dot([1,2,0],       color = BLUE)
        zn  = Dot([3,0,0],       color = BLUE)
        xnn = Dot([5/3, 4/3, 0], color = BLUE)
        xnyn = Line(xn, yn)
        xnzn = Line(xn, zn)
        ynzn = Line(zn, yn)
        xnny = Line(xnn, yn)
        xnl  = Tex(r"$x_n$", font_size = 34).next_to(xn , LEFT ).shift(0.2*UP   + 0.2*RIGHT)
        ynl  = Tex(r"$y$"  , font_size = 34).next_to(yn , LEFT ).shift(0.2*DOWN + 0.2*RIGHT)
        znl  = Tex(r"$z$"  , font_size = 34).next_to(zn , RIGHT).shift(0.2*DOWN + 0.2*LEFT )
        xnnl = Tex(r"$x_{n+1}$", font_size = 34).next_to(xnn , RIGHT).shift(0.2*UP + 0.2*LEFT )
        newf = VGroup(xn, yn, zn, xnyn, xnzn, ynzn, xnl, ynl, znl, xnn, xnnl, xnny)
        newf.shift(2*RIGHT + 0.55 *UP)
        ###########################################################
        self.add(still, newf)
        self.wait()
        self.play(Create(pr9), Create(pr99))
        self.wait()
        self.play(Create(pr10), Create(box))
        self.wait()



class topo1(Scene):
    def construct(self):
        topo0 = Tex(r"Theorem (Toponogov)", color = BLUE, font_size = 34).to_edge(UL).shift(RIGHT + DOWN)
        topo1 = Tex(r"Let $\Sigma$ be a complete surface such that for each $p \in \Sigma$,", font_size = 34).next_to(topo0, DOWN)
        topo2 = Tex(r"there is $r > 0 $ such that for any $a, b, c \in B(p, r)$,", font_size = 34).next_to(topo1, DOWN)
        topo3 = Tex(r"$\measuredangle (a_b^c) \geq \tilde{\measuredangle}(a_b^c)$.", font_size = 34).next_to(topo2, DOWN)
        topo4 = Tex(r"Then for each $x,y,z \in \Sigma $, one has", font_size = 34).next_to(topo3, DOWN)
        topo5 = Tex(r"$\measuredangle (x_y^z) \geq \tilde{\measuredangle}(x_y^z)$.", font_size = 34).next_to(topo4, DOWN)
        l = topo0.get_left()[0]
        topo1.shift((l - topo1.get_left()[0])*RIGHT)
        topo2.shift((l - topo2.get_left()[0])*RIGHT)
        topo3.shift( topo3.get_center()[0]*LEFT)
        topo4.shift((l - topo4.get_left()[0])*RIGHT)
        topo5.shift( topo5.get_center()[0]*LEFT)
        all = VGroup(topo0, topo1, topo2, topo3, topo4, topo5)
        self.play(Create(all))
        self.wait(2)
        pr0 = Tex(r"Proof", color = BLUE, font_size = 34).to_edge(UL).shift(0.5*RIGHT + 0.5*DOWN)
        pr1 = Tex(r"If the conclusion fails for a triangle ", r"$[x_1,y_1,z_1]$", r" with", font_size = 34).next_to(pr0, DOWN)
        pr2 = Tex(r"$\ell_1 = d(x_1,y_1) + d(x_1,z_1)$, ", font_size = 34).next_to(pr1, DOWN)
        pr3 = Tex(r"by the induction lemma, it also fails for a triangle ", r"$[x_2, y_2, z_2]$", r" with", font_size = 34).next_to(pr2, DOWN)
        pr4 = Tex(r"$ \ell_2 = d(x_2,y_2) + d(x_2,z_2)  < (2  / 3) \ell_1 $", font_size = 34).next_to(pr3, DOWN)
        pr5 = Tex(r"and $d (x_2, x_1) \leq 10 \ell_1 $.", r" By induction, the conlcusion fails for", font_size = 34).next_to(pr4, DOWN)
        pr6 = Tex(r"a sequence of triangles ", r"$[x_n, y_n, z_n]$", r" with ", font_size = 34).next_to(pr5, DOWN)
        pr7 = Tex(r"$ \ell_n = d(x_n,y_n) + d(x_n,z_n)  < (2 / 3) \ell_{n-1} $", font_size = 34).next_to(pr6, DOWN)
        pr8 = Tex(r"and $d(x_n, x_{n+1}) \leq  10 (2/3)^{n-1} \ell_1 $.", font_size = 34).next_to(pr7, DOWN)
        pr9 = Tex(r"The sequence $x_n$ converges, giving a contradiction.", font_size = 34).next_to(pr8, DOWN)
        l = pr0.get_left()[0]
        pr1.shift((l - pr1.get_left()[0])*RIGHT)
        pr3.shift((l - pr3.get_left()[0])*RIGHT)
        pr5.shift((l - pr5.get_left()[0])*RIGHT)
        pr6.shift((l - pr6.get_left()[0])*RIGHT)
        pr8.shift((l - pr8.get_left()[0])*RIGHT)
        pr9.shift((l - pr9.get_left()[0])*RIGHT)
        pr2.shift(( - pr2.get_center()[0] - 1)*RIGHT)
        pr4.shift(( - pr4.get_center()[0] - 1)*RIGHT)
        pr7.shift(( - pr7.get_center()[0] - 1)*RIGHT)
        blu = Line([0,0,0], [1.3,0,0], color = BLUE  ).next_to(pr1[1], DOWN).shift(0.15*UP)
        ora = Line([0,0,0], [1.3,0,0], color = ORANGE).next_to(pr3[1], DOWN).shift(0.15*UP)
        yel = Line([0,0,0], [1.3,0,0], color = YELLOW).next_to(pr6[1], DOWN).shift(0.15*UP)
        #####################################################
        rad = 0.05
        fo = 25
        a = 2
        x = Dot([0,0,0], color = RED, radius = rad)
        y = Dot([1/2,sqrt(3),0], color = RED, radius = rad)
        z = Dot([1,0,0], color = RED, radius = rad)
        xy = Line(x.get_center(), y.get_center(), stroke_width = a)
        xz = Line(x.get_center(), z.get_center(), stroke_width = a)
        yz = Line(z.get_center(), y.get_center(), stroke_width = a)
        xl = Tex(r"$x_1$", font_size = fo).next_to(x, DOWN).shift(0.3*UP + 0.2*LEFT)
        yl = Tex(r"$y_1$", font_size = fo).next_to(y, UP).shift(0.3*DOWN + 0.2*LEFT)
        zl = Tex(r"$z_1$", font_size = fo).next_to(z, DOWN).shift(0.3*UP + 0.2*RIGHT)
        x2 = Dot([0,0,0], color = RED, radius = rad)
        y2 = Dot([1/4,sqrt(3)/2,0], color = RED, radius = rad)
        z2 = Dot([1/2,0,0], color = RED, radius = rad)
        xy2 = Line(x2.get_center(), y2.get_center(), stroke_width = a)
        xz2 = Line(x2.get_center(), z2.get_center(), stroke_width = a)
        yz2 = Line(z2.get_center(), y2.get_center(), stroke_width = a)
        xl2 = Tex(r"$x_2$", font_size = fo).next_to(x2, DOWN).shift(0.3*UP + 0.2*LEFT)
        yl2 = Tex(r"$y_2$", font_size = fo).next_to(y2, UP).shift(0.3*DOWN + 0.2*LEFT)
        zl2 = Tex(r"$z_2$", font_size = fo).next_to(z2, DOWN).shift(0.3*UP + 0.2*RIGHT)
        x3 = Dot([0,0,0], color = RED, radius = rad)
        y3 = Dot([1/4,sqrt(3)/4,0], color = RED, radius = rad)
        z3 = Dot([1/4,0,0], color = RED, radius = rad)
        xy3 = Line(x3.get_center(), y3.get_center(), stroke_width = a)
        xz3 = Line(x3.get_center(), z3.get_center(), stroke_width = a)
        yz3 = Line(z3.get_center(), y3.get_center(), stroke_width = a)
        xl3 = Tex(r"$x_3$", font_size = fo).next_to(x3, DOWN).shift(0.3*UP + 0.2*LEFT)
        yl3 = Tex(r"$y_3$", font_size = fo).next_to(y3, UP).shift(0.3*DOWN + 0.2*LEFT)
        zl3 = Tex(r"$z_3$", font_size = fo).next_to(z3, DOWN).shift(0.3*UP + 0.2*RIGHT)
        g1 = VGroup(x,y,z,xy,xz,yz,xl,yl,zl)
        g2 = VGroup(x2,y2,z2,xy2,xz2,yz2,xl2,yl2,zl2)
        g3 = VGroup(x3,y3,z3,xy3,xz3,yz3,xl3,yl3,zl3)
        g1.shift(2.5*DOWN + 3*RIGHT)
        g2.shift(1.5*DOWN + 4*RIGHT)
        g3.shift(DOWN + 4.3*RIGHT)
        t1 = Polygon(x.get_center() , y.get_center() , z.get_center() , fill_opacity = 0.3 , fill_color = BLUE  , stroke_width = 0)
        t2 = Polygon(x2.get_center(), y2.get_center(), z2.get_center(), fill_opacity = 0.3 , fill_color = ORANGE, stroke_width = 0)
        t3 = Polygon(x3.get_center(), y3.get_center(), z3.get_center(), fill_opacity = 0.3 , fill_color = YELLOW, stroke_width = 0)
        #######################################################
        self.play(FadeOut(all), Create(pr0), Create(pr1), Create(pr2), Create(g1), Create(blu), Create(t1))
        self.wait()
        self.play(Create(pr3), Create(pr4), Create(pr5[0]), Create(g2), Create(ora), Create(t2))
        self.wait()
        self.play(Create(pr5[1]), Create(pr6), Create(pr7), Create(pr8), Create(g3), Create(yel), Create(t3))
        self.wait()
        self.play(Create(pr9))
        self.wait()

class kpl2(Scene):
    def construct(self):
        lem = Tex(r"Lemma", font_size = 34, color = BLUE).shift(5*LEFT +  UP)
        le1 = Tex(r"If $K \geq 0$, then for each $p \in \Sigma$ there is $r > 0$ such that", font_size = 34).next_to(lem, DOWN)
        le2 = Tex(r" if $x,y,z \in B(p, r)$, then ",  font_size = 34).next_to(le1, DOWN)
        le3 = Tex(r"$\measuredangle (x_y^z)$", r" $ \geq$ ", r"$\tilde{\measuredangle }(x_y^z)$", r".", font_size = 34).next_to(le2, DOWN)
        l = lem.get_left()[0]
        le1.shift((l - le1.get_left()[0])*RIGHT)
        le2.shift((l - le2.get_left()[0])*RIGHT)
        le3.shift(   le3.get_center()[0] *LEFT)
        old = VGroup(lem, le1, le2, le3)
        self.play(Create(old))
        self.wait(2)