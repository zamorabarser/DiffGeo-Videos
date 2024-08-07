from asyncio import threads
from this import d
from tkinter import E
from manim import *
from numpy import sqrt
import math



class ti(Scene):
    def construct(self):
        t1 = Text("Flat surfaces and", font_size=60).shift(0.5*UP)
        t2 = Text("comparison", font_size=60).shift(0.5*DOWN)
        self.play(Write(t1), Write(t2))
        self.wait()        


class fst(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/16)
        t = Tex(r"Theorem", font_size = 34, color = BLUE).to_edge(UL).shift(0.6*DOWN+RIGHT)
        t1 = Tex(r"If ", r"$\Sigma$", r"  $\subset \mathbb{R}^3$ is a proper surface with zero Gauss curvature, then", font_size = 34).next_to(t, DOWN)
        t2 = Tex(r"it is a cylinder. That is, there is a line ",r"$L$", r" such that for all ", r"$p$", r" $\in \Sigma$,", font_size = 34).next_to(t1, DOWN)
        t3 = Tex(r"the line parallel to $L$ passing through $p$ lies entirely in $\Sigma$.", font_size = 34).next_to(t2, DOWN)
        l = t.get_left()[0]
        t1.shift([t1.get_left()[0]- l ] *LEFT)
        t2.shift([t2.get_left()[0]- l ] *LEFT)
        t3.shift([t3.get_left()[0]- l ] *LEFT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(t1[1], DOWN).shift(0.15*UP)
        lli = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(t2[1], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(t2[3], DOWN).shift(0.15*UP)
        op = 0.2
        h = - 1.5
        hei = 1.5
        w = 1/2
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift( [ 0, - w * PI / 3 - 2 , hei * sqrt(3) /2 + h ] )
        p2 = p.copy()
        line = ParametricFunction(lambda t : [ - t , 0.5 * t + 4, h + 1 ] , t_range = [-4, 4], color = BLUE)
        line2 = ParametricFunction(lambda t : [ - t , 0.5 * t - w * PI / 3  -2 , hei * sqrt(3) /2  + h ] , t_range = [-3, 3], color = BLUE)
        lines = []
        linegroup = VGroup()
        for i in range(0,13):
            lines.append( ParametricFunction(lambda t: [ - t , 0.5 * t + w*PI*(-1 + i/6)  -2 ,  - hei * np.sin (PI*(-1 + i/6)) + h ] , t_range = [-3,3], color = BLUE ) )
            linegroup.add(lines[i])
        s = Surface(
            lambda u, v:  [ - v , 0.5 * v + w*u  -2 ,  - hei * np.sin (u) + h ],
            u_range = [ -PI , PI ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            #resolution = [ 32, 16]
        )
        self.add_fixed_in_frame_mobjects(t, t1, t2, t3, sli, lli, pli)
        self.remove(t, t1, t2, t3, sli, lli, pli)
        self.play(Create(t), Create(t1), Create(t2), Create(t3), Create(sli), Create(lli), Create(pli), Create(s))
        self.wait()
        self.play(Create(p), Create(line))
        self.wait()
        self.play(Create(line2))
        self.add(p2)
        self.wait()
        self.play(Create(linegroup))
        self.add(p2)
        self.wait()


class adef(Scene):
    def construct(self):
        t = Tex(r"Theorem", font_size = 34, color = BLUE).to_edge(UL).shift(0.6*DOWN+RIGHT)
        t1 = Tex(r"If ", r"$\Sigma$", r" $\subset \mathbb{R}^3$ is a proper surface with zero Gauss curvature, then", font_size = 34).next_to(t, DOWN)
        t2 = Tex(r"it is a cylinder. That is, there is a line ",r"$L$", r" such that for all ", r"$p$", r" $\in \Sigma$,", font_size = 34).next_to(t1, DOWN)
        t3 = Tex(r"the line parallel to $L$ passing through $p$ lies entirely in $\Sigma$.", font_size = 34).next_to(t2, DOWN)
        l = t.get_left()[0]
        t1.shift([t1.get_left()[0]- l ] *LEFT)
        t2.shift([t2.get_left()[0]- l ] *LEFT)
        t3.shift([t3.get_left()[0]- l ] *LEFT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(t1[1], DOWN).shift(0.15*UP)
        lli = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(t2[1], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(t2[3], DOWN).shift(0.15*UP)
        a1 = Tex(r"$A : = \{ x \in \Sigma \, \vert \, S_x \neq 0 \} $", font_size = 34).next_to(t3, DOWN).shift(0.2*DOWN)
        a2 = Tex(r"$S_x : T_x \Sigma \to T_x \Sigma $ is diagonalizable with orthogonal eigenspaces and", font_size = 34).next_to(a1, DOWN).shift(0.2*DOWN)
        a3 = Tex(r"det$(S_x) = K (x) = 0$.", font_size = 34).next_to(a2, DOWN)
        a1.shift([a1.get_center()[0] ] *LEFT)
        a2.shift([a2.get_left()[0]- l ] *LEFT)
        a3.shift([a3.get_center()[0] ] *LEFT)
        ali = Line([0,0,0], [3,0,0], color = YELLOW).next_to(a1, DOWN).shift(0.1*UP)
        detli = Line([0,0,0], [3,0,0], color = YELLOW).next_to(a3, DOWN).shift(0.1*UP)
        self.add(t, t1, t2, t3, sli, lli, pli)
        self.wait()
        self.play(Write(a1), Create(ali))
        self.wait()
        self.play(Write(a2), Write(a3), Create(detli))
        self.wait()


class seg(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta = 0 )
        t = Tex(r"Lemma", font_size = 34, color = BLUE).to_edge(UL).shift(0.3*DOWN+ 0.6*RIGHT)
        t1 = Tex(r"If ", r"$\Sigma$", r"  $\subset \mathbb{R}^3$ is a surface with zero Gauss curvature, then for each", font_size = 34).next_to(t, DOWN)
        t2 = Tex(r"$p$", r" $\in A$, there is a line segment ", r"$L$", r" passing through $p$ contained in $A$.", font_size = 34).next_to(t1, DOWN)
        t3 = Tex(r"Moreover, if ", r"$L^{\prime}$", r" is another segment passing through $p$ and contained in $\Sigma$,", font_size = 34).next_to(t2, DOWN)
        t4 = Tex(r"the ", r"line containing $L$", r" contains also $L^{\prime}$.", font_size = 34).next_to(t3, DOWN)
        l = t.get_left()[0]
        t1.shift([t1.get_left()[0]- l ] *LEFT)
        t2.shift([t2.get_left()[0]- l ] *LEFT)
        t3.shift([t3.get_left()[0]- l ] *LEFT)
        t4.shift([t4.get_left()[0]- l ] *LEFT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(t1[1], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(t2[0], DOWN).shift(0.15*UP)
        lli = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(t2[2], DOWN).shift(0.15*UP)
        op = 0.2
        th = - PI /6
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift( [ 2 * np.cos( th )  ,  2 * np.sin( th ) ,  -1 ] )
        line = ParametricFunction(lambda t : [ (2 - t ) * np.cos(th) , (2 - t) * np.sin(th) , t  -1 ] , t_range = [-0.5, 0.5], color = BLUE)
        s = Surface(
            lambda u, v:  [ u * np.cos(v)  ,  u * np.sin(v) ,  - u + 1   ],
            u_range = [  1 , 3  ],
            v_range = [ -PI /2 , PI / 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        self.add_fixed_in_frame_mobjects(t, t1, t2, t3, t4, sli, lli, pli)
        self.remove(t, t1, t2, t3, t4, sli, lli, pli)
        self.play(Create(t), Create(t1), Create(t2), Create(t3), Create(t4), Create(sli), Create(lli), Create(pli), Create(s), Create(p))
        self.wait()
        self.play(Create(line))
        self.wait()
        self.play(FadeOut(line))
        w =  3
        self.play(s.animate.shift([0,w,0]), p.animate.shift([0,w,0]))
        self.wait()
        



class segp(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta = 0 )
        t = Tex(r"Lemma", font_size = 34, color = BLUE).to_edge(UL).shift(0.3*DOWN+ 0.6*RIGHT)
        t1 = Tex(r"If ", r"$\Sigma$", r"  $\subset \mathbb{R}^3$ is a surface with zero Gauss curvature, then for each", font_size = 34).next_to(t, DOWN)
        t2 = Tex(r"$p$", r" $\in A$, there is a line segment ", r"$L$", r" passing through $p$ contained in $A$.", font_size = 34).next_to(t1, DOWN)
        t3 = Tex(r"Moreover, if ", r"$L^{\prime}$", r" is another segment passing through $p$ and contained in $\Sigma$,", font_size = 34).next_to(t2, DOWN)
        t4 = Tex(r"the ", r"line containing $L$", r" contains also $L^{\prime}$.", font_size = 34).next_to(t3, DOWN)
        l = t.get_left()[0]
        t1.shift([t1.get_left()[0]- l ] *LEFT)
        t2.shift([t2.get_left()[0]- l ] *LEFT)
        t3.shift([t3.get_left()[0]- l ] *LEFT)
        t4.shift([t4.get_left()[0]- l ] *LEFT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(t1[1], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(t2[0], DOWN).shift(0.15*UP)
        lli = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(t2[2], DOWN).shift(0.15*UP)
        op = 0.2
        th = - PI /6
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift( [ 2 * np.cos( th )  ,  2 * np.sin( th ) ,  -1 ] )
        s = Surface(
            lambda u, v:  [ u * np.cos(v)  ,  u * np.sin(v) ,  - u + 1   ],
            u_range = [  1 , 3  ],
            v_range = [ -PI /2 , PI / 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        w = 3
        s.shift([0,w,0])
        p.shift([0,w,0])
        self.add_fixed_in_frame_mobjects(t, t1, t2, t3, t4, sli, lli, pli)
        self.add(s, p)
        self.wait()
        ch1 = Tex(r"Eigenspaces of $S_x$ give a", font_size = 34).next_to(t4, DOWN).shift(DOWN)
        ch2 = Tex(r"decomposition $T_x \Sigma = $ ", r"$E_1$", r" $\oplus$ ", r"$E_2$", font_size = 34).next_to(ch1, DOWN)
        ch3 = Tex(r"with $E_1 = $ Ker$(S_x)$", font_size = 34).next_to(ch2, DOWN)
        ch1.shift([ch1.get_left()[0]- l - 0.5 ] *LEFT)
        ch2.shift([ch2.get_left()[0]- l - 0.5 ] *LEFT)
        ch3.shift([ch3.get_left()[0]- l - 0.5 ] *LEFT)
        e1li = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(ch2[1], DOWN).shift(0.15*UP)
        e2li = Line([0,0,0], [0.4,0,0], color = RED).next_to(ch2[3], DOWN).shift(0.15*UP)
        e1 = []
        e2 = []
        e12 = VGroup()
        for i in range(0,5):
            e1.append([])
            e2.append([])
            for j in range(0,5):
                e1[i].append(  ParametricFunction( lambda t : [ (4/3 + j /3 + t ) * np.cos(-5*PI/12 + i * PI / 8 )  , (4/3 + j /3 + t ) * np.sin(-5*PI/12 + i * PI / 8) + w ,  - (4/3 + j /3 + t ) + 1   ] , t_range = [-0.1 , 0.1 ] , color = YELLOW )   )
                e2[i].append(  ParametricFunction( lambda t : [ (4/3 + j /3 ) * np.cos(-5*PI/12 + i * PI / 8)  + t * np.sin(-5*PI/12 + i * PI / 8)  , (4/3 + j /3  ) * np.sin(-5*PI/12 + i * PI / 8) + w - t * np.cos(-5*PI/12 + i * PI / 8) ,  - (4/3 + j /3 ) + 1   ] , t_range = [-0.15 , 0.15 ] , color = RED )   )
                e12.add(e1[i][j], e2[i][j])
        self.add_fixed_in_frame_mobjects(ch1, ch2, ch3, e1li, e2li)
        self.remove(ch1, ch2, ch3, e1li, e2li)
        self.play(Create(ch1), Write(ch2), Write(ch3), Create(e1li), Create(e2li))
        self.wait()
        self.play(Create(e12))
        self.wait()
        block1 = VGroup(t, t1, t2, t3, t4, sli, pli, lli) 
        gap1 = -0.3
        gap2 = gap1 + ch1.get_right()[0] - ch2.get_left()[0] + 0.2 
        gap3 = gap2 + ch2.get_right()[0] - ch3.get_left()[0] + 0.2
        vgap1 = t1.get_center()[1] - ch1.get_center()[1]
        vgap2 = vgap1 + ch1.get_center()[1] - ch2.get_center()[1]
        vgap3 = vgap2 + ch2.get_center()[1] - ch3.get_center()[1]
        self.play( FadeOut(block1), ch1.animate.shift( [ gap1 , vgap1 , 0 ] ), ch2.animate.shift( [ gap2 , vgap2 , 0 ] ), e1li.animate.shift( [ gap2 , vgap2 , 0 ] ), e2li.animate.shift( [ gap2 , vgap2 , 0 ] ), ch3.animate.shift( [ gap3 , vgap3 , 0 ] ) )
        self.wait()
        char1 = Tex(r"There is a chart $s : $ ", r"$U$", r" $ \to \Sigma$ with ", r"$s(0,0) = p$", r", $\vert s_u \vert = 1 $ along $\{v = 0\}$,", font_size = 34).next_to(t1, DOWN)
        char2 = Tex(r"$S(s_u) = 0$", r", and  ", r"$ S (s_v) = - k (s_v) $", r".", font_size = 34).next_to(char1, DOWN).shift(0.1*DOWN)
        lnew = ch1.get_left()[0]
        char1.shift([char1.get_left()[0]- lnew ]*LEFT)
        char2.shift([char2.get_center()[0] ]*LEFT)
        uli = Line([0,0,0], [0.4,0,0], color = GREEN).next_to(char1[1], DOWN).shift(0.15*UP)
        p0li = Line([0,0,0], [1.3,0,0], color = ORANGE).next_to(char1[3], DOWN).shift(0.15*UP)
        suli = Line([0,0,0], [1.5,0,0], color = YELLOW).next_to(char2[0], DOWN).shift(0.15*UP)
        svli = Line([0,0,0], [2.3,0,0], color = RED).next_to(char2[2], DOWN).shift(0.15*UP)
        u = Surface(
            lambda u, v:  [ u* sqrt(3) /2 + v /2 , - u  /2 +  v* sqrt(3)/2 - 3.3 , -1.5 ],
            u_range = [  -2 , 2  ],
            v_range = [ -2 , 2 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = op, 
            resolution = [ 10 , 10 ]
        )
        ucli = []
        vcli = []
        cli = VGroup()
        for i in range(0,5):
            ucli.append( Line( [ -2 * sqrt(3) /2 + (-1.5 + 0.75 * i ) /2 ,  2  /2 + ( -1.5 + 0.75 * i )* sqrt(3)/2 - 3.3 , -1.5]  , [ 2* sqrt(3) /2 + (-1.5 + 0.75 * i ) /2 , - 2  /2 +  (-1.5 + 0.75 * i )* sqrt(3)/2 - 3.3 , -1.5]  , color = YELLOW ))
            vcli.append( Line( [ (-1.5 + 0.75 * i)* sqrt(3) /2  -2 /2 , - (-1.5 + 0.75 * i)  /2 -  2 * sqrt(3)/2 - 3.3 , -1.5 ] , [ (-1.5 + 0.75 * i)* sqrt(3) /2 + 2 /2 , - (-1.5 + 0.75 * i) /2 + 2 * sqrt(3)/2 - 3.3 , -1.5]   , color = RED))
            cli.add(ucli[i], vcli[i])
        sar = CurvedArrow([-1.5,-0.5,0], [0.5,-0.5,0], radius = -4)
        sarli = Tex(r"$s$", font_size = 34).next_to(sar, UP)
        p2 = Sphere(radius = 0.05)
        p2.set_color(ORANGE).shift( [ 0 , -3.3 , -1.5 ] )
        self.add_fixed_in_frame_mobjects(sar, sarli, char1, char2, suli, svli, uli, p0li)
        self.remove(sar, sarli, char1, char2, suli, svli, uli, p0li)
        self.play(Write(char1), Write(char2), Create(u), Create(sar), Create(sarli), Create(uli), Create(p2), Create(p0li))
        self.wait()
        self.play(Create(cli), Create(suli), Create(svli))
        self.wait()
        self.play(FadeOut(cli), FadeOut(u), FadeOut(sar), FadeOut(sarli), FadeOut(p2), FadeOut(e12))
        self.wait()




class segp2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta = 0 )
        t = Tex(r"Lemma", font_size = 34, color = BLUE).to_edge(UL).shift(0.3*DOWN+ 0.8*RIGHT)
        t1 = Tex(r"If ", r"$\Sigma$", r"  $\subset \mathbb{R}^3$ is a surface with zero Gauss curvature, then for each", font_size = 34).next_to(t, DOWN)
        l = t.get_left()[0]
        op = 0.2
        th = - PI /6
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift( [ 2 * np.cos( th )  ,  2 * np.sin( th ) + 3 ,  -1 ] )
        s = Surface(
            lambda u, v:  [ u * np.cos(v)  ,  u * np.sin(v) + 3 ,  - u + 1   ],
            u_range = [  1 , 3  ],
            v_range = [ -PI /2 , PI / 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        ch1 = Tex(r"Eigenspaces of $S_x$ give a", font_size = 34)
        ch2 = Tex(r"decomposition $T_x \Sigma = $ ", r"$E_1$", r" $\oplus$ ", r"$E_2$", font_size = 34)
        ch3 = Tex(r"with $E_1 = $ Ker$(S_x)$", font_size = 34)
        gap1 = l - ch1.get_left()[0]
        gap2 = gap1 + ch1.get_right()[0] - ch2.get_left()[0] + 0.2 
        gap3 = gap2 + ch2.get_right()[0] - ch3.get_left()[0] + 0.2
        vgap = t1.get_center()[1] - ch1.get_center()[1]
        char1 = Tex(r"There is a chart $s : $ ", r"$U$", r" $ \to \Sigma$ with ", r"$s(0,0) = p$", r", $\vert s_u \vert = 1 $ along $\{v = 0\}$,", font_size = 34).next_to(t1, DOWN)
        char2 = Tex(r"$S(s_u) = 0$", r", and  ", r"$ S (s_v) = - k (s_v) $", r".", font_size = 34).next_to(char1, DOWN).shift(0.1*DOWN)
        ch1.shift([gap1, vgap, 0])
        ch2.shift([gap2, vgap, 0])
        ch3.shift([gap3, vgap, 0])
        char1.shift([char1.get_left()[0]- l ]*LEFT)
        char2.shift([char2.get_center()[0] ]*LEFT)
        e1li = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(ch2[1], DOWN).shift(0.15*UP)
        e2li = Line([0,0,0], [0.4,0,0], color = RED).next_to(ch2[3], DOWN).shift(0.15*UP)
        uli = Line([0,0,0], [0.4,0,0], color = GREEN).next_to(char1[1], DOWN).shift(0.15*UP)
        p0li = Line([0,0,0], [1.3,0,0], color = ORANGE).next_to(char1[3], DOWN).shift(0.15*UP)
        suli = Line([0,0,0], [1.5,0,0], color = YELLOW).next_to(char2[0], DOWN).shift(0.15*UP)
        svli = Line([0,0,0], [2.3,0,0], color = RED).next_to(char2[2], DOWN).shift(0.15*UP)
        u = Surface(
            lambda u, v:  [ u* sqrt(3) /2 + v /2 , - u  /2 +  v* sqrt(3)/2 - 3.3 , -1.5 ],
            u_range = [  -2 , 2  ],
            v_range = [ -2 , 2 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = op, 
            resolution = [ 10 , 10 ]
        )
        xvf = []
        yvf = []
        nvf = []
        frame = VGroup()
        for i in range(0,5):
            xvf.append([])
            yvf.append([])
            nvf.append([])
            for j in range(0,5):
                xvf[i].append(  
                    Arrow(  
                        [ (4/3 +  j /3 + 0.165 ) * np.cos(-5*PI/12 + i * PI / 8 )  , (4/3 +  j /3  + 0.165 ) * np.sin(-5*PI/12 + i * PI / 8) + 3 ,  - (4/3 +  j /3 + 0.165 ) + 1  ] , 
                        [ (4/3 +  j /3 - 0.5 ) * np.cos(-5*PI/12 + i * PI / 8 )  , (4/3 +  j /3 - 0.5  ) * np.sin(-5*PI/12 + i * PI / 8) + 3 ,  - (4/3 +  j /3 - 0.5) + 1  ] , 
                        color = YELLOW, 
                        max_tip_length_to_length_ratio = 0.2
                    ) 
                )
                yvf[i].append(  
                    Arrow(  
                        [ (4/3 + j /3 ) * np.cos(-5*PI/12 + i * PI / 8 ) + (0.25 )*np.sin(-5*PI/12 + i * PI / 8) , (4/3 + j /3  ) * np.sin(-5*PI/12 + i * PI / 8) + 3 - (0.25 )*np.cos(-5*PI/12 + i * PI / 8),  - (4/3 +  j /3 ) + 1  ] , 
                        [ (4/3 + j /3 ) * np.cos(-5*PI/12 + i * PI / 8 ) - (0.75 )*np.sin(-5*PI/12 + i * PI / 8) , (4/3 + j /3  ) * np.sin(-5*PI/12 + i * PI / 8) + 3 + (0.75 )*np.cos(-5*PI/12 + i * PI / 8),  - (4/3 +  j /3 ) + 1  ] , 
                        color = RED, 
                        max_tip_length_to_length_ratio = 0.2
                    ) 
                )
                nvf[i].append(  
                    Arrow(  
                        [ (4/3 +  j /3 ) * np.cos(-5*PI/12 + i * PI / 8 ) - 0.165* (4/3 +  j /3 ) * np.cos(-5*PI/12 + i * PI / 8 ) / sqrt( 1 + (4/3 +  j /3 )**2 ) , (4/3 +  j /3 ) * np.sin(-5*PI/12 + i * PI / 8) + 3 - 0.165*(4/3 +  j /3 ) * np.sin(-5*PI/12 + i * PI / 8 ) / sqrt( 1 + (4/3 +  j /3 )**2 ),  - (4/3 +  j /3 ) + 1 - 0.165* (4/3 +  j /3 ) / sqrt( 1 + (4/3 +  j /3 )**2 ) ] , 
                        [ (4/3 +  j /3 ) * np.cos(-5*PI/12 + i * PI / 8 ) + 0.7  * (4/3 +  j /3 ) * np.cos(-5*PI/12 + i * PI / 8 ) / sqrt( 1 + (4/3 +  j /3 )**2 ) , (4/3 +  j /3 ) * np.sin(-5*PI/12 + i * PI / 8) + 3 + 0.7 * (4/3 +  j /3 ) * np.sin(-5*PI/12 + i * PI / 8 ) / sqrt( 1 + (4/3 +  j /3 )**2 ),  - (4/3 +  j /3 ) + 1  +  0.7* (4/3 +  j /3 ) / sqrt( 1 + (4/3 +  j /3 )**2 ) ] , 
                        color = PINK, 
                        max_tip_length_to_length_ratio = 0.2
                    ) 
                )
                frame.add(xvf[i][j], yvf[i][j], nvf[i][j])
        x = Tex(r"$X : = s_u / \vert s_u \vert$", r", ", r"$Y : = s_v / \vert s_v \vert$", r", ", r"$N : = X \times Y $",font_size = 34).next_to(char2, DOWN).shift(0.2*DOWN)
        nu = Tex(r"$N_u = - S(s_u) = 0 $", r", $ \, $ $N_v = - S(s_v) = -k s_v$." , font_size = 34).next_to(x, DOWN).shift(0.1*DOWN)
        x.shift((x.get_left()[0] - l )*LEFT)
        nu.shift((l - nu.get_left()[0])*RIGHT)
        claim = Tex(r"Claim", r": $X$, $Y$, $N$ don't depend on $u$." , font_size = 34).next_to(nu, DOWN)
        claim.shift((l - claim.get_left()[0])*RIGHT)
        xli = Line([0,0,0], [1.5,0,0], color = YELLOW).next_to(x[0], DOWN).shift(0.15*UP)
        yli = Line([0,0,0], [1.5,0,0], color = RED).next_to(x[2], DOWN).shift(0.15*UP)
        nli = Line([0,0,0], [1.5,0,0], color = PINK).next_to(x[4], DOWN).shift(0.15*UP)
        cli = Line([0,0,0], [0.7,0,0]).next_to(claim[0], DOWN).shift(0.15*UP)
        self.add_fixed_in_frame_mobjects(ch1, ch2, ch3, char1, char2, e1li, e2li, uli, suli, svli, p0li, x, nu, xli, yli, nli, claim, cli)
        self.add(s, p)
        self.remove(x, nu, xli, yli, nli, claim, cli)
        self.wait()
        self.play(Create(frame), Create(x), Create(xli), Create(yli), Create(nli))
        self.wait()
        self.play(Create(nu))
        self.wait()
        self.play(Create(claim), Create(cli))
        self.wait()
        self.play(FadeOut(s), FadeOut(frame), FadeOut(p), claim.animate.shift( claim.get_center()[0]*LEFT ), cli.animate.shift( claim.get_center()[0]*LEFT ) )
        self.wait()
        


class segp3(Scene):
    def construct(self):
        t = Tex(r"Lemma", font_size = 34, color = BLUE).to_edge(UL).shift(0.3*DOWN+ 0.8*RIGHT)
        t1 = Tex(r"If ", r"$\Sigma$", r"  $\subset \mathbb{R}^3$ is a surface with zero Gauss curvature, then for each", font_size = 34).next_to(t, DOWN)
        l = t.get_left()[0]
        ch1 = Tex(r"Eigenspaces of $S_x$ give a", font_size = 34)
        ch2 = Tex(r"decomposition $T_x \Sigma = $ ", r"$E_1$", r" $\oplus$ ", r"$E_2$", font_size = 34)
        ch3 = Tex(r"with $E_1 = $ Ker$(S_x)$", font_size = 34)
        gap1 = l - ch1.get_left()[0]
        gap2 = gap1 + ch1.get_right()[0] - ch2.get_left()[0] + 0.2 
        gap3 = gap2 + ch2.get_right()[0] - ch3.get_left()[0] + 0.2
        vgap = t1.get_center()[1] - ch1.get_center()[1]
        char1 = Tex(r"There is a chart $s : $ ", r"$U$", r" $ \to \Sigma$ with ", r"$s(0,0) = p$", r", $\vert s_u \vert = 1 $ along $\{ v = 0\} $,", font_size = 34).next_to(t1, DOWN)
        char2 = Tex(r"$S(s_u) = 0$", r", and  ", r"$ S (s_v) = - k (s_v) $", r".", font_size = 34).next_to(char1, DOWN).shift(0.1*DOWN)
        ch1.shift([gap1, vgap, 0])
        ch2.shift([gap2, vgap, 0])
        ch3.shift([gap3, vgap, 0])
        char1.shift([char1.get_left()[0]- l ]*LEFT)
        char2.shift([char2.get_center()[0] ]*LEFT)
        e1li = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(ch2[1], DOWN).shift(0.15*UP)
        e2li = Line([0,0,0], [0.4,0,0], color = RED).next_to(ch2[3], DOWN).shift(0.15*UP)
        uli = Line([0,0,0], [0.4,0,0], color = GREEN).next_to(char1[1], DOWN).shift(0.15*UP)
        p0li = Line([0,0,0], [1.3,0,0], color = ORANGE).next_to(char1[3], DOWN).shift(0.15*UP)
        suli = Line([0,0,0], [1.5,0,0], color = YELLOW).next_to(char2[0], DOWN).shift(0.15*UP)
        svli = Line([0,0,0], [2.3,0,0], color = RED).next_to(char2[2], DOWN).shift(0.15*UP)
        x = Tex(r"$X : = s_u / \vert s_u \vert$", r", ", r"$Y : = s_v / \vert s_v \vert$", r", ", r"$N : = X \times Y $",font_size = 34).next_to(char2, DOWN).shift(0.2*DOWN)
        nu = Tex(r"$N_u = - S(s_u) = 0 $", r", $ \, $ $N_v = - S(s_v) = -k s_v$." , font_size = 34).next_to(x, DOWN).shift(0.1*DOWN)
        x.shift((x.get_left()[0] - l )*LEFT)
        nu.shift((l - nu.get_left()[0])*RIGHT)
        claim = Tex(r"Claim", r": $X$, $Y$, $N$ don't depend on $u$." , font_size = 34).next_to(nu, DOWN)
        claim.shift((claim.get_center()[0])*LEFT)
        yu = Tex(r"$0 = N_{uv} $", r" $= (N_v)_u$", r" $ = (-ks_v)_u$", r" $ = (-k\vert s_v \vert Y )_u$", r" $ = (-k \vert s_v \vert )_u Y - k\vert s_v \vert Y_u$." , font_size = 34).next_to(claim, DOWN)
        yu.shift((l - yu.get_left()[0])*RIGHT)
        yu2 = Tex(r"$ \vert Y \vert  \equiv 1$ $\Rightarrow$ $Y \perp Y_u $ ", r" $\Rightarrow $ ", r" $ k \vert s_v \vert Y_u = 0 $ ", r" $\Rightarrow$ ", r"$Y_u = 0$.", font_size = 34).next_to(yu, DOWN)
        yu2.shift((l - yu2.get_left()[0])*RIGHT)
        xu = Tex(r"$ X_u = (Y \times N)_u = $ ", r"$Y_u$", r" $ \times N + Y \times  $ ", r"$N_u$", r" $ =0 $.", font_size = 34).next_to(yu2, DOWN)
        xu.shift((l - xu.get_left()[0])*RIGHT)
        xli = Line([0,0,0], [1.5,0,0], color = YELLOW).next_to(x[0], DOWN).shift(0.15*UP)
        yli = Line([0,0,0], [1.5,0,0], color = RED).next_to(x[2], DOWN).shift(0.15*UP)
        nli = Line([0,0,0], [1.5,0,0], color = PINK).next_to(x[4], DOWN).shift(0.15*UP)
        cli = Line([0,0,0], [0.7,0,0]).next_to(claim[0], DOWN).shift(0.15*UP)
        yuzli = Line([0,0,0], [0.4,0.2,0], color = RED).next_to(xu[1], DOWN).shift(0.5*UP)
        nuzli = Line([0,0,0], [0.4,0.2,0], color = RED).next_to(xu[3], DOWN).shift(0.5*UP)
        nubox = SurroundingRectangle(nu[0], color = BLUE, buff = 0.1)
        yubox = SurroundingRectangle(yu2[4], color = BLUE, buff = 0.1)
        xubox = SurroundingRectangle(xu, color = BLUE, buff = 0.1)
        self.add(ch1, ch2, ch3, char1, char2, e1li, e2li, uli, suli, svli, p0li, x, nu, xli, yli, nli, claim, cli)
        self.wait()
        self.play(Write(yu[0]))
        self.wait()
        self.play(Write(yu[1]))
        self.wait()
        self.play(Write(yu[2]))
        self.wait()
        self.play(Write(yu[3]))
        self.wait()
        self.play(Write(yu[4]))
        self.wait()
        self.play(Write(yu2[0]))
        self.wait()
        self.play(Write(yu2[1]), Write(yu2[2]))
        self.wait()
        self.play(Write(yu2[3]), Write(yu2[4]))
        self.wait()
        self.play(Write(xu[0]))
        self.wait()
        self.play(Write(xu[1]), Write(xu[2]), Write(xu[3]))
        self.wait()
        self.play(Create(yuzli), Create(nuzli))
        self.wait()
        self.play(Write(xu[4]))
        self.wait()
        self.play(Create(nubox), Create(yubox), Create(xubox))
        self.wait()
        


        
class segp4(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/16)
        op = 0.2
        xsu = Tex(r"$X$ $=$ direction of curve ", r"$( u \mapsto s(u,0) )$", r", so this curve is a line.", font_size = 34).to_edge(UL).shift(0.6*DOWN+ 0.3*RIGHT)
        gam = Tex(r"Assume ", r"$\gamma $", r" $: ( - \varepsilon, \varepsilon ) \to \Sigma $ is another unit-speed segment with $\gamma (0 ) = p$.", font_size = 34).next_to(xsu, DOWN)
        gamp0 = Tex(r"$\gamma ^{\prime} (0) = \alpha X + \beta Y$", font_size = 34).next_to(gam, DOWN)
        gam1 = Tex(r"$S(\gamma ^{\prime } (0)) \cdot \gamma ^{\prime }(0) $", r" $ =  S(\alpha X  + \beta Y ) \cdot (\alpha X  + \beta Y)  $", r" $ = ( k \beta Y) (\alpha X  + \beta Y) $", r" $= k \beta ^2 $.", font_size = 34).next_to(gamp0, DOWN)
        gam2 = Tex(r"$ S ( \gamma ^{\prime} (0)) \cdot \gamma ^{\prime }(0)  =  \gamma ^{\prime \prime }(0) \cdot N  $",  r" $ = 0$.",  font_size = 34).next_to(gam1, DOWN)
        gam3 = Tex(r"$ \Rightarrow \, $ $ \beta = 0$.",  font_size = 34).next_to(gam2, DOWN)
        l = xsu.get_left()[0]        
        gam.shift((gam.get_left()[0] - l)*LEFT)
        gamp0.shift((gamp0.get_center()[0])*LEFT)
        gam1.shift((gam1.get_left()[0] - l)*LEFT)
        gam2.shift((gam2.get_left()[0] - l)*LEFT)
        gam3.shift((gam3.get_center()[0])*LEFT)
        linli = Line([0,0,0], [1.3,0,0], color = BLUE).next_to(xsu[1], DOWN).shift(0.15*UP)
        gli = Line([0,0,0], [0.4,0,0], color = GREEN).next_to(gam[1], DOWN).shift(0.15*UP)
        h = 1.8
        s = Surface(
            lambda u, v:  [ u , v , ((u + v)**2 - ( u - v )**2 ) /40 - h   ],
            u_range = [ -3 , 3 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift( [ 0 , 0 , - h ] )
        line = ParametricFunction(lambda t : [ t ,  0  , - h ] , t_range = [-1, 1], color = BLUE)
        g    = ParametricFunction(lambda t : [ 0 ,  t , - h ] ,   t_range = [-1, 1], color = GREEN)
        self.add_fixed_in_frame_mobjects(xsu, gam, gamp0, gam1, gli, linli, gam2, gam3)
        self.remove(xsu, gam, gamp0, gam1, gli, linli, gam2, gam3)
        self.play(Write(xsu), Create(s), Create(p), Create(line), Create(linli))
        self.wait()
        self.play(Write(gam), Create(g), Create(gli))
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/4)
        self.wait(2)
        self.stop_ambient_camera_rotation()
        self.wait()
        self.play(Write(gamp0))
        self.wait()
        self.play(Write(gam1[0]))
        self.wait()
        self.play(Write(gam1[1]))
        self.wait()
        self.play(Write(gam1[2]))
        self.wait()
        self.play(Write(gam1[3]))
        self.wait()
        self.play(Write(gam2[0]))
        self.wait()
        self.play(Write(gam2[1]))
        self.wait()
        self.play(Write(gam3))
        self.wait()




        
class keq(Scene):
    def construct(self):
        le = Tex(r"Lemma", font_size = 34, color = BLUE).to_edge(UL).shift( 0.1 * DOWN + 0.1 * RIGHT ) 
        le1 = Tex(r"With $s: U \to \Sigma $ defined as before, one has", font_size = 34).next_to(le, DOWN)
        le2 = Tex(r"$ \partial_u ^2  \left[ \frac{1}{k} \right] = 0 .  $", font_size = 34).next_to(le1, DOWN)
        pr = Tex(r"Proof", font_size = 34, color = BLUE).next_to(le2, DOWN)
        pr1 = Tex(r"$\partial_u^2 \left[ \frac{1}{k} \right] =  \partial_u \left[ - \frac{ 1}{k^2} k_u \right]  $", r" $ = \, \frac{1}{k^3} $", r"$\, [ 2 (k_u)^2 - k k _{uu} ]  $", font_size = 34).next_to(pr, DOWN)
        pr2 = Tex(r"$ s_{uu} = (s_u)_u $", r" $=  (\vert s_u \vert X) _u $", r" $ = \vert s_u \vert_u X + \vert s_u \vert$", r"$\, X_u$", r" $ = \vert s_u \vert_u X  $ ", r" $ \Rightarrow $ ", r"$ s_{uu} \parallel s_u \perp s_v$", font_size = 34).next_to(pr1, DOWN)
        pr3 = Tex(r"$ \partial _v \vert s_u \vert ^2 = 2 s_u \cdot s_{uv} $", r" $= 2 [  ( $ ", r"$s_u \cdot s_v$", r" $ )_u  - $ ", r"$ s_{uu} \cdot s_ v$", r" $ ]$", r" $  = 0$", font_size = 34).next_to(pr2, DOWN)
        pr4 = Tex(r"$\vert s_u \vert $ is constant along $\{ v= 0 \}$  $ \Rightarrow $  $\vert s_u \vert $ is constant over $U$", font_size = 34).next_to(pr3, DOWN)
        pr5 = Tex(r"$ 0  = \partial_u \vert s_u \vert ^2  = 2 s_u \cdot s_{uu} $ ", r" $\Rightarrow$ ", r"$s_{uu} = 0$", font_size = 34).next_to(pr4, DOWN)
        pr6 = Tex(r"$ 0  = - N _{uv} = -(N_v)_u$", r" $ = (k s_v ) _u$", r" $ =  k_u s_v + k s_{uv} $ " , r" $\Rightarrow$ ", r"$s_{uv} = - \frac{k_u}{k} s_v$", font_size = 34).next_to(pr5, DOWN)
        pr7 = Tex(r"$ 0  = - N _{uuv} = -(N_v)_{uu}$", r" $ = (k s_v ) _{uu} $", r" $ =  k_{uu} s_v + 2 k_u s_{uv}  + k $ ", r"$ s_{uuv}$", r" $ = $ ", r"$[ k_{uu} - \frac{2(k_u)^2}{k} ]$", r"$ s_v$ ", font_size = 34).next_to(pr6, DOWN)
        l = le.get_left()[0]        
        le1.shift((le1.get_left()[0] - l)*LEFT)
        le2.shift((le2.get_center()[0])*LEFT)
        pr.shift((pr.get_left()[0] - l)*LEFT)
        pr1.shift((pr1.get_left()[0] - l)*LEFT)
        pr2.shift((pr2.get_left()[0] - l)*LEFT)
        pr3.shift((pr3.get_left()[0] - l)*LEFT)
        pr4.shift((pr4.get_center()[0])*LEFT)
        pr5.shift((pr5.get_left()[0] - l)*LEFT)
        pr6.shift((pr6.get_left()[0] - l)*LEFT)
        pr7.shift((pr7.get_left()[0] - l)*LEFT)
        sdli = Line([0,0,0], [1.8,0,0], color = PINK ).next_to(pr1[2], DOWN).shift(0.15*UP)
        suusvli = Line([0,0,0], [1.5,0,0], color = YELLOW ).next_to(pr2[6], DOWN).shift(0.15*UP)
        kkli = Line([0,0,0], [1.8,0,0], color = GREEN ).next_to(pr6[4], DOWN).shift(0.15*UP)
        xuli = Line([0,0,0], [0.4,0.2,0], color = RED).next_to(pr2[3], DOWN).shift(0.5*UP)
        suli = Line([0,0,0], [0.8,0.2,0], color = RED).next_to(pr3[2], DOWN).shift(0.5*UP)
        suuli = Line([0,0,0], [0.8,0.2,0], color = YELLOW).next_to(pr3[4], DOWN).shift(0.5*UP + 0.1*LEFT)
        suuzli = Line([0,0,0], [0.8,0,0], color = ORANGE).next_to(pr5[2], DOWN).shift(0.15*UP)
        svuuli = Line([0,0,0], [0.4,0.2,0], color = ORANGE).next_to(pr7[3], DOWN).shift(0.5*UP + 0.1*LEFT)
        kkreli = Line([0,0,0], [0.3,0,0], color = GREEN ).next_to(pr7[4], DOWN).shift(0.15*UP)
        kkkli = Line([0,0,0], [1.5,0,0], color = PINK).next_to(pr7[5], DOWN).shift(0.15*UP)
        self.play(Create(le), Create(le1), Create(le2))
        self.wait()
        self.play(Create(pr), Create(pr1[0]))
        self.wait()
        self.play(Create(pr1[1]), Create(pr1[2]))
        self.wait()
        self.play(Create(sdli))
        self.wait()
        self.play(Create(pr2[0]), Create(pr2[1]))
        self.wait()
        self.play(Create(pr2[2]), Create(pr2[3]))
        self.wait()
        self.play(Create(xuli))
        self.wait()
        self.play(Create(pr2[4]))
        self.wait()
        self.play(Create(pr2[5]), Create(pr2[6]), Create(suusvli))
        self.wait()
        self.play(Create(pr3[0]))
        self.wait()
        self.play(Create(pr3[1]), Create(pr3[2]), Create(pr3[3]), Create(pr3[4]), Create(pr3[5]))
        self.wait()
        self.play(Create(suli), Create(suuli),Create(pr3[6]))
        self.wait()
        self.play(Create(pr4))
        self.wait()
        self.play(Create(pr5[0]), Create(pr5[1]))
        self.wait()
        self.play(Create(pr5[2]), Create(suuzli))
        self.wait()
        self.play(Create(pr6[0]))
        self.wait()
        self.play(Create(pr6[1]))
        self.wait()
        self.play(Create(pr6[2]))
        self.wait()
        self.play(Create(pr6[3]), Create(pr6[4]), Create(kkli))
        self.wait()
        self.play(Create(pr7[0]))
        self.wait()
        self.play(Create(pr7[1]))
        self.wait()
        self.play(Create(pr7[2]), Create(pr7[3]))
        self.wait()
        self.play(Create(pr7[4]), Create(pr7[5]), Create(pr7[6]), Create(kkkli),  Create(svuuli), Create(kkreli))
        self.wait()
        



        
class long(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/16)
        le = Tex(r"Lemma", font_size = 34, color = BLUE).to_edge(UL).shift(0.7*RIGHT + 0.6*DOWN)
        le1 = Tex(r"Let ", r"$\gamma$", r" $: \mathbb{R} \to \mathbb{R}^3$  be the line with $\gamma (t) = s(t,0)$ for $t$ small.", font_size = 34).next_to(le, DOWN)
        le2 = Tex(r"If ", r"$\Sigma$", r" is proper, then $\gamma (t) \in A$ for all $t \in \mathbb{R}$.", font_size = 34).next_to(le1, DOWN)
        l = le.get_left()[0]
        le1.shift([le1.get_left()[0]- l ] *LEFT)
        le2.shift([le2.get_left()[0]- l ] *LEFT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(le2[1], DOWN).shift(0.15*UP)
        lli = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(le1[1], DOWN).shift(0.15*UP)
        op = 0.2
        h = - 1
        hei = 1.5
        w = 1/2
        line = ParametricFunction(lambda t : [ - t , 0.5 * t - w * PI / 3   , hei * sqrt(3) /2  + h ] , t_range = [-3, 3], color = BLUE)
        s = Surface(
            lambda u, v:  [ - v , 0.5 * v + w*u  ,  - hei * np.sin (u) + h ],
            u_range = [ -PI , PI ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            #resolution = [ 32, 16]
        )
        self.add_fixed_in_frame_mobjects(le, le1, le2, sli, lli)
        self.remove(le, le1, le2, sli, lli)
        self.play(Create(le), Create(le1), Create(le2), Create(sli), Create(lli), Create(s), Create(line))
        self.wait()
        self.play(FadeOut(s), FadeOut(line))
        self.wait()





        
class longp(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI / 6 )
        le = Tex(r"Lemma", font_size = 34, color = BLUE).to_edge(UL).shift(0.7*RIGHT + 0.6*DOWN)
        le1 = Tex(r"Let ", r"$\gamma$", r" $: \mathbb{R} \to \mathbb{R}^3$  be the line with $\gamma (t) = s(t,0)$ for $t$ small.", font_size = 34).next_to(le, DOWN)
        le2 = Tex(r"If ", r"$\Sigma$", r" is proper, then $\gamma (t) \in A$ for all $t \in \mathbb{R}$.", font_size = 34).next_to(le1, DOWN)
        l = le.get_left()[0]
        le1.shift([le1.get_left()[0]- l ] *LEFT)
        le2.shift([le2.get_left()[0]- l ] *LEFT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(le2[1], DOWN).shift(0.15*UP)
        lli = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(le1[1], DOWN).shift(0.15*UP)
        state = VGroup(le, le1, le2, sli, lli)
        pr = Tex(r"Proof", font_size = 34, color = BLUE).next_to(le2, DOWN)
        pr1 = Tex(r"Assume the result is false; there is $t > 0$ with $\gamma (t) $ not in $A$.", font_size = 34).next_to(pr, DOWN)
        pr2 = Tex(r"Set $T = \sup \{ t  > 0  \, \vert \, \gamma (u) \in A  $ for all $u \in [0,t] \}$.", font_size = 34).next_to(pr1, DOWN)
        pr3 = Tex(r"Case 1", r" : ", r"$\gamma (T)$", r" $ \in A$.", font_size = 34).next_to(le, RIGHT)
        pr4 = Tex(r"From the first lemma, there is a segment ", r"$L$", r" $\subset A$ with $\gamma (T)$ in its interior.", font_size = 34).next_to(pr3, DOWN)
        pr5 = Tex(r"By uniqueness, $L$ is an extension of $\gamma $, contradicting the choice of $T$.", font_size = 34).next_to(pr4, DOWN)
        pr6 = Tex(r"Case 2", r" : ", r"$\gamma (T)$", r" is not in $ A$.", font_size = 34).next_to(pr5, DOWN)
        pr7 = Tex(r"Since $\gamma (t) \in A \subset \Sigma$ for all $t < T$, then $\gamma (T) \in \Sigma $ by properness.", font_size = 34).next_to(pr6, DOWN)
        pr8 = Tex(r"$\Rightarrow $ the shape operator at $\gamma (T)$ is identically zero.", font_size = 34).next_to(pr7, DOWN)
        pr.shift([pr.get_left()[0]- l ] *LEFT)
        pr1.shift([pr1.get_left()[0]- l ] *LEFT)
        pr2.shift([pr2.get_left()[0]- l ] *LEFT)
        pr3.shift([pr3.get_center()[0]] *LEFT)
        pr4.shift([pr4.get_left()[0]- l ] *LEFT)
        pr5.shift([pr5.get_left()[0]- l ] *LEFT)
        pr6.shift([pr6.get_center()[0]] *LEFT)
        pr7.shift([pr7.get_left()[0]- l ] *LEFT)
        pr8.shift([pr8.get_center()[0]] *LEFT)
        c1li = Line([0,0,0], [1,0,0]).next_to(pr3[0], DOWN).shift(0.15*UP)
        c2li = Line([0,0,0], [1,0,0]).next_to(pr6[0], DOWN).shift(0.15*UP)
        gtli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(pr3[2], DOWN).shift(0.15*UP)
        gtli2 = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(pr6[2], DOWN).shift(0.15*UP)
        lli2 = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(pr4[1], DOWN).shift(0.15*UP)
        begp = VGroup(pr, pr1, pr2)
        h = - 0.5
        s = Surface(
            lambda u, v:  [ u * np.cos(v)  ,  u * np.sin(v) ,  - u +h   ],
            u_range = [ 2/3 , 2  ],
            v_range = [ -PI /2 , PI / 2 + PI/3  ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = 0.2, 
            resolution = [ 16, 16]
        )
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift( [ 2/3 + 4/9 , 0 ,  -2/3 - 4/9 + h ] )
        g = ParametricFunction(lambda t : [  - t , 0 , t + h ] , t_range = [ - 2 + 2/9  ,  -2/3 -4/9  ], color = BLUE)
        g2 = ParametricFunction(lambda t : [  - t , 0 , t + h ] , t_range = [  -2/3 -6/9 , -2/3 -2/9 ], color = YELLOW)
        g3 = ParametricFunction(lambda t : [  - t , 0 , t + h ] , t_range = [  -2/3 -6/9 , -2/3 -2/9 ], color = BLUE)
        self.add_fixed_in_frame_mobjects(le, le1, le2, pr, pr1, pr2, pr3, pr4, pr5, sli, lli, c1li, gtli, lli2, pr6, pr7, gtli2, c2li, pr8)
        self.remove(pr, pr1, pr2, pr3, pr4, pr5, gtli, lli2, c1li, pr6, pr7, gtli2, c2li, pr8)
        self.wait()
        self.play(Write(pr), Write(pr1))
        self.wait()
        self.play(Write(pr2))
        self.wait()
        self.play(Create(s), Create(p), Create(g))
        self.wait()
        self.play(Write(pr3), Create(c1li), Create(gtli), FadeOut(state))
        self.wait()
        self.play(Write(pr4), Create(lli2))
        self.wait()
        self.play(Create(g2))
        self.wait()
        self.play(Write(pr5))
        self.wait()
        self.play(g2.animate.become(g3))
        self.wait()
        self.play(Write(pr6), Create(c2li), Create(gtli2), FadeOut(begp), FadeOut(g3), FadeOut(g2))
        self.wait()
        self.play(Write(pr7))
        self.wait()
        self.play(Write(pr8))
        self.wait()



class graph(ThreeDScene):
    def construct(self):
        ax = Axes(
            x_range=[0, 10], y_range=[0, 10]
        )
        labels = ax.get_axis_labels(x_label="t", y_label=" 1 / k ")
        gra = Line([-5,-2.45,0], [2,0,0], color = PINK)
        g = ParametricFunction(lambda t : [  t + 2 + sqrt(20/7) ,  - 1 / t - sqrt(7/20) , 0  ] , t_range = [ - sqrt(20 / 7) , -0.1  ], color = PINK)
        self.play(Create(ax), Create(labels))
        self.wait()
        self.play(Create(gra))
        self.wait()
        self.play(Create(g))
        self.wait()



class acor(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/16)
        t = Tex(r"Corollary", font_size = 34, color = BLUE).to_edge(UL).shift(0.8*DOWN+0.3*RIGHT)
        t1 = Tex(r"If ", r"$\Sigma$", r"  $\subset \mathbb{R}^3$ is a proper surface with zero Gauss curvature, then for all ", r"$p$", r" $\in A$, ", font_size = 34).next_to(t, DOWN)
        t2 = Tex(r"there is a unique line ", r"$\ell (p)$", r" passing through $p$ entirely contained in $A$.", font_size = 34).next_to(t1, DOWN)
        l = t.get_left()[0]
        t1.shift([t1.get_left()[0]- l ] *LEFT)
        t2.shift([t2.get_left()[0]- l ] *LEFT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(t1[1], DOWN).shift(0.15*UP)
        lli = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(t2[1], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(t1[3], DOWN).shift(0.15*UP)
        op = 0.2
        h = - 1.3
        hei = 1.5
        w = 1/2
        p = Sphere(radius = 0.07)
        p.set_color(ORANGE).shift( [ 0, - w * PI / 3  , hei * sqrt(3) /2 + h ] )
        lp = ParametricFunction(lambda t : [ - t , 0.5 * t - w * PI / 3  , hei * sqrt(3) /2  + h ] , t_range = [-3, 3], color = BLUE)
        s = Surface(
            lambda u, v:  [ - v , 0.5 * v + w*u  ,  - hei * np.sin (u) + h ],
            u_range = [ -PI , PI ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op
        )
        self.add_fixed_in_frame_mobjects(t, t1, t2, sli, lli, pli)
        self.remove(t, t1, t2, sli, lli, pli)
        self.play(Create(t), Create(t1), Create(t2), Create(sli), Create(lli), Create(pli), Create(s))
        self.wait()
        self.play(Create(p))
        self.wait()
        self.play(Create(lp))
        self.wait()




class comp(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/6)
        p = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift(0.3*DOWN+0.4*RIGHT)
        p1 = Tex(r"Let ", r"$\Sigma$", r"  $\subset \mathbb{R}^3$ be a complete surface, ", r"$p$", r" $\in \Sigma$, ", font_size = 34).next_to(p, DOWN)
        p2 = Tex(r"$K : \Sigma \to \mathbb{R}$ its Gauss curvature, and $r_0 > 0 $ such that $\text{exp}_p : $ ", r"$T_p \Sigma$", r" $ \to \Sigma$ is", font_size = 34).next_to(p1, DOWN)
        p3 = Tex(r"non-singular in $B(0, r_0) \subset T_p \Sigma $. Then", font_size = 34).next_to(p2, DOWN)
        p4 = Tex(r"$\bullet $ If $K \leq 0$ in the ball of radius $r_0$ around $p$, then", font_size = 34).next_to(p3, DOWN)
        p5 = Tex(r"$ \vert d _x (\text{exp}_p) (V) \vert \geq \vert V \vert $ for all ", r"$x$", r" $ \in B(0, r_0)$, ", r"$V$", r" $ \in T_x(T_p \Sigma )$.", font_size = 34).next_to(p4, DOWN)
        l = p.get_left()[0]
        p1.shift([p1.get_left()[0]- l ] *LEFT)
        p2.shift([p2.get_left()[0]- l ] *LEFT)
        p3.shift([p3.get_left()[0]- l ] *LEFT)
        p4.shift([p4.get_left()[0]- l ] *LEFT)
        p5.shift([p5.get_left()[0]- l - 0.3 ] *LEFT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(p1[1], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(p1[3], DOWN).shift(0.15*UP)
        tpsli = Line([0,0,0], [0.6,0,0], color = GREEN).next_to(p2[1], DOWN).shift(0.15*UP)
        xli = Line([0,0,0], [0.4,0,0], color = RED).next_to(p5[1], DOWN).shift(0.15*UP)
        vli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(p5[3], DOWN).shift(0.15*UP)
        h = - 2
        op = 0.2
        y = 1.7
        xdi = 1
        vle = 1.7
        s = Surface(
            lambda u, v:  [ u - y , v + y * sqrt(3) , (v**2 - u**2)/12 + h ],
            u_range = [ - 2 , 2 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        tps = Surface(
            lambda u, v:  [ u + y  , v - y * sqrt(3) , h ],
            u_range = [ - 2 , 2 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        sco = s.copy()
        tpsco = tps.copy()
        po = Sphere(radius = 0.05)
        z = po.copy()
        po.set_color(ORANGE).shift( [ - y , y * sqrt(3) , h ] )
        z.set_color(ORANGE).shift( [ y , - y * sqrt(3) , h ] )
        x = Sphere(radius = 0.05)
        x.set_color(RED).shift( [ y , - y * sqrt(3) + xdi , h ] )
        ex = Sphere(radius = 0.05)
        ex.set_color(RED).shift( [ - y ,  y * sqrt(3) + xdi ,  xdi**2 / 12 + h ] )
        v = Arrow( 
            [ y - (0.1) * vle , - y * sqrt(3) + xdi , h ], 
            [ y + vle , - y * sqrt(3) + xdi , h ],
            color = YELLOW, 
            max_tip_length_to_length_ratio = 0.2 
        ) 
        dv = Arrow( 
            [ - y - (0.1) * 1.1 * vle ,  y * sqrt(3) + xdi , xdi**2 / 12 + h ], 
            [ - y + 1.1 * vle ,  y * sqrt(3) + xdi , xdi**2 / 12 + h ],
            color = YELLOW, 
            max_tip_length_to_length_ratio = 0.2 
        ) 
        ar = Arrow([-1,-1.4,0], [1,-1.4,0], max_tip_length_to_length_ratio = 0.1 , stroke_width = 3)
        arla = Tex(r"$\text{exp}_p$", font_size = 34).next_to(ar, UP).shift(0.15*DOWN)
        allfig = VGroup(s, tps, po, z, x, ex, v, dv, ar, arla)
        self.add_fixed_in_frame_mobjects(p, p1, p2, p3, p4, p5, sli, pli, xli, vli, tpsli, ar, arla)
        self.remove(p, p1, p2, p3, p4, p5, sli, pli, xli, vli, tpsli, ar, arla)
        self.play(Create(p), Create(p1), Create(p2), Create(p3), Create(p4), Create(p5), Create(sli), Create(pli), Create(xli), Create(vli), Create(tpsli))
        self.wait()
        self.play(Create(s), Create(po), Create(tps), Create(ar), Create(arla), Create(z))
        self.wait()
        self.add(tpsco)
        self.play(Transform(tpsco, sco))
        self.remove(tpsco, sco)
        self.wait()
        self.play(Create(x), Create(ex), Create(v), Create(dv))
        self.wait()
        self.play(FadeOut(allfig))
        self.wait()




class comp2(Scene):
    def construct(self):
        p = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift(0.3*DOWN+0.4*RIGHT)
        p1 = Tex(r"Let ", r"$\Sigma$", r"  $\subset \mathbb{R}^3$ be a complete surface, ", r"$p$", r" $\in \Sigma$, ", font_size = 34).next_to(p, DOWN)
        p2 = Tex(r"$K : \Sigma \to \mathbb{R}$ its Gauss curvature, and $r_0 > 0 $ such that $\text{exp}_p : $ ", r"$T_p \Sigma$", r" $ \to \Sigma$ is", font_size = 34).next_to(p1, DOWN)
        p3 = Tex(r"non-singular in $B(0, r_0) \subset T_p \Sigma $. Then", font_size = 34).next_to(p2, DOWN)
        p4 = Tex(r"$\bullet $ If $K \leq 0$ in the ball of radius $r_0$ around $p$, then", font_size = 34).next_to(p3, DOWN)
        p5 = Tex(r"$ \vert d _x (\text{exp}_p) (V) \vert \geq \vert V \vert $ for all ", r"$x$", r" $ \in B(0, r_0)$, ", r"$V$", r" $ \in T_x(T_p \Sigma )$.", font_size = 34).next_to(p4, DOWN)
        p6 = Tex(r"$\bullet $ If $K \geq 0$ in the ball of radius $r_0$ around $p$, then", font_size = 34).next_to(p5, DOWN)
        p7 = Tex(r"$ \vert d _x (\text{exp}_p) (V) \vert \leq \vert V \vert $ for all ", r"$x$", r" $ \in B(0, r_0)$, ", r"$V$", r" $ \in T_x(T_p \Sigma )$.", font_size = 34).next_to(p6, DOWN)
        l = p.get_left()[0]
        p1.shift([p1.get_left()[0]- l ] *LEFT)
        p2.shift([p2.get_left()[0]- l ] *LEFT)
        p3.shift([p3.get_left()[0]- l ] *LEFT)
        p4.shift([p4.get_left()[0]- l ] *LEFT)
        p5.shift([p5.get_left()[0]- l - 0.3 ] *LEFT)
        p6.shift([p6.get_left()[0]- l ] *LEFT)
        p7.shift([p7.get_left()[0]- l - 0.3 ] *LEFT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(p1[1], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.4,0,0], color = ORANGE).next_to(p1[3], DOWN).shift(0.15*UP)
        tpsli = Line([0,0,0], [0.6,0,0], color = GREEN).next_to(p2[1], DOWN).shift(0.15*UP)
        xli = Line([0,0,0], [0.4,0,0], color = RED).next_to(p5[1], DOWN).shift(0.15*UP)
        vli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(p5[3], DOWN).shift(0.15*UP)
        xli2 = Line([0,0,0], [0.4,0,0], color = RED).next_to(p7[1], DOWN).shift(0.15*UP)
        vli2 = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(p7[3], DOWN).shift(0.15*UP)
        self.add(p, p1, p2, p3, p4, p5, sli, pli, xli, vli, tpsli)
        self.wait()
        self.play(Create(p6), Create(p7), Create(xli2), Create(vli2))
        self.wait()



class polar(Scene):
    def construct(self):
        h = 3
        stroke = 1
        num = 4
        dis = 2
        hl = []
        vl = []
        sr = []
        sr2 = []
        st = []
        st2 = []
        vfr = VGroup()
        vft = VGroup()
        vft2 = VGroup()
        grid = VGroup()
        for i in range(0,2*num + 1):
            hl.append(Line([ - h , - h + i * h / num , 0 ], [ h , -h + i * h / num , 0 ], stroke_width = stroke ))
            vl.append(Line([ - h + i * h / num , - h , 0 ], [ -h + i * h / num , h , 0 ], stroke_width = stroke ))
            grid.add(hl[i], vl[i])
        def scl(i):
            return 1 - 1 / ( 3 * sqrt( 4 + ( - 2 + i) ** 2 ))
        def scb(i):
            return 1 - 1 / ( 3 * sqrt( 1 + ( - 2 + i) ** 2 ))
        def lenl(i):
            return 1 + sqrt(2)/sqrt( 4 + ( - 2 + i ) ** 2)
        def lenb(i):
            return 1 + sqrt(2)/sqrt( 1 + ( - 2 + i ) ** 2)
        for i in range(0,5):
            sr.append( Arrow( 
                [ ( - h / 2 + i * h / num ) * scl(i)  , h * scl (i) / 2 , 0 ] ,   
                [ ( - h / 2 + i * h / num ) * lenl(i) , h * lenl(i) / 2 , 0 ] , 
                color = BLUE ) 
                ) # top ..... length = h * sqrt( 4 + ( 2 + j )**2 )
            sr.append( Arrow( 
                [ ( - h / 2 + i * h / num ) * scl (i)  , ( - h / 2 ) * scl (i)  , 0 ] , 
                [ ( - h / 2 + i * h / num ) * lenl (i) , - h * lenl(i) / 2 , 0 ] , 
                color = BLUE ) 
                ) # bottom
            sr.append( Arrow( 
                [ ( - h / 2 + i * h / num ) * scb(i) , ( h / 4 ) * scb(i) , 0 ] ,   
                [ ( - h / 2 + i * h / num ) * lenb(i), h * lenb(i) / 4 , 0 ] , 
                color = BLUE ) 
                ) # up ..... length = h * sqrt( 1 + ( 2 + j )**2 )
            sr.append( Arrow( 
                [ ( - h / 2 + i * h / num ) * scb(i) , ( - h / 4 ) * scb(i) , 0 ] , 
                [ ( - h / 2 + i * h / num ) * lenb(i), - h * lenb(i) / 4 , 0 ] , 
                color = BLUE ) 
                ) # down
        for i in range (0,20):#number is weird here
                vfr.add(sr[ i ])
        sr2.append( Arrow( [ - h / 2 + h /12 , 0 , 0 ] , [ - h / 2  - sqrt(2) * h / 4 , 0 , 0] , color = BLUE ) )
        sr2.append( Arrow( [   h / 2 - h /12 , 0 , 0 ] , [   h / 2  + sqrt(2) * h / 4 , 0 , 0] , color = BLUE ) )
        sr2.append( Arrow( [ - h / 4 + h /12 , 0 , 0 ] , [ - h / 4  - sqrt(2) * h / 4 , 0 , 0] , color = BLUE ) )
        sr2.append( Arrow( [   h / 4 - h /12 , 0 , 0 ] , [   h / 4  + sqrt(2) * h / 4 , 0 , 0] , color = BLUE ) )
        vfr.add(sr2[0], sr2[1], sr2[2], sr2[3])
        for i in range(0,4):
            for j in range(0,3):
                st.append( Arrow(
                    [ ( j + 1 ) * ( np.cos( PI * i / 2 ) - np.sin( PI * i / 2 ) / (3 * (j + 1) ) ) * h / 4 , ( j + 1 ) * ( np.sin( PI * i / 2 ) + np.cos( PI * i / 2 ) / (3 * (j + 1) ) ) * h / 4 , 0 ], 
                    [ ( j + 1 ) * ( np.cos( PI * i / 2 ) + np.sin( PI * i / 2 )                  ) * h / 4 , ( j + 1 ) * ( np.sin( PI * i / 2 ) - np.cos( PI * i / 2 )                  ) * h / 4 , 0 ],
                    color = YELLOW )
                    )
                st2.append( Arrow(
                    [ ( j + 1 ) * ( np.cos( PI * i / 2 ) - np.sin( PI * i / 2 ) / (3 * (j + 1) )    ) * h / 4 , ( j + 1 ) * ( np.sin( PI * i / 2 ) + np.cos( PI * i / 2 ) / (3 * (j + 1) )     ) * h / 4 , 0 ], 
                    [ ( j + 1 ) * ( np.cos( PI * i / 2 ) + np.sin( PI * i / 2 ) * sqrt(2) / (j + 1) ) * h / 4 , ( j + 1 ) * ( np.sin( PI * i / 2 ) - np.cos( PI * i / 2 ) * sqrt(2) / (j + 1 ) ) * h / 4 , 0 ],
                    color = GREEN )
                    )
                vft.add( st [ 5 * i + j ])
                vft2.add(st2[ 5 * i + j ])
            for j in range(0,2):
                st.append( Arrow( 
                    [ sqrt(2) * ( j + 1 ) * ( np.cos( PI * ( i / 2 + 1 / 4 ) ) - np.sin( PI * ( i / 2 + 1 / 4 ) ) / ( sqrt(2) * 3 * (j + 1) ) ) * h / 4 , sqrt(2) * ( j + 1 ) * ( np.sin( PI * ( i / 2 + 1 / 4 ) ) + np.cos( PI * ( i / 2 + 1 / 4 ) ) / ( sqrt(2) * 3 * (j + 1) ) ) * h / 4 , 0 ], 
                    [ sqrt(2) * ( j + 1 ) * ( np.cos( PI * ( i / 2 + 1 / 4 ) ) + np.sin( PI * ( i / 2 + 1 / 4 ) )                             ) * h / 4 , sqrt(2) * ( j + 1 ) * ( np.sin( PI * ( i / 2 + 1 / 4 ) ) - np.cos( PI * ( i / 2 + 1 / 4 ) )                             ) * h / 4 , 0 ] , 
                    color = YELLOW )
                    )
                st2.append( Arrow( 
                    [ sqrt(2) * ( j + 1 ) * ( np.cos( PI * ( i / 2 + 1 / 4 ) ) - np.sin( PI * ( i / 2 + 1 / 4 ) ) / ( sqrt(2) * 3 * (j + 1) ) ) * h / 4 , sqrt(2) * ( j + 1 ) * ( np.sin( PI * ( i / 2 + 1 / 4 ) ) + np.cos( PI * ( i / 2 + 1 / 4 ) ) / ( sqrt(2) * 3 * (j + 1) ) ) * h / 4 , 0 ], 
                    [ sqrt(2) * ( j + 1 ) * ( np.cos( PI * ( i / 2 + 1 / 4 ) ) + np.sin( PI * ( i / 2 + 1 / 4 ) ) / ( j + 1 )                 ) * h / 4 , sqrt(2) * ( j + 1 ) * ( np.sin( PI * ( i / 2 + 1 / 4 ) ) - np.cos( PI * ( i / 2 + 1 / 4 ) ) / ( j + 1 )                 ) * h / 4 , 0 ] , 
                    color = GREEN )
                    )
                vft.add( st [ 5 * i + 3 + j ])
                vft2.add(st2[ 5 * i + 3 + j ])
        grid.shift([ - dis , 0 , 0 ])
        vfr.shift([ - dis , 0, 0 ])
        vft.shift([ - dis , 0, 0 ])
        vft2.shift([ - dis , 0, 0 ])
        partr = Tex(r"$\partial _r $", font_size = 34).shift(3.2*RIGHT + 1.45 * UP)
        partt = Tex(r"$\partial _{\theta} $", font_size = 34).next_to(partr, RIGHT)
        lpr = Tex(r"$\vert \partial _r \vert = 1 $", font_size = 34).next_to(partt, DOWN).shift(0.25*LEFT + 0.3 * DOWN)
        lpt = Tex(r"$\vert \partial _{\theta} \vert = r $", font_size = 34).next_to(lpr, DOWN).shift(0.3*DOWN)
        npt = Tex(r"$\vert \frac{1}{r} \partial _{\theta} \vert = 1 $", font_size = 34).next_to(lpt, DOWN).shift(0.3*DOWN)
        prli = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(partr, DOWN).shift(0.15*UP)
        ptli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(partt, DOWN).shift(0.15*UP)
        nli = Line([0,0,0], [1.3,0,0], color = GREEN).next_to(npt, DOWN).shift(0.15*UP)
        self.play(Create(grid))
        self.wait()
        self.play(Create(vfr), Create(vft), Create(partr), Create(partt), Create(prli), Create(ptli))
        self.wait()
        self.play(Create(lpr), Create(lpt) )
        self.wait()
        self.play(FadeOut(vfr), FadeOut(vft), Create(vft2), Create(npt), Create(nli))
        self.wait()
        


class compp(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES , theta =  PI / 6 )
        p = Tex(r"Proof", font_size = 34, color = BLUE).to_edge(UL).shift(0.3*DOWN+0.5*RIGHT)
        p1 = Tex(r"Let ", r"$U$", r" $ = (0, r_0) \times (- \pi , \pi ) $, and $s : U \to \Sigma$ given by", font_size = 34).next_to(p, DOWN)
        p2 = Tex(r"$ s ( r , \theta ) : = \text{exp}_p ( r \cos \theta , r \sin \theta ).$", font_size = 34).next_to(p1, DOWN)
        p3 = Tex(r"$ b ( r , \theta ) : = \vert s_{\theta} \vert $, ", r" $b_{rr} + K b = 0 $, ", r" $ b(r) = b(r, 0) $.", font_size = 34).next_to(p2, DOWN)
        p4 = Tex(r"$ \dfrac{b (r)}{r} = \dfrac{ \vert  \text{d exp}_p  (  \partial _{\theta}) \vert }{r}  $ ", r"$ = \vert \text{d exp}_p \left(  \frac{\partial_{\theta} }{r}  \right) \vert $", font_size = 34).next_to(p3, DOWN)
        p5 = Tex(r"$\lim$ ", r" $\dfrac{b (r)}{r} = 1$", font_size = 34).next_to(p4, DOWN)
        l = p.get_left()[0]
        p1.shift([p1.get_left()[0]- l ] *LEFT)
        p2.shift([p2.get_center()[0]]*LEFT)
        p3.shift([p3.get_center()[0]]*LEFT)
        p4.shift([p4.get_center()[0]]*LEFT)
        p5.shift([p5.get_center()[0]]*LEFT)
        p5[0].shift(0.11*UP + 0.1*LEFT)
        p6 = Tex(r"$r \to 0 $", font_size = 18).next_to(p5[0], DOWN).shift(0.15*UP)
        p7 = VGroup(p5, p6)
        uli = Line([0,0,0], [0.4,0,0], color = GREEN).next_to(p1[1], DOWN).shift(0.15*UP)
        tli = Line([0,0,0], [0.3,0,0], color = YELLOW).shift( 0.45 * UP + 0.36 * LEFT )
        trli = Line([0,0,0], [0.3,0,0], color = GREEN).shift( 2.13* RIGHT + 0.05 * UP )
        box = SurroundingRectangle(p7, color = PINK , buff = 0.1)
        h = - 2.2
        op = 0.2
        sh = 2.6
        w = 0.83
        a = 1
        s = Surface(
            lambda u, v:  [ w * u  - sh , w * v + sh * sqrt(3) , w* (v**2 - u**2) / 12 + h ],
            u_range = [ - 2 , 2 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        tps = Surface(
            lambda u, v:  [ w * u  , w * v , h ],
            u_range = [ - 2 , 2 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        rtp = Surface(
            lambda u, v:  [ w * u  +  sh  , w * v - sh * sqrt(3) , h ],
            u_range = [  -2 , 2 ],
            v_range = [  0 , 2 ],
            checkerboard_colors = [GREEN, GREEN ],
            fill_opacity = 2 * op, 
            resolution = [ 16, 8]
        )
        rlines = []
        rcurves = []
        thlines = []
        thcurves = []
        rgroup = VGroup()
        thgroup = VGroup()
        for i in range(0, 9):
            rlines.append( Line ( [ sh + w * ( - 2 +  i / 2 )  , - sh * sqrt(3) , h ] , [ sh + w * ( - 2 + i / 2) , - sh * sqrt(3) + 2 * w , h ] , color = BLUE))
            rgroup.add(rlines[i])
        for i in range(0,8):
            rcurves.append( ParametricFunction(lambda t : [ t * np.cos( PI * i / 4 ) , t * np.sin( PI * i / 4 ) , h ] , t_range = [ 0 , 2 * w ], color = BLUE) )
            rgroup.add(rcurves[i])
        for i in range(0, 3):    
            thlines.append( Line ( [ sh - 2 * w , - sh * sqrt(3) + w * 2 * ( i + 1 ) / 3 , h ] , [ sh + 2 * w , - sh * sqrt(3) + w * 2 * ( i + 1 ) / 3 , h ] , color = YELLOW) )  
            thcurves.append( ParametricFunction( lambda t : [ 2 * w * ( i + 1 ) * np.cos(t) / 3 , 2 * w * ( i + 1 ) * np.sin(t) / 3 , h ], t_range = [0, TAU] , color = YELLOW) )
            thgroup.add(thlines[i], thcurves[i])       
        axu = ThreeDAxes(x_range = [-3,3,1], y_range = [-3,3,1], z_range = [0,3,1], x_length = 4 * w , y_length = 4 * w  , z_length = 0.0001 , tips = False ).shift( [ sh , - sh * sqrt(3) , h ] )
        axt = ThreeDAxes(x_range = [-3,3,1], y_range = [-3,3,1], z_range = [0,3,1], x_length = 4 * w , y_length = 4 * w  , z_length = 0.0001 , tips = False ).shift( [ 0 , 0  , h ] )
        axs = ThreeDAxes(x_range = [-3,3,1], y_range = [-3,3,1], z_range = [-2,2,1], x_length = 4 * w , y_length = 4 * w , z_length = 8*w/3  , tips = False ).shift( [ - sh , sh * sqrt(3) , h ] )
        uco = rtp.copy()
        tco = tps.copy()
        sco = s.copy()
        po = Sphere(radius = 0.05)
        z = po.copy()
        po.set_color(ORANGE).shift( [ -sh , sh *sqrt(3) , h ] )
        z.set_color(ORANGE).shift( [ 0 , 0 , h ] )
        thv = []
        pth = []
        sth = []
        theta = VGroup()
        for i in range(0, 6):
            thv.append(Arrow( 
                [ sh + 0.2 , -sh *sqrt(3) +  a * w * ( i + 1 ) / 4 , h ], 
                [ sh - 1.5 * w * ( i + 1 ) / 4 , -sh *sqrt(3) + a * w * ( i + 1 ) / 4, h ],
                color = YELLOW, 
                max_tip_length_to_length_ratio = 0.1
            ))
            pth.append(Arrow( 
                [ 0.2 , a * w * ( i + 1 ) / 4 , h ], 
                [ - 1.5 * w * ( i + 1 ) / 4 , a * w * ( i + 1 ) / 4 , h ],
                color = YELLOW, 
                max_tip_length_to_length_ratio = 0.1 
            ))
            sth.append(Arrow( 
                [ 0.2 - sh , a * w * ( i + 1 ) / 4 + sh * sqrt(3) , w * ( a * ( i + 1 ) / 4 ) **2 / 12 + h ], 
                [ - 1.5 * w * ( i + 1 ) / 4 - sh , a * w * ( i + 1 ) / 4 + sh * sqrt(3) , w * ( a * ( i + 1 ) / 4 ) ** 2 / 12 + h ],
                color = YELLOW, 
                max_tip_length_to_length_ratio = 0.1 
            ))
            theta.add(thv[i], pth[i], sth[i])
        thvr = []
        pthr = []
        sthr = []
        thetar = VGroup()
        for i in range(0, 6):
            thvr.append(Arrow( 
                [ sh + 0.2 , -sh *sqrt(3) +  a * w * ( i + 1 ) / 4 , h ], 
                [ sh - 1.8 * w , -sh *sqrt(3) + a * w * ( i + 1 ) / 4, h ],
                color = GREEN, 
                max_tip_length_to_length_ratio = 0.2 
            ))
            pthr.append(Arrow( 
                [ 0.2 , a * w * ( i + 1 ) / 4 , h ], 
                [ - 1.8 * w , a * w * ( i + 1 ) / 4 , h ],
                color = GREEN, 
                max_tip_length_to_length_ratio = 0.2 
            ))
            sthr.append(Arrow( 
                [ 0.2 - sh , a * w * ( i + 1 ) / 4 + sh * sqrt(3) , w * ( a * ( i + 1 ) / 4 ) **2 / 12 + h ], 
                [ - 1.8 * w - sh , a * w * ( i + 1 ) / 4 + sh * sqrt(3) , w * ( a * ( i + 1 ) / 4 ) ** 2 / 12 + h ],
                color = GREEN, 
                max_tip_length_to_length_ratio = 0.2 
            ))
            thetar.add(thvr[i], pthr[i], sthr[i])
        self.add_fixed_in_frame_mobjects(p, p1, p2, p3, p4, uli, tli, trli, p5, box, p6 )
        self.remove(p, p1, p2, p3, p4, uli, tli, trli, p5, box, p6 )
        self.play(Create(p), Create(p1), Create(p2), Create(s), Create(po), Create(z), Create(tps), Create(rtp), Create(axu), Create(axt), Create(axs), Create(uli))
        self.wait()
        self.play(Create(rgroup), Create(thgroup))
        self.wait()
        self.play(Create(p3[0]), Create(p3[1]))
        self.wait()
        self.play(FadeOut(rgroup), FadeOut(thgroup))
        self.wait()
        self.play(Create(theta))
        self.wait()
        self.play(Create(p3[2]))
        self.wait()
        self.play(Create(p4[0]), Create(tli))
        self.wait()
        self.play(FadeOut(theta))
        self.wait()
        self.play(Create(p4[1]), Create(thetar), Create(trli))
        self.wait()
        self.play(Create(p5), Create(box), Create(p6))
        self.wait()


class graph2(ThreeDScene):
    def construct(self):
        w = 2
        a = 0.7
        jac1 = Tex(r"If $K \leq 0$, then  ", r" $b _{rr} = - K b \geq 0$ ", r" $ \Rightarrow $ $ b $ is convex", font_size = 34).shift( 1.1 * UP + 1.4 * w * LEFT )
        jac2 = Tex(r"$ \Rightarrow b (r) \geq  r  $ for all $r \leq r_0 $", font_size = 34).next_to(jac1, DOWN)
        jac3 = Tex(r"If not,  $ b (\xi)  < \xi $ for some $\xi \in (0, r_0 )$, and ", font_size = 34).next_to(jac2, DOWN)
        jac4 = Tex(r" $ \lim $ ", r"  $ \frac{b (r)}{r} < 1 $ ", font_size = 34).next_to(jac3, DOWN)
        jac3.shift(( jac3.get_left()[0] - jac1.get_left()[0] )* LEFT)
        jac4[0].shift(0.11 * UP + 0.1*LEFT)
        jac5 = Tex(r"$ r \to 0 $", font_size = 17).next_to(jac4[0], DOWN).shift(0.15*UP)
        lim = VGroup(jac4, jac5)
        lin = Line([0,0,0], [1.5,0,0], color = RED).next_to(lim, DOWN).shift(0.15*UP)
        ax = Axes(
            x_range=[0, 10], y_range=[0, 10] , x_length = 2 * w , y_length = 2 * w
        ).shift( [ 2 * w , - a * w + 1, 0 ])
        labels = ax.get_axis_labels(x_label="r", y_label="b(r)")
        labels[0].set_font_size(34)
        labels[1].set_font_size(34)
        g = ParametricFunction(lambda t : [  w * t  + w , w * t + w * t ** 3 / 16 - ( 1 + a ) * w + 1 , 0  ] , t_range = [ 0 , 1.9 ], color = PINK)
        g2 = ParametricFunction(lambda t : [  w * t  + w , 3 * w * t / 4 + w * t ** 3 / 6 - ( 1 + a ) * w + 1 , 0  ] , t_range = [ 0 , 1.8 ], color = RED)
        r = DashedLine([ w , - ( 1 + a ) * w + 1 , 0 ], [ (2.6) * w , ( 1 - a - 0.4 ) * w + 1 , 0 ] )# dash_length=2.0).shift(UP*2)
        self.play(Create(ax), Create(labels), Create(jac1))
        self.wait()
        self.play(Create(jac2), Create(r), Create(g))
        self.wait()
        self.play(Create(jac3))
        self.wait()
        self.play(Create(jac4), Create(jac5), FadeOut(g), Create(g2), Create(lin))
        self.wait()



class kneg1(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/6)
        kn1 = Tex(r"If $K \leq 0$, then  ", r" $b (r , \theta ) \geq r$  for all $ r \in (0, r_0 ) $, $\theta \in ( - \pi , \pi )$", font_size = 34).to_edge(UL).shift(0.3*RIGHT + DOWN)
        kn2 = Tex(r"$ \Rightarrow $ ", r" $\vert s_{\theta} \vert \geq \vert \partial _{\theta } \vert $ ", font_size = 34).next_to(kn1, DOWN)
        kn3 = Tex(r"If ", r"$V$", r" $ = \alpha \cdot $ ", r"$\partial _r$", r" $ + \beta \cdot $ ", r"$\partial _{\theta }$", r",  then ", font_size = 34).next_to(kn2, DOWN)
        kn4 = Tex(r" $ \vert V \vert ^2 $", r" $ = \vert \alpha \cdot \partial _ r \vert ^2 + \vert \beta \cdot \partial_{\theta} \vert ^2 =  \alpha ^2 + \beta ^2 r ^2 $", font_size = 34).next_to(kn3, DOWN)
        kn5 = Tex(r" $ \leq \alpha ^2 + \beta ^2 \cdot \vert  \text{d exp}_p (\partial _{\theta} ) \vert ^2 $", font_size = 34).next_to(kn4, DOWN)
        kn6 = Tex(r" $ =  \alpha ^2 \cdot \vert  \text{d exp}_p (\partial _r ) \vert ^2  + \beta ^2 \cdot \vert  \text{d exp}_p (\partial _{\theta} ) \vert ^2   $", font_size = 34).next_to(kn5, DOWN)
        kn7 = Tex(r" $ = \vert \text{d exp}_p  ( \alpha \cdot  \partial_r  + \beta \cdot  \partial _{\theta} ) \vert ^2  = \vert \text{d exp}_p (V) \vert ^2$", font_size = 34).next_to(kn6, DOWN)
        l = kn1.get_left()[0]
        kn2.shift(kn2.get_center()[0]*LEFT)
        kn3.shift((l - kn3.get_left()[0])*RIGHT)
        kn4.shift((l - kn4.get_left()[0] + 2.5 )*RIGHT)
        kn5.shift((kn4[1].get_left()[0] - kn5.get_left()[0])*RIGHT)
        kn6.shift((kn4[1].get_left()[0] - kn6.get_left()[0])*RIGHT)
        kn7.shift((kn4[1].get_left()[0] - kn7.get_left()[0])*RIGHT)
        vli = Line([0,0,0], [0.4,0,0], color = PINK ).next_to(kn3[1], DOWN).shift(0.15*UP)
        rli = Line([0,0,0], [0.4,0,0], color = BLUE ).next_to(kn3[3], DOWN).shift(0.15*UP)
        tli = Line([0,0,0], [0.4,0,0], color = YELLOW ).next_to(kn3[5], DOWN).shift(0.15*UP)
        h = - 2
        op = 0.2
        y = 1.7
        xdi = 1
        vle = 1.7
        s = Surface(
            lambda u, v:  [ u - y , v + y * sqrt(3) , (v**2 - u**2)/12 + h ],
            u_range = [ - 2 , 2 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        tps = Surface(
            lambda u, v:  [ u + y  , v - y * sqrt(3) , h ],
            u_range = [ - 2 , 2 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        sco = s.copy()
        tpsco = tps.copy()
        po = Sphere(radius = 0.05)
        z = po.copy()
        po.set_color(ORANGE).shift( [ - y , y * sqrt(3) , h ] )
        z.set_color(ORANGE).shift( [ y , - y * sqrt(3) , h ] )
        x = Sphere(radius = 0.05)
        x.set_color(RED).shift( [ y , - y * sqrt(3) + xdi , h ] )
        ex = Sphere(radius = 0.05)
        ex.set_color(RED).shift( [ - y ,  y * sqrt(3) + xdi ,  xdi**2 / 12 + h ] )
        part = Arrow( 
            [ y - (0.1) * vle , - y * sqrt(3) + xdi , h ], 
            [ y + vle , - y * sqrt(3) + xdi , h ],
            color = YELLOW, 
            max_tip_length_to_length_ratio = 0.2 
        ) 
        dpart = Arrow( 
            [ - y - (0.1) * 1.1 * vle ,  y * sqrt(3) + xdi , xdi**2 / 12 + h ], 
            [ - y + 1.1 * vle ,  y * sqrt(3) + xdi , xdi**2 / 12 + h ],
            color = YELLOW, 
            max_tip_length_to_length_ratio = 0.2 
        )
        parr = Arrow( 
            [ y , - y * sqrt(3) + xdi - (0.1) * vle  , h ], 
            [ y  , - y * sqrt(3) + xdi + 0.7 * vle , h ],
            color = BLUE, 
            max_tip_length_to_length_ratio = 0.2 
        ) 
        dparr = Arrow( 
            [ - y  ,  y * sqrt(3) + xdi - (0.1) * vle , xdi**2 / 12 + h ], 
            [ - y ,  y * sqrt(3) + xdi + 0.7 * vle , xdi**2 / 12 + h + 0.3 ],
            color = BLUE, 
            max_tip_length_to_length_ratio = 0.2 
        )
        v = Arrow( 
            [ y - (0.08) * vle , - y * sqrt(3) + xdi - (0.048) * vle, h ], 
            [ y + vle , - y * sqrt(3) + xdi + 0.6 * vle , h ],
            color = PINK, 
            max_tip_length_to_length_ratio = 0.2 
        ) 
        dv = Arrow( 
            [ - y - (0.08) * 1.1 * vle ,  y * sqrt(3) + xdi - (0.0435) * vle , xdi**2 / 12 + h ], 
            [ - y + 1.1 * vle ,  y * sqrt(3) + xdi + 0.6 * vle , xdi**2 / 12 + h + 0.3 ],
            color = PINK, 
            max_tip_length_to_length_ratio = 0.2 
        ) 
        ar = Arrow([-1,-1.4,0], [1,-1.4,0], max_tip_length_to_length_ratio = 0.1 , stroke_width = 3)
        arla = Tex(r"$\text{exp}_p$", font_size = 34).next_to(ar, UP).shift(0.15*DOWN)
        allfig = VGroup(s, tps, po, z, x, ex, part, dpart, ar, arla, v, dv, parr, dparr)
        self.add_fixed_in_frame_mobjects(kn1, kn2, kn3, kn4, ar, arla, kn5, kn6, kn7, vli, rli, tli)
        self.remove(kn1, kn2, kn3, kn4, ar, arla, kn5, kn6, kn7, vli, tli, rli)
        self.play(Create(kn1))
        self.wait()
        self.play(Create(kn2))
        self.wait()
        self.play(Create(kn3), Create(kn4), Create(vli), Create(rli), Create(tli), Create(allfig))
        self.wait()
        self.play(Create(kn5))
        self.wait()
        self.play(Create(kn6), FadeOut(allfig))
        self.wait()
        self.play(Create(kn7))
        self.wait()
        


class graph3(ThreeDScene):
    def construct(self):
        w = 2
        a = 0.7
        jac1 = Tex(r"If $K \geq 0$, then  ", r" $b _{rr} = - K b \leq 0$ ", r" $ \Rightarrow $ $ b $ is concave", font_size = 34).shift( 1.1 * UP + 1.4 * w * LEFT )
        jac2 = Tex(r"$ \Rightarrow b (r) \leq  r  $ for all $r \leq r_0 $", font_size = 34).next_to(jac1, DOWN)
        jac3 = Tex(r"If not,  $ b (\xi)  > \xi $ for some $\xi \in (0, r_0 )$, and ", font_size = 34).next_to(jac2, DOWN)
        jac4 = Tex(r" $ \lim $ ", r"  $ \frac{b (r)}{r} > 1 $ ", font_size = 34).next_to(jac3, DOWN)
        jac3.shift(( jac3.get_left()[0] - jac1.get_left()[0] )* LEFT)
        jac4[0].shift(0.11 * UP + 0.1*LEFT)
        jac5 = Tex(r"$ r \to 0 $", font_size = 17).next_to(jac4[0], DOWN).shift(0.15*UP)
        lim = VGroup(jac4, jac5)
        lin = Line([0,0,0], [1.5,0,0], color = RED).next_to(lim, DOWN).shift(0.15*UP)
        ax = Axes(
            x_range=[0, 10], y_range=[0, 10] , x_length = 2 * w , y_length = 2 * w
        ).shift( [ 2 * w , - a * w + 1, 0 ])
        labels = ax.get_axis_labels(x_label="r", y_label="b(r)")
        labels[0].set_font_size(34)
        labels[1].set_font_size(34)
        g = ParametricFunction(lambda t : [  w * t  + w , w * t - w * t ** 3 / 16 - ( 1 + a ) * w + 1 , 0  ] , t_range = [ 0 , 1.9 ], color = PINK)
        g2 = ParametricFunction(lambda t : [  w * t  + w , 4 * w * t / 3 - w * t ** 3 / 6 - ( 1 + a ) * w + 1 , 0  ] , t_range = [ 0 , 1.8 ], color = RED)
        r = DashedLine([ w , - ( 1 + a ) * w + 1 , 0 ], [ (2.6) * w , ( 1 - a - 0.4 ) * w + 1 , 0 ] )# dash_length=2.0).shift(UP*2)
        self.play(Create(ax), Create(labels), Create(jac1))
        self.wait()
        self.play(Create(jac2), Create(r), Create(g))
        self.wait()
        self.play(Create(jac3))
        self.wait()
        self.play(Create(jac4), Create(jac5), FadeOut(g), Create(g2), Create(lin))
        self.wait()


class kneg2alt(Scene):
    def construct(self):
        ex = Tex(r"Exercise ", r" (Cartan--Hadamard Theorem)", font_size = 34, color = BLUE).to_edge(UL).shift( 0.8*RIGHT + 1.5 * DOWN)
        ex1 = Tex(r"If $\Sigma \subset \mathbb{R}^3$ is a proper surface of non-positive Gauss curvature,", font_size = 34).next_to(ex, DOWN)
        ex2 = Tex(r"then for all $p \in \Sigma$, the map $\text{exp}_p : T_p \Sigma \to \Sigma $ is non-singular.", font_size = 34).next_to(ex1, DOWN)
        ex3 = Tex(r"Moreover, it is a covering map. ", font_size = 34).next_to(ex2, DOWN)
        hi = Tex(r"Hint", font_size = 34, color = ORANGE).next_to(ex2, DOWN).shift(0.2*DOWN)
        hi1 = Tex(r"Let $r_1 : =  \sup \{ r > 0 \vert \text{exp}_p $ is non-singular in $B(0, r) \}$.", font_size = 34).next_to(hi, DOWN)
        hi2 = Tex(r"Show that $\text{exp}_p$ is non-singular in the ball $B(0, r_1 + \varepsilon )$ for some $\varepsilon > 0 $.", font_size = 34).next_to(hi1, DOWN)
        l = ex.get_left()[0]
        ex1.shift((l - ex1.get_left()[0])*RIGHT)
        ex2.shift((l - ex2.get_left()[0])*RIGHT)
        ex3.shift((l - ex3.get_left()[0])*RIGHT)
        hi.shift((l - hi.get_left()[0])*RIGHT)
        hi1.shift((l - hi1.get_left()[0])*RIGHT)
        hi2.shift((l - hi2.get_left()[0])*RIGHT)
        self.play(Write(ex[0]), Write(ex1), Write(ex2))
        self.wait()
        self.play(Write(hi), Write(hi1))
        self.wait()
        self.play(Write(hi2))
        self.wait()
        self.play(FadeOut(hi), FadeOut(hi1), FadeOut(hi2), Write(ex3))
        self.wait()
        self.play(Write(ex[1]))
        self.wait()
        

def bump0(t):
    if t <= 0:
        return 0
    return 2 ** ( 1 - 1 / t )


def bump1(t):
    return  bump0(t) * bump0(1 - t) 


def bump2(t):
    return  bump1(t+0.5) + bump1(t + 2.5) / 2 + bump1(t - 1.5)/2



class lifts(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=82 * DEGREES, theta=  PI/6)
        le = Tex(r"Let ", r"$\Sigma$", r" $ \subset \mathbb{R}^3$ be a proper surface with zero Gauss curvature and ", r"$p$", r" $ \in \Sigma$.", font_size = 34).to_edge(UL).shift(0.2*RIGHT + 0.4*DOWN)
        le1 = Tex(r"If ",  r"$x$", r" $ \in $ ", r"$T_p \Sigma $", r" is such that $\text{exp}_p (x) \in $ ", r"$A $", r", then there is a line ", r"$\ell ( \exp_p (x) )$", r" $ \subset A$",  font_size = 34).next_to(le, DOWN)
        le2 = Tex(r"containing $\exp (x)$. ", r" Denote the lift by ", r"$\tilde{\ell} ( x )$", r" $ \subset T_p \Sigma$.", font_size = 34).next_to(le1, DOWN)
        le3 = Tex(r"If ",  r"$y$", r" $ \in $ ", r"$T_p \Sigma $", r" is such that $\text{exp}_p (y) \in $ ", r"$A $", r", then there is a line ", r"$\tilde{\ell} ( y )$", r" $ \subset T_p \Sigma$",  font_size = 34).next_to(le2, DOWN)
        l = le.get_left()[0]
        le1.shift([le1.get_left()[0]- l ] *LEFT)
        le2.shift([le2.get_left()[0]- l ] *LEFT)
        le3.shift([le3.get_left()[0]- l ] *LEFT)
        sli1 = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(le[1], DOWN).shift(0.15*UP)
        sli2 = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(le[1], DOWN).shift(0.1*UP)
        pli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(le[3], DOWN).shift(0.15*UP)
        tpsli = Line([0,0,0], [0.5,0,0], color = GREEN).next_to(le1[3], DOWN).shift(0.15*UP)
        xli = Line([0,0,0], [0.3,0,0], color = RED).next_to(le1[1], DOWN).shift(0.15*UP)
        ali = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(le1[5], DOWN).shift(0.15*UP)
        lxli = Line([0,0,0], [1.3,0,0], color = YELLOW).next_to(le1[7], DOWN).shift(0.15*UP)
        tlxli = Line([0,0,0], [0.5,0,0], color = YELLOW).next_to(le2[2], DOWN).shift(0.15*UP)
        yli = Line([0,0,0], [0.3,0,0], color = TEAL).next_to(le3[1], DOWN).shift(0.15*UP)
        tpsli2 = Line([0,0,0], [0.5,0,0], color = GREEN).next_to(le3[3], DOWN).shift(0.15*UP)
        ali2 = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(le3[5], DOWN).shift(0.15*UP)
        tlyli = Line([0,0,0], [0.5,0,0], color = PINK).next_to(le3[7], DOWN).shift(0.15*UP)
        op = 0.3
        h = - 2.5
        hh = -0.7
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift([0,0,h + 3 * bump2(0)])
        z = Sphere(radius = 0.05)
        z.set_color(ORANGE).shift([0,0,hh])
        x = Sphere(radius = 0.05)
        x.set_color(RED).shift([ - 2.3 , 1 ,hh])
        exx = Sphere(radius = 0.05)
        exx.set_color(RED).shift([ - 2.2 * 0.8 , 1 , 3 * bump2( 2.2 ) + h ])
        y = Sphere(radius = 0.05)
        y.set_color(TEAL).shift([ - 2.75 , -0.6 ,hh])
        exy = Sphere(radius = 0.05)
        exy.set_color(TEAL).shift([ - 2.4 * 0.8 , -0.6 , 3 * bump2( 2.4 ) + h ])
        lx = ParametricFunction(lambda t : [  - 2.2 * 0.8 , t , 3 * bump2(2.2)  + h ] , t_range = [-3,3], color = YELLOW)
        tlx = ParametricFunction(lambda t : [  - 2.3 , t , hh ] , t_range = [-3, 3], color = YELLOW)
        yline = ParametricFunction(lambda t : [  - 2.4 * 0.8 , t , 3 * bump2(2.4)  + h ] , t_range = [-3,3], color = PINK)
        yline2 = ParametricFunction(lambda t : [  - 2.75 , t , hh ] , t_range = [-3, 3], color = PINK)
        tps = Surface(lambda u,v: [u,v,hh], u_range = [-3,3], v_range = [-3,3], checkerboard_colors = [GREEN, GREEN], fill_opacity = op, resolution = [16,16])
        s0 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ -3 , -2.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op, 
            resolution = [2, 16]
        )
        s1 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ -2.5 , -1.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s2 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ -1.5 , -0.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s3 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ -0.5 , 0.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s4 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ 0.5 , 1.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s5 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ 1.5 , 2.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s6 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ 2.5 , 3 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op, 
            resolution = [2, 16]
        )
        s = VGroup(s1, s2, s3, s4, s5, p, z, tps, s0, s6)
        self.add_fixed_in_frame_mobjects(le, le1, le2, sli1, sli2, xli, pli, lxli, tlxli, ali, tpsli, le3, yli, tpsli2, tlyli, ali2)
        self.remove(le, le1, le2, sli1, sli2, xli, pli, lxli, tlxli, ali, tpsli, le3, yli, tpsli2, tlyli, ali2)
        self.play(Create(le), Create(s), Create(sli1), Create(sli2), Create(pli))
        self.wait()
        self.begin_ambient_camera_rotation(rate=PI/4)
        self.wait(2)
        self.stop_ambient_camera_rotation()
        self.wait()
        self.play(Create(le1), Create(xli), Create(tpsli), Create(ali), Create(le2[0]), Create(x), Create(exx))
        self.wait()
        self.play(Create(lxli), Create(lx))
        self.wait()
        self.play(Create(le2[1]), Create(le2[2]), Create(le2[3]), Create(tlxli), Create(tlx))
        self.wait()
        self.play(Create(le3), Create(yli), Create(tpsli2), Create(tlyli), Create(ali2), Create(y), Create(exy))
        self.wait()
        self.play(Create(yline), Create(yline2))
        self.wait()






class para1(Scene):
    def construct(self):
        h = 3
        stroke = 1
        num = 4
        dis = 2
        hl = []
        vl = []
        grid = VGroup()
        for i in range(0,2*num + 1):
            hl.append(Line([ - h , - h + i * h / num , 0 ], [ h , -h + i * h / num , 0 ], stroke_width = stroke ))
            vl.append(Line([ - h + i * h / num , - h , 0 ], [ -h + i * h / num , h , 0 ], stroke_width = stroke ))
            grid.add(hl[i], vl[i])
        grid.shift([ - dis , 0 , 0 ])
        x = Dot ( [ -dis - 2 * h / num , - 2 * h / num , 0 ] , color = RED )
        y = Dot ( [ -dis +  h / num , - h / num , 0 ] , color = TEAL )
        z = Dot ( [ -dis -  h / num , 3 * h / num , 0 ] , color = GREEN )
        tlx = Line([ -dis - 2.4 * h / num , - 4 * h / num , 0 ], [ -dis - 0.8 * h / num , 4 * h / num , 0 ], color = YELLOW)
        tly = Line([ -dis + 2.5 * h / num , - 4 * h / num , 0 ], [ -dis - 1.5 * h / num , 4 * h / num , 0 ], color = PINK)
        inte = Tex(r"$\tilde{\ell} (x)$", r" $\cap $ ", r"$\tilde{\ell}(y)$", r" $ = \{ $ ", r"$z$", r" $\}$", font_size = 34).shift(4*RIGHT + 0.5 *UP)
        zla = Tex(r"$ \exp_p ( z ) \in A$", font_size = 34).next_to(inte, DOWN)
        lxli = Line([0,0,0], [0.5,0,0], color = YELLOW).next_to(inte[0], DOWN).shift(0.15*UP)
        lyli = Line([0,0,0], [0.5,0,0], color = PINK).next_to(inte[2], DOWN).shift(0.15*UP)
        zli = Line([0,0,0], [0.3,0,0], color = GREEN).next_to(inte[4], DOWN).shift(0.15*UP)
        self.play(Create(grid), Create(tlx), Create(tly), Create(x), Create(y))
        self.wait()
        self.play(Create(z), Create(inte), Create(lxli), Create(lyli), Create(zli))
        self.wait()
        self.play(Create(zla))
        self.wait()
        









class para2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/6)
        p = Tex(r"Claim", r" : The lines ", r"$\ell ( \exp ( x ) )$", r" and ", r"$\ell (\exp_p (y))$", r" are parallel", font_size = 34).shift( 2.5 * UP )
        p1 = Tex(r"Let ", r"$\tilde{a}$", r" $\in  \tilde{\ell }  (x)$, ", r"$\tilde{b}$", r" $\in \tilde{\ell}(y)$, and ", r"$a$", r" $= \exp_p (\tilde{a})$,  ", r"$b$" ,r" $ = \exp_p (\tilde{b})$.", r" $\Rightarrow$ ", r" $d ( a, b ) \leq d (\tilde{a}, \tilde{b})$", font_size = 34).next_to(p, DOWN)
        p11 = VGroup(p1[0], p1[1], p1[2], p1[3], p1[4],  p1[5],  p1[6],  p1[7],  p1[8])
        p2 = Tex(r"If  $\ell ( \exp _p (x) )$ and $\ell (\exp_p (y )) $ are not parallel, and $\tilde{a}, \tilde{b} \to \infty$ together, then ", font_size = 34).next_to(p1, DOWN)
        p3 = Tex(r"$d ( \tilde{a}, \tilde{b} )$ remains constant, ", font_size = 34).next_to(p2, DOWN)
        p4 = Tex(r" $d (a,b) \to \infty $", font_size = 34).next_to(p3, RIGHT)
        p34 = VGroup(p3,p4)
        p1.shift(0.35*LEFT)
        l = p1.get_left()[0]
        p2.shift([p2.get_left()[0]- l ] *LEFT)
        p34.shift([p34.get_center()[0]] *LEFT)
        cli = Line([0,0,0], [0.4,0,0]).next_to(p[0], DOWN).shift(0.15*UP)
        lxli = Line([0,0,0], [0.7,0,0], color = YELLOW).next_to(p[2], DOWN).shift(0.15*UP)
        lyli = Line([0,0,0], [0.7,0,0], color = PINK).next_to(p[4], DOWN).shift(0.15*UP)
        tali = Line([0,0,0], [0.3,0,0], color = RED).next_to(p1[1], DOWN).shift(0.15*UP)
        tbli = Line([0,0,0], [0.3,0,0], color = TEAL).next_to(p1[3], DOWN).shift(0.15*UP)
        ali = Line([0,0,0], [0.3,0,0], color = RED).next_to(p1[5], DOWN).shift(0.15*UP)
        bli = Line([0,0,0], [0.3,0,0], color = TEAL).next_to(p1[7], DOWN).shift(0.15*UP)
        box1 = SurroundingRectangle(p1[10], color = RED , buff = 0.1 )
        box2 = SurroundingRectangle(p34, color = RED , buff = 0.1 )
        h = - 2
        op = 0.2
        y = 1.7
        xdi = 1
        vle = 1.7
        s = Surface(
            lambda u, v:  [ u - y , v +  y * sqrt(3) ,  u * v / 4 + h ],
            u_range = [ - 2 , 2 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        tps = Surface(
            lambda u, v:  [ u + y  , v - y * sqrt(3) , h ],
            u_range = [ - 2 , 2 ],
            v_range = [ - 2 , 2 ],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = op, 
            resolution = [ 16, 16]
        )
        xx = 1.5
        tlx = Line( [ - xx + y , -2 - y * sqrt(3) , h ] , [ -xx + y , 2 - y * sqrt(3) , h ] , color = YELLOW )
        tly = Line( [   xx + y , -2 - y * sqrt(3) , h ] , [  xx + y , 2 - y * sqrt(3) , h ] , color = PINK )
        tlx2 = Line( [ - xx + y , -1.5 - y * sqrt(3) , h ] , [ -xx + y , 1.5 - y * sqrt(3) , h ] , color = YELLOW )
        tly2 = Line( [   xx + y , -1.5 - y * sqrt(3) , h ] , [  xx + y , 1.5 - y * sqrt(3) , h ] , color = PINK )
        lx = Line( [ - xx - y , -2 + y * sqrt(3) , h - xx * ( - 2 ) / 4 ] , [ -xx - y , 2 + y * sqrt(3) , h - xx * ( 2 ) / 4] , color = YELLOW )
        ly = Line( [   xx - y , -2 + y * sqrt(3) , h + xx * ( - 2 ) / 4 ] , [  xx - y , 2 + y * sqrt(3) , h + xx * ( 2 ) / 4] , color = PINK )
        lx2 = Line( [ - xx - y , -1.5 + y * sqrt(3) , h - xx * ( - 1.5 ) / 4 ] , [ -xx - y , 1.5 + y * sqrt(3) , h - xx * ( 1.5 ) / 4] , color = YELLOW )
        ly2 = Line( [   xx - y , -1.5 + y * sqrt(3) , h + xx * ( - 1.5 ) / 4 ] , [  xx - y , 1.5 + y * sqrt(3) , h + xx * ( 1.5 ) / 4] , color = PINK )
        po = Sphere(radius = 0.05)
        z = po.copy()
        ta = Sphere(radius = 0.1)
        tb = ta.copy()
        a = ta.copy()
        b = ta.copy()
        po.set_color(ORANGE).shift( [ - y , y * sqrt(3) , h ] )
        z.set_color(ORANGE).shift( [ y , - y * sqrt(3) , h ] )
        ta.set_color(RED).shift( [ -xx +  y ,  - 1.5  - y * sqrt(3) , h ] )
        a.set_color(RED).shift( [ - xx -y , -1.5 +  y * sqrt(3) ,  h + xx * (1.5) / 4 ] )
        tb.set_color(TEAL).shift( [ xx + y , -1.5 - y * sqrt(3) , h ] )
        b.set_color(TEAL).shift( [ xx - y , -1.5 + y *sqrt(3) , h - xx * (1.5) / 4 ] )
        ar = Arrow([-1,-1.4,0], [1,-1.4,0], max_tip_length_to_length_ratio = 0.1 , stroke_width = 3)
        arla = Tex(r"$\text{exp}_p$", font_size = 34).next_to(ar, UP).shift(0.15*DOWN)
        allfig = VGroup(s, tps, po, z, ar, arla, tlx, tly, lx, ly)
        self.add_fixed_in_frame_mobjects(p, p1, p2, p3, p4, cli, lxli, lyli, ar, arla, tali, tbli, ali, bli, box1, box2)
        self.remove(p, p1, p2, p3, p4, cli, lxli, lyli, ar, arla, tali, tbli, ali, bli, box1, box2)
        self.play(Create(allfig), Create(p), Create(cli), Create(lxli), Create(lyli))
        self.wait()
        self.play(Create(p11), Create(ali), Create(bli), Create(tali), Create(tbli), Create(a), Create(b), Create(ta), Create(tb))
        self.wait()
        self.play(Create(p1[9]), Create(p1[10]))
        self.wait()
        self.play(Create(p2), Create(p34))
        self.wait()
        self.play(MoveAlongPath(a,lx2), MoveAlongPath(b,ly2), MoveAlongPath(ta, tlx2), MoveAlongPath(tb, tly2), run_time = 3 )
        self.wait()
        self.play(Create(box1), Create(box2))
        self.wait()



class acor2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=80 * DEGREES, theta=  PI/3)
        t = Tex(r"Corollary", font_size = 34, color = BLUE).to_edge(UL).shift(0.8*DOWN+0.3*RIGHT)
        t1 = Tex(r"If ", r"$\Sigma$", r"  $\subset \mathbb{R}^3$ is a proper surface with zero Gauss curvature, there is a line ", r"$L$", r" $\subset \mathbb{R}^3$, ", font_size = 34).next_to(t, DOWN)
        t2 = Tex(r"such that for all ", r"$x$", r" $\in $ ", r"$A$", r", the line ", r"$\ell (p)$", r" $ \subset $ ", r"$A$",  r" is parallel to ",  r"$L$", font_size = 34).next_to(t1, DOWN)
        l = t.get_left()[0]
        t1.shift([t1.get_left()[0]- l ] *LEFT)
        t2.shift([t2.get_left()[0]- l ] *LEFT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(t1[1], DOWN).shift(0.15*UP)
        lli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(t1[3], DOWN).shift(0.15*UP)
        xli = Line([0,0,0], [0.4,0,0], color = RED).next_to(t2[1], DOWN).shift(0.15*UP)
        ali = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(t2[3], DOWN).shift(0.15*UP)
        lxli = Line([0,0,0], [0.6,0,0], color = YELLOW).next_to(t2[5], DOWN).shift(0.15*UP)
        ali2 = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(t2[7], DOWN).shift(0.15*UP)
        lli2 = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(t2[9], DOWN).shift(0.15*UP)
        op = 0.2
        h = - 1.8
        yra = 3.5
        x = Sphere(radius = 0.07).set_color(RED).shift( [ 1.5 * 2.125 , 1.5 , 6 * bump2(2.125) + h ] )
        lx = ParametricFunction(lambda t : [ 1.5 * 2.125 , t , 6 * bump2(2.125) + h ] , t_range = [-yra, yra], color = YELLOW)
        ll = []
        llg = VGroup()
        for i in range(0, 3):
            ll.append([])
            for j in range(0, 7):
                ll[i].append( ParametricFunction(lambda t : [ 1.5 * (-2.5 + 2 * i + ( j + 1 ) / 8 ) , t , 6 * bump2(-2.5 + 2 * i + ( j + 1 ) / 8 ) + h ] , t_range = [-yra, yra], color = YELLOW) )
                llg.add(ll[i][j])
        s1 = Surface(
            lambda u, v:  [ 1.5 * u , v , 6 * bump2(u) + h ],
            u_range = [ -2.5 , -1.5 ],
            v_range = [ -yra , yra ],
            checkerboard_colors = [BLUE , BLUE],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s2 = Surface(
            lambda u, v:  [ 1.5 * u , v , 6 * bump2(u) + h ],
            u_range = [ -1.5 , -0.5 ],
            v_range = [ -yra , yra ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s3 = Surface(
            lambda u, v:  [ 1.5 * u , v , 6 * bump2(u) + h ],
            u_range = [ -0.5 , 0.5 ],
            v_range = [ -yra , yra ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s4 = Surface(
            lambda u, v:  [ 1.5 * u , v , 6 * bump2(u) + h ],
            u_range = [ 0.5 , 1.5 ],
            v_range = [ -yra , yra ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s5 = Surface(
            lambda u, v:  [ 1.5 * u , v , 6 * bump2(u) + h ],
            u_range = [ 1.5 , 2.5 ],
            v_range = [ -yra , yra ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s = VGroup(s1, s2, s3, s4, s5)
        figures = VGroup( s, x  )
        self.add_fixed_in_frame_mobjects(t, t1, t2, sli, lli, xli, ali, lxli, ali2, lli2)
        self.remove(t, t1, t2, sli, lli, xli, ali, lxli, ali2, lli2)
        self.play(Create(t), Create(t1), Create(t2), Create(sli), Create(lli), Create(s), Create(xli), Create(ali), Create(lxli), Create(ali2), Create(lli2))
        self.wait()
        self.play(Create(x))
        self.wait()
        self.play(Create(lx))
        self.wait()
        self.play(Create(llg))
        self.wait()



class acor3(Scene):
    def construct(self):
        t = Tex(r"Corollary", font_size = 34, color = BLUE).to_edge(UL).shift(0.8*DOWN+0.3*RIGHT)
        t1 = Tex(r"If ", r"$\Sigma$", r"  $\subset \mathbb{R}^3$ is a proper surface with zero Gauss curvature, there is a line ", r"$L$", r" $\subset \mathbb{R}^3$, ", font_size = 34).next_to(t, DOWN)
        t2 = Tex(r"such that for all ", r"$x$", r" $\in $ ", r"$A$", r", the line ", r"$\ell (p)$", r" $ \subset $ ", r"$A$",  r" is parallel to ",  r"$L$", font_size = 34).next_to(t1, DOWN)
        c0 = Tex(r"Corollary", font_size = 34, color = BLUE).next_to(t2, DOWN)
        c1 = Tex(r"For each ", r"$y$", r" $ \in \partial A$, the line parallel to $L$ passing through $y$ is contained in $\partial A$.", font_size = 34).next_to(c0, DOWN)
        pr0 = Tex(r"Proof", font_size = 34, color = BLUE).next_to(c1, DOWN)
        pr1 = Tex(r"Take ", r"$x_n$", r" $ \in A$ with $x_n \to y$. Then the lines ", r"$\ell (x_n)$", r" converge to ", r"the desired line.", font_size = 34).next_to(pr0, DOWN)
        l = t.get_left()[0]
        t1.shift([t1.get_left()[0]- l ] *LEFT)
        t2.shift([t2.get_left()[0]- l ] *LEFT)
        c0.shift([c0.get_left()[0]- l ] *LEFT)
        c1.shift([c1.get_left()[0]- l ] *LEFT)
        pr0.shift([pr0.get_left()[0]- l ] *LEFT)
        pr1.shift([pr1.get_left()[0]- l ] *LEFT)
        sli = Line([0,0,0], [0.4,0,0], color = PURPLE).next_to(t1[1], DOWN).shift(0.15*UP)
        lli = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(t1[3], DOWN).shift(0.15*UP)
        xli = Line([0,0,0], [0.4,0,0], color = RED).next_to(t2[1], DOWN).shift(0.15*UP)
        ali = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(t2[3], DOWN).shift(0.15*UP)
        lxli = Line([0,0,0], [0.6,0,0], color = YELLOW).next_to(t2[5], DOWN).shift(0.15*UP)
        ali2 = Line([0,0,0], [0.4,0,0], color = BLUE).next_to(t2[7], DOWN).shift(0.15*UP)
        lli2 = Line([0,0,0], [0.4,0,0], color = YELLOW).next_to(t2[9], DOWN).shift(0.15*UP)
        yyli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(c1[1], DOWN).shift(0.15*UP)
        xnli = Line([0,0,0], [0.3,0,0], color = RED).next_to(pr1[1], DOWN).shift(0.15*UP)
        lxnli = Line([0,0,0], [0.7,0,0], color = YELLOW).next_to(pr1[3], DOWN).shift(0.15*UP)
        lyli = Line([0,0,0], [0.7,0,0], color = GREEN).next_to(pr1[5], DOWN).shift(0.15*UP)
        a = 25
        el = 4
        h = -2
        sl = 10
        y = Dot( [ 0 , h , 0] , color = ORANGE)
        ly = Line( [-el , h  - el / sl, 0 ], [el, h + el /sl , 0], color = GREEN )
        x = []
        lx = []
        xg = VGroup()
        lxg = VGroup()
        for i in range(1,5):
            x.append( Dot([ i / 10 , h - i / 8 - i ** 2 / a, 0 ], color = RED))
            lx.append(Line( [ - el + i /10 , h  -i / 8 - i ** 2 / a - el / sl , 0 ] , [ el + i /10 , h -i / 8 - i ** 2 / a + el /sl , 0 ] , color = YELLOW ) )
            xg.add(x[i-1])
            lxg.add(lx[i-1])
        self.add(t, t1, t2, sli, lli, xli, ali, lxli, ali2, lli2)
        self.wait()
        self.play(Create(t), Create(t1), Create(t2), Create(sli), Create(lli), Create(xli), Create(ali), Create(lxli), Create(ali2), Create(lli2))
        self.wait()
        self.play(Create(c0), Create(c1), Create(yyli))
        self.wait()
        self.play(Create(y))
        self.wait()
        self.play(Create(pr0), Create(pr1), Create(xg), Create(lxnli), Create(lyli), Create(xnli))
        self.wait()
        self.play(Create(lxg))
        self.wait()
        self.play(Create(ly))
        self.wait()


class lifts2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=78 * DEGREES, theta=  PI/3)
        le = Tex(r"$A \cup \partial A$ is a union of parallel lines ", r"$\Rightarrow$  $\Sigma \backslash  A  $ is a  union of ", r"strips", font_size = 34).shift(2.5*UP)
        stli = Line([0,0,0] , [0.6, 0, 0], color = GOLD   ).next_to(le[2], DOWN).shift(0.15*UP)
        stli1 = Line([0,0,0] , [0.2, 0, 0], color = YELLOW).next_to(le[2], DOWN).shift(0.15*UP + 0.3 * LEFT )
        stli2 = Line([0,0,0] , [0.2, 0, 0], color = YELLOW).next_to(le[2], DOWN).shift(0.15*UP + 0.3 * RIGHT)
        stlii = VGroup(stli, stli1, stli2)
        op = 0.3
        h = -0.5
        #hh = 0.7
        sc1 = 1.4
        sc2 = 6
        w = sc1 * 1.5
        #p = Sphere(radius = 0.05)
        #p.set_color(ORANGE).shift([0,0,h + sc2 * bump2(0)])
        #z = Sphere(radius = 0.05)
        #z.set_color(ORANGE).shift([0,0,hh])
        #tps1  = Surface(lambda u,v: [w * u,w * v,hh], u_range = [-2.5, 0.5], v_range = [-3,3], checkerboard_colors = [GREEN, GREEN], fill_opacity = op, resolution = [12,20])
        #tps2  = Surface(lambda u,v: [w * u,w * v,hh], u_range = [ 0.5, 1.5], v_range = [-3,3], checkerboard_colors = [GREEN, GREEN], fill_opacity = op, resolution = [4,20])
        #tps3  = Surface(lambda u,v: [w * u,w * v,hh], u_range = [ 1.5, 2.5], v_range = [-3,3], checkerboard_colors = [GREEN, GREEN], fill_opacity = op, resolution = [4,20])
        #tps22 = Surface(lambda u,v: [w * u,w * v,hh], u_range = [ 0.5, 1.5], v_range = [-3,3], checkerboard_colors = [GOLD, GOLD], fill_opacity = 2*op, resolution = [4,20])
        #tps = VGroup(tps1, tps2, tps3)
        s0 = Surface(
            lambda u, v:  [ sc1 * u , w *  v , sc2 * bump2(u) + h ],
            u_range = [ -3 , -2.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op, 
            resolution = [2, 16]
        )
        s1 = Surface(
            lambda u, v:  [ sc1 * u , w *  v , sc2 * bump2(u) + h ],
            u_range = [ -2.5 , -1.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s2 = Surface(
            lambda u, v:  [ sc1 * u ,w *  v , sc2 * bump2(u) + h ],
            u_range = [ -1.5 , -0.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s44 = Surface(
            lambda u, v:  [ sc1 * u , w *  v , sc2 * bump2(u) + h ],
            u_range = [ 0.5 , 1.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [ GOLD, GOLD ],
            fill_opacity = 2*op, 
            resolution = [4, 16]
        )
        s3 = Surface(
            lambda u, v:  [ sc1 * u ,w *  v , sc2 * bump2(u) + h ],
            u_range = [ -0.5 , 0.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s4 = Surface(
            lambda u, v:  [ sc1 * u ,w *  v , sc2 * bump2(u) + h ],
            u_range = [ 0.5 , 1.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s5 = Surface(
            lambda u, v:  [ sc1 * u , w * v , sc2 * bump2(u) + h ],
            u_range = [ 1.5 , 2.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s6 = Surface(
            lambda u, v:  [ sc1 * u , w *  v , sc2 * bump2(u) + h ],
            u_range = [ 2.5 , 3 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op, 
            resolution = [2, 16]
        )
        s = VGroup(s0, s1, s2, s3, s4, s5, s6)
        ll  = []
        #llt = []
        llg = VGroup()
        #lltg = VGroup()
        llg2 = VGroup()
        #lltg2 = VGroup()
        for i in range(0, 3):
            ll.append([])
        #    llt.append([])
        #    for j in range(0, 9):
        #        llt[i].append( ParametricFunction(lambda t : [ w * (  -2.5 + 2 * i + j / 8 ) , w * t , hh ] , t_range = [-3, 3], color = YELLOW) )
        #        lltg.add(llt[i][j])
        #        if (i != 1 or j !=8) and (i != 2 or j != 0 ):
        #            lltg2.add(llt[i][j])
            for j in range(0,17):
                ll[i].append(  ParametricFunction(lambda t : [ sc1 * (-2.5 + 2 * i + j / 16 ) , w * t , sc2 * bump2(-2.5 + 2 * i + j / 16 ) + h ] , t_range = [-3, 3], color = YELLOW) )
                llg.add(ll[i][j])
                if (i != 1 or j !=16) and (i != 2 or j != 0 ):
                    llg2.add( ll[i][j] )
        sothers = VGroup(s0,s1,s2,s3,s5,s6)
        self.add_fixed_in_frame_mobjects(le, stlii)
        self.remove(le, stlii)
        self.play(Create(s))
        self.wait()
        self.play(Create(llg), Write(le[0]))
        self.wait()
        self.play(FadeOut(llg2), Transform(s4, s44), Create(stlii), Create(le[1]), Create(le[2]), FadeOut(sothers))
        self.wait()





class thumbnail(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=74 * DEGREES, theta=  PI/6 + PI /2)
        op = 0.3
        hh = 0.6
        h = -1.5
        p = Sphere(radius = 0.05)
        p.set_color(ORANGE).shift([0,0,h + 3 * bump2(0)])
        z = Sphere(radius = 0.05)
        z.set_color(ORANGE).shift([0,0,hh])
        x = Sphere(radius = 0.05)
        x.set_color(RED).shift([ - 2.3 , 1 ,hh])
        exx = Sphere(radius = 0.05)
        exx.set_color(RED).shift([ - 2.2 * 0.8 , 1 , 3 * bump2( 2.2 ) + h ])
        y = Sphere(radius = 0.05)
        y.set_color(TEAL).shift([ - 2.75 , -0.6 ,hh])
        exy = Sphere(radius = 0.05)
        exy.set_color(TEAL).shift([ - 2.4 * 0.8 , -0.6 , 3 * bump2( 2.4 ) + h ])
        lx = ParametricFunction(lambda t : [  - 2.2 * 0.8 , t , 3 * bump2(2.2)  + h ] , t_range = [-3,3], color = YELLOW)
        tlx = ParametricFunction(lambda t : [  - 2.3 , t , hh ] , t_range = [-3, 3], color = YELLOW)
        yline = ParametricFunction(lambda t : [  - 2.4 * 0.8 , t , 3 * bump2(2.4)  + h ] , t_range = [-3,3], color = PINK)
        yline2 = ParametricFunction(lambda t : [  - 2.75 , t , hh ] , t_range = [-3, 3], color = PINK)
        tps = Surface(lambda u,v: [u,v,hh], u_range = [-3,3], v_range = [-3,3], checkerboard_colors = [GREEN, GREEN], fill_opacity = op, resolution = [16,16])
        s0 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ -3 , -2.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op, 
            resolution = [2, 16]
        )
        s1 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ -2.5 , -1.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s2 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ -1.5 , -0.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s3 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ -0.5 , 0.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s4 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ 0.5 , 1.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B,PURPLE_B],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s5 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ 1.5 , 2.5 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = op, 
            resolution = [4, 16]
        )
        s6 = Surface(
            lambda u, v:  [ u * 0.8 , v , 3 * bump2(u) + h ],
            u_range = [ 2.5 , 3 ],
            v_range = [ -3 , 3 ],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = op, 
            resolution = [2, 16]
        )
        s = VGroup(s1, s2, s3, s4, s5, p, z, tps, s0, s6, exx, exy, yline, yline2, lx, tlx, x, y)
        self.add(s)
        self.wait()

