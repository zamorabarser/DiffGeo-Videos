from asyncio import threads
from this import d
from tkinter import E
from manim import *
from numpy import sqrt
import math



class ti(Scene):
    def construct(self):
        t1 = Text("Semigeodesic charts", font_size=60)
        self.play(Write(t1))
        self.wait()        


class geextt(Scene):
    def construct(self):
        gee = Tex(r"Theorem", color = BLUE,  font_size = 34).to_edge(UL).shift(0.4*DOWN)
        l = gee.get_left()[0]
        gee1 = Tex(r"For a parametrization $\phi : U \to \Sigma$, and $p = \phi (0)$, there is an open set $A \subset \mathbb{R}^4$", font_size = 34).next_to(gee, DOWN)
        gee1.shift( ( l - gee1.get_left()[0])*RIGHT)
        gee2 = Tex(r"with $0 \in A$ and $\varepsilon > 0 $ such that for all $(u,v,a,b) \in A$, there is a geodesic ", font_size = 34).next_to(gee1, DOWN)
        gee2.shift( ( l - gee2.get_left()[0])*RIGHT)
        gee3 = Tex(r"$\gamma_{u,v,a,b} : ( - \varepsilon , \varepsilon ) \to \Sigma $ with $\gamma_{u,v,a,b}(0) = \phi (u,v)$, $\gamma_{u,v,a,b}^{\prime}(0) = a \phi_u + b \phi_v $,", font_size = 34).next_to(gee2, DOWN)
        gee3.shift( ( l - gee3.get_left()[0])*RIGHT)
        gee4 = Tex(r"and the point $\gamma_{u,v,a,b}(t)$ depends smoothly on $(u,v,a,b,t) \in A \times (-\varepsilon , \varepsilon)$.", font_size = 34).next_to(gee3, DOWN)
        gee4.shift( ( l - gee4.get_left()[0])*RIGHT)
        self.add(gee, gee1, gee2, gee3, gee4)
        self.wait()
        co = Tex(r"Corollary", color = BLUE,  font_size = 34).next_to(gee4, DOWN).shift(0.2*DOWN)
        co1 = Tex(r"For each $p \in \Sigma$, there is a ball $B \subset T_p \Sigma$ and a map $ \text{exp}_p : B \to \Sigma$", font_size = 34).next_to(co, DOWN)
        co2 = Tex(r"such that $\text{exp}_p (v) = \gamma _v (1)$, where $ \gamma _v $ is a geodesic with $\gamma_v(0) = p$, $\gamma_v^{\prime}(0) = v$.", font_size = 34).next_to(co1, DOWN)
        co3 = Tex(r"Moreover, $\text{exp}_p$ is smooth, and $d_p \text{exp} = \text{Id}_{T_p \Sigma}$.", font_size = 34).next_to(co2, DOWN)
        co4 = Tex(r"If $B$ is small enough, $s = \text{exp}_p : B \to \Sigma$ is a chart.", font_size = 34).next_to(co3, DOWN)
        co.shift( ( l - co.get_left()[0])*RIGHT)
        co1.shift( ( l - co1.get_left()[0])*RIGHT)
        co2.shift( ( l - co2.get_left()[0])*RIGHT)
        co3.shift( ( l - co3.get_left()[0])*RIGHT)
        co4.shift( ( l - co4.get_left()[0])*RIGHT)
        self.play(Write(co), Write(co1), Write(co2), Write(co3))
        self.wait()
        self.play(Write(co4))
        self.wait()



class polar(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/6)
        sdef = Tex(r"$s: $ ", r"$B$", r" $\to$ ", r"$\Sigma$", font_size = 34).shift(2.4*UP)
        srt = Tex(r"$s_r$", r", ", r"$s_{\theta}$", r" $: B \backslash \{ 0 \} \to \mathbb{R}^3$", font_size = 34).next_to(sdef, DOWN).shift(0.3*DOWN)
        gle = Tex(r"Gauss Lemma", font_size = 34, color = BLUE).to_edge(UL).shift(RIGHT + 2.2 * DOWN)
        gle1 = Tex(r"The vector fields $s_r$ and $s_{\theta}$ are orthogonal.", font_size = 34).next_to(gle, DOWN)
        l = gle.get_left()[0]
        gle1.shift((l-gle1.get_left()[0])*RIGHT)
        bli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(sdef[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(sdef[3], DOWN).shift(0.15*UP)
        srli = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(srt[0], DOWN).shift(0.15*UP)
        stli = Line([0,0,0], [0.3,0,0], color = YELLOW).next_to(srt[2], DOWN).shift(0.15*UP)
        ts = Surface(
            lambda u, v: [ u , v , 0],
            u_range = [-2.6, 2.6],
            v_range = [-2.6, 2.6],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.3,
            resolution = 24
        )
        b = Surface(
            lambda u, v: [1.5*u * np.sin(v),1.5* u* np.cos(v) , 0],
            u_range = [0, PI/2],
            v_range = [0, TAU],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.6,
            resolution = 12
        )
        pb = Surface(
            lambda u, v: [1.5*np.sin(u) * np.sin(v),1.5* np.sin(u)* np.cos(v) ,1.5* np.cos(u) - 1.5],
            u_range = [0, PI/2],
            v_range = [0, TAU],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.6,
            resolution = 12
        )
        s = Surface(
            lambda u, v: [1.5*np.sin(u) * np.sin(v), 1.5*np.sin(u)* np.cos(v) , 1.5*np.cos(u)- 1.5],
            u_range = [0, PI],
            v_range = [0, TAU],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.2,
            resolution = 24
        )
        sr = []
        sr0 = []
        st = []
        st0 = []
        vsr = VGroup()
        vst = VGroup()
        a = 8
        c = 12
        def x(i,j):
            return 1.5*np.sin(PI*j/c)*np.cos(2*PI*i/a) - 1.5
        def y(i,j):
            return 1.5*np.sin(PI*j/c)*np.sin(2*PI*i/a) + sqrt(27)/2
        def x0(i,j):
            return 1.5*PI*(j/c)*np.cos(2*PI*i/a) + 1.5
        def y0(i,j):
            return 1.5*(PI*j/c)*np.sin(2*PI*i/a) - sqrt(27)/2
        def x00(i,j):
            return 1.5*np.cos(2*PI*i/a) 
        def y00(i,j):
            return 1.5*np.sin(2*PI*i/a) 
        def xp(i,j):
            return - 1.6*np.sin(PI*j/c)*np.sin(2*PI*i/a)
        def yp(i,j):
            return 1.6*np.sin(PI*j/c)*np.cos(2*PI*i/a)
        def xp0(i,j):
            return - 1.6*(PI*j/c)*np.sin(2*PI*i/a)
        def yp0(i,j):
            return 1.6*(PI*j/c)*np.cos(2*PI*i/a)
        def xpp(i,j):
            return np.cos(PI*j/c)*np.cos(2*PI*i/a) 
        def ypp(i,j):
            return np.cos(PI*j/c)*np.sin(2*PI*i/a)
        def zpp(j):
            return - np.sin(PI*j/c)
        def z(j):
            return 1.5*np.cos(PI*j/c) - 1.5 
        for j in range(1,6):
            sr.append([])
            sr0.append([])
            st.append([])
            st0.append([])
            for i in range(0,8):
                st[j-1].append( Arrow( [x(i,j) - 0.15*xp(i,j), y(i,j)- 0.15*yp(i,j), z(j)], [x(i,j) + xp(i,j), y(i,j)+ yp(i,j),z(j)] , color = YELLOW , max_tip_length_to_length_ratio = 0.2) )
                st0[j-1].append( Arrow( [x0(i,j) - 0.08*xp0(i,j) , y0(i,j)- 0.08*yp0(i,j), -1.5], [x0(i,j) + xp0(i,j), y0(i,j)+ yp0(i,j) ,-1.5] , color = YELLOW , max_tip_length_to_length_ratio = 0.2) )
                sr[j-1].append( Arrow( [x(i,j) - 0.15*xpp(i,j), y(i,j)- 0.15*ypp(i,j), z(j) - 0.15*zpp(j)], [x(i,j) + xpp(i,j), y(i,j)+ ypp(i,j),z(j) + zpp(j)] , color = BLUE , max_tip_length_to_length_ratio = 0.2) )
                sr0[j-1].append( Arrow( [x0(i,j) , y0(i,j) , -1.5], [x0(i,j) + x00(i,j) , y0(i,j) + y00(i,j), -1.5 ] , color = BLUE ) )
                vsr.add(sr[j-1][i], sr0[j-1][i])
                vst.add(st[j-1][i], st0[j-1][i])
        srg = sr[3][1].copy()
        srg0 = sr0[3][1].copy()
        stg = st[3][1].copy()
        stg0 = st0[3][1].copy()
        self.add_fixed_in_frame_mobjects(sdef, srt, sli, bli, srli, stli, gle, gle1)
        self.remove(sdef, srt, sli, bli, srli, stli, gle, gle1)
        self.play(Create(s), Create(ts), Create(sdef), Create(sli), Create(bli))
        self.wait()
        self.play(Create(b))
        self.wait()
        b1 = b.copy()
        self.play(Transform(b1, pb))
        self.wait()
        self.remove(pb, b1)
        self.add(pb)
        self.play(pb.animate.shift([-1.5, sqrt(27)/2, 0]), s.animate.shift([-1.5, sqrt(27)/2, 0]), ts.animate.shift([1.5, - sqrt(27)/2, -1.5]), b.animate.shift([1.5, - sqrt(27)/2, -1.5]) )
        self.wait()
        self.play(Create(srt), Create(srli), Create(stli))
        self.wait()
        self.play(Create(vsr))
        self.wait()
        self.play(FadeOut(vsr))
        self.wait()
        self.play(Create(vst))
        self.wait()
        self.add(stg, stg0)
        len = 0.3
        self.play(FadeOut(vst), Create(srg), Create(srg0), srt.animate.shift(len*UP), srli.animate.shift(len*UP), stli.animate.shift(len*UP), Create(gle), Create(gle1))
        self.wait()
        #self.begin_ambient_camera_rotation(rate=PI/8)
        




class gl(Scene):
    def construct(self):
        sdef = Tex(r"$s: $ ", r"$B$", r" $\to$ ", r"$\Sigma$", font_size = 34).shift(2.4*UP)
        srt = Tex(r"$s_r$", r", ", r"$s_{\theta}$", r" $: B \backslash \{ 0 \} \to \mathbb{R}^3$", font_size = 34).next_to(sdef, DOWN).shift(0.3*DOWN)
        gle = Tex(r"Gauss Lemma", font_size = 34, color = BLUE).to_edge(UL).shift(RIGHT + 2.2 * DOWN)
        gle1 = Tex(r"The vector fields $s_r$ and $s_{\theta}$ are orthogonal.", font_size = 34).next_to(gle, DOWN)
        glp = Tex(r"Proof", font_size = 34, color = BLUE).next_to(gle1, DOWN)
        glp1 = Tex(r"$\vert s_r \vert \equiv 1 $", r" $ \Rightarrow $ $  0 =  \partial_{\theta} \vert s_r \vert ^2 $", r" $ = 2 s_r \cdot s_{r\theta}$.", font_size = 34).next_to(glp, DOWN)
        glp2 = Tex(r"$\partial _r (s_r \cdot s_{\theta}) $ ", r"$ = $ ", r"$s_{rr} \cdot s_{\theta}$", r" $ + $ ", r"$s_r \cdot s_{r \theta}$ ", r"$= 0$.", font_size = 34).next_to(glp1, DOWN)
        glp3 = Tex(r"$ \lim_{r \to 0} ( s_r \cdot s_{\theta} )  = 0$", r" $ \Rightarrow$ $s_r \cdot s_{\theta} \equiv 0 $.", font_size = 34).next_to(glp2, DOWN)
        l = gle.get_left()[0]
        gle1.shift((l-gle1.get_left()[0])*RIGHT)
        glp.shift((l-glp.get_left()[0])*RIGHT)
        glp1.shift((l-glp1.get_left()[0])*RIGHT)
        glp2.shift((l-glp2.get_left()[0])*RIGHT)
        glp3.shift((l-glp3.get_left()[0])*RIGHT)
        len = 0.3
        bli = Line([0,0,0], [0.3,0,0], color = ORANGE).next_to(sdef[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.3,0,0], color = PURPLE_B).next_to(sdef[3], DOWN).shift(0.15*UP)
        srli = Line([0,0,0], [0.3,0,0], color = BLUE).next_to(srt[0], DOWN).shift(0.15*UP)
        stli = Line([0,0,0], [0.3,0,0], color = YELLOW).next_to(srt[2], DOWN).shift(0.15*UP)
        srt.shift(len*UP)
        srli.shift(len*UP)
        stli.shift(len*UP)
        zer1 = Line([0,0,0], [0.9,0.3,0], color = RED).next_to(glp2[2], DOWN).shift(0.5*UP)
        zer2 = Line([0,0,0], [0.9,0.3,0], color = RED).next_to(glp2[4], DOWN).shift(0.5*UP)
        self.add(sdef, sli, bli, srt, srli, stli, gle, gle1)
        self.wait()
        self.play(Create(glp), Create(glp1[0]))
        self.wait()
        self.play(Write(glp1[1]))
        self.wait()
        self.play(Write(glp1[2]))
        self.wait()
        self.play(Write(glp2[0]))
        self.wait()
        self.play(Write(glp2[1]), Write(glp2[2]), Write(glp2[3]), Write(glp2[4]))
        self.wait()
        self.play(Create(zer1), Create(zer2))
        self.wait()
        self.play(Write(glp2[5]))
        self.wait()
        self.play(Write(glp3[0]))
        self.wait()
        self.play(Write(glp3[1]))
        self.wait()
        #self.begin_ambient_camera_rotation(rate=PI/8)
        

class sgc(Scene):
    def construct(self):
        d0 = Tex(r"Definition", font_size = 34, color = BLUE).to_edge(UL).shift(RIGHT +  DOWN)
        d1 = Tex(r"A chart $s : U \to \Sigma$ is called a ", r"semigeodesic chart", r" if ", font_size = 34).next_to(d0, DOWN)
        d2 = Tex(r"$\bullet$ For each $v$, the map $u \mapsto s(u,v)$ is a unit speed geodesic.", font_size = 34).next_to(d1, DOWN)
        d3 = Tex(r"$\bullet$ The vectors $s_u$ and $s_v$ are orthogonal.", font_size = 34).next_to(d2, DOWN)
        d1[1].set_color(BLUE)
        e0 = Tex(r"Exercise", font_size = 34, color = BLUE).next_to(d3, DOWN).shift(0.5*DOWN)
        e1 = Tex(r"For each $p \in \Sigma $ there is a semigeodesic chart", font_size = 34).next_to(e0, DOWN)
        e2 = Tex(r"$s : U \to \Sigma$ with $p \in s(U)$.", font_size = 34).next_to(e1, DOWN)
        l = d0.get_left()[0]
        d1.shift((l-d1.get_left()[0])*RIGHT)
        d2.shift((l-d2.get_left()[0])*RIGHT)
        d3.shift((l-d3.get_left()[0])*RIGHT)
        e0.shift((l-e0.get_left()[0])*RIGHT)
        e1.shift((l-e1.get_left()[0])*RIGHT)
        e2.shift((l-e2.get_left()[0])*RIGHT)
        self.play(Write(d0), Write(d1), Write(d2), Write(d3))
        self.wait()
        self.play(Write(e0), Write(e1), Write(e2))
        self.wait()
        #self.begin_ambient_camera_rotation(rate=PI/8)
        


class gmc(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/6)
        p0 = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift( 0.5*DOWN )
        p1 = Tex(r"A unit-speed minimizing curve $\gamma : [0, \ell ] \to \Sigma $ is a geodesic.", font_size = 34).next_to(p0, DOWN)
        p2 = Tex(r"Conversely, if $\gamma : (a,b) \to \Sigma$ is a geodesic, for each $t_0 \in (a,b)$", font_size = 34).next_to(p1, DOWN)
        p3 = Tex(r"there is $\varepsilon > 0 $ such that $\gamma \vert _{[ t_0 - \varepsilon, t_0 + \varepsilon ]}$ is minimizing.", font_size = 34).next_to(p2, DOWN)
        pr0 = Tex(r"Proof", font_size = 34, color = BLUE).next_to(p3, DOWN)
        pr1 = Tex(r"Can assume ", r"$\gamma (t)$", r" $ = s(r(t), \theta (t))$ in ", r"polar coordinates around $\gamma (0)$.", font_size = 34).next_to(pr0, DOWN)
        pr2 = Tex(r"If $\gamma $ is smooth, then $\gamma ^{\prime} = r^{\prime } s_r + \theta ^{\prime} s_{\theta}$.", font_size = 34).next_to(pr1, DOWN)
        pr3 = Tex(r"length$(\gamma) = \int_0 ^{\ell} \vert \gamma ^{\prime} \vert $", r" $\geq \int_0 ^{\ell} r^{\prime} \vert s_r \vert $", r" $= \int_0^{\ell} r^{\prime} = r(\ell)$.", font_size = 34).next_to(pr2, DOWN)
        pr4 = Tex(r"If $\alpha : [0, r(\ell)]  \to \Sigma$  is given by ", r"$ \alpha (t)$", r" $ = s(t, \theta (\ell ))$, then", font_size = 34).next_to(pr3, DOWN)
        pr5 = Tex(r"length$(\alpha) = r(\ell) \leq $ length$(\gamma) \leq $ length$(\alpha)$ ", r"$\Rightarrow$ $\theta ^{\prime} \equiv 0$, $\alpha \equiv \gamma$.", font_size = 34).next_to(pr4, DOWN)
        l = p0.get_left()[0]
        p1.shift((l-p1.get_left()[0])*RIGHT)
        p2.shift((l-p2.get_left()[0])*RIGHT)
        p3.shift((l-p3.get_left()[0])*RIGHT)
        pr0.shift((l-pr0.get_left()[0])*RIGHT)
        pr1.shift((l-pr1.get_left()[0])*RIGHT)
        pr2.shift((l-pr2.get_left()[0])*RIGHT)
        pr3.shift((l-pr3.get_left()[0])*RIGHT)
        pr4.shift((l-pr4.get_left()[0])*RIGHT)
        pr5.shift((l-pr5.get_left()[0])*RIGHT)
        gamli = Line([0,0,0], [0.5, 0, 0], color = BLUE).next_to(pr1[1], DOWN).shift(0.15*UP)
        alpli = Line([0,0,0], [0.5, 0, 0], color = YELLOW).next_to(pr4[1], DOWN).shift(0.15*UP)
        ncli = Line([0,0,0], [4, 0, 0], color = ORANGE).next_to(pr1[3], DOWN).shift(0.15*UP)
        dis = 2.7
        pb = Surface(
            lambda u, v: [np.sin(u) * np.sin(v) - dis/2, np.sin(u)* np.cos(v) + dis*sqrt(3) ,  np.cos(u) - 1.5],
            u_range = [0, PI/2],
            v_range = [0, TAU],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.6,
            resolution = 12
        )
        gam = ParametricFunction(lambda t : [ np.sin(t) *np.cos( PI/4 +   np.sin( 3 * t ) /2 ) - dis/2 , np.sin(t) * np.sin( PI/4 +  np.sin( 3 * t ) /2  ) + dis * sqrt(3) , np.cos(t) - 1.5  ] , t_range = [0, PI/3], color = BLUE)
        alp = ParametricFunction(lambda t : [ np.sin(t) / sqrt(2) - dis/2 , np.sin(t) / sqrt(2) + dis * sqrt(3) , np.cos(t) -1.5 ] , t_range = [0, PI/3], color = YELLOW)
        self.add_fixed_in_frame_mobjects(p0, p1, p2, p3, pr0, pr1, pr2, pr3, pr4, pr5, gamli, alpli, ncli)
        self.remove(p0, p1, p2, p3, pr0, pr1, pr2, pr3, pr4, pr5, gamli, alpli, ncli)
        self.play(Write(p0), Write(p1), Write(p2), Write(p3))
        self.wait()
        self.play(Write(pr0), Write(pr1), Create(pb), Create(gam), Create(gamli), Create(ncli))
        self.wait()
        self.play(Write(pr2))
        self.wait()
        self.play(Write(pr3[0]))
        self.wait()
        self.play(Write(pr3[1]))
        self.wait()
        self.play(Write(pr3[2]))
        self.wait()
        self.play(Write(pr4), Create(alp), Create(alpli))
        self.wait()
        self.play(Write(pr5[0]))
        self.wait()
        self.play(Write(pr5[1]))
        self.wait()
        self.play(FadeOut(pr1), FadeOut(pr2), FadeOut(pr3), FadeOut(pr4), FadeOut(pr5), FadeOut(gam), FadeOut(alp), FadeOut(alpli), FadeOut(gamli), FadeOut(ncli), FadeOut(pb))
        self.wait()
        
        #self.begin_ambient_camera_rotation(rate=PI/8)
        


class gmc2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/6)
        p0 = Tex(r"Proposition", font_size = 34, color = BLUE).to_edge(UL).shift( 0.5*DOWN )
        p1 = Tex(r"A unit-speed minimizing curve $\gamma : [0, \ell ] \to \Sigma $ is a geodesic.", font_size = 34).next_to(p0, DOWN)
        p2 = Tex(r"Conversely, if $\gamma : (a,b) \to \Sigma$ is a geodesic, for each $t_0 \in (a,b)$", font_size = 34).next_to(p1, DOWN)
        p3 = Tex(r"there is $\varepsilon > 0 $ such that $\gamma \vert _{[ t_0 - \varepsilon, t_0 + \varepsilon ]}$ is minimizing.", font_size = 34).next_to(p2, DOWN)
        pr0 = Tex(r"Proof", font_size = 34, color = BLUE).next_to(p3, DOWN)
        pr1 = Tex(r"For $t_0 \in (a,b)$, pick $\varepsilon > 0 $ such that for $t \in [0, 2 \varepsilon ]$,", font_size = 34).next_to(pr0, DOWN)
        pr2 = Tex(r"$\gamma (t_0 - \varepsilon + t)  = s(t , 0  )$", r" in ", r"polar coordinates", r" around ", r"$ p : = \gamma (t_0 - \varepsilon )$.", font_size = 34).next_to(pr1, DOWN)
        pr3 = Tex(r"For $r_0 < $ radius$(B)$, if ", r"$\beta$", r" $ : [0, \ell ] \to \Sigma $ connects ", r"$p$", r" to  ", r"$r = r_0$,", font_size = 34).next_to(pr2, DOWN)
        pr4 = Tex(r"length$(\beta) \geq r_0$. ", r"$\Rightarrow$ if $\beta$ connects $p$ to $\gamma (t_0 + \varepsilon ),$ length$(\beta) \geq 2 \varepsilon$.", font_size = 34).next_to(pr3, DOWN)
        pr5 = Tex(r"length$(\alpha) = r(\ell) \leq $ length$(\gamma) \leq $ length$(\alpha)$ ", r"$\Rightarrow$ $\theta ^{\prime} \equiv 0$, $\alpha \equiv \gamma$.", font_size = 34).next_to(pr4, DOWN)
        l = p0.get_left()[0]
        p1.shift((l-p1.get_left()[0])*RIGHT)
        p2.shift((l-p2.get_left()[0])*RIGHT)
        p3.shift((l-p3.get_left()[0])*RIGHT)
        pr0.shift((l-pr0.get_left()[0])*RIGHT)
        pr1.shift((l-pr1.get_left()[0])*RIGHT)
        pr2.shift((l-pr2.get_left()[0])*RIGHT)
        pr3.shift((l-pr3.get_left()[0])*RIGHT)
        pr4.shift((l-pr4.get_left()[0])*RIGHT)
        pr5.shift((l-pr5.get_left()[0])*RIGHT)
        gamli = Line([0,0,0], [3, 0, 0], color = BLUE).next_to(pr2[0], DOWN).shift(0.15*UP)
        pcli = Line([0,0,0], [2, 0, 0], color = ORANGE).next_to(pr2[2], DOWN).shift(0.15*UP)
        betli = Line([0,0,0], [0.5, 0, 0], color = PINK).next_to(pr3[1], DOWN).shift(0.15*UP)
        pli = Line([0,0,0], [0.5, 0, 0], color = GREEN).next_to(pr3[3], DOWN).shift(0.15*UP)
        pli0 = Line([0,0,0], [2, 0, 0], color = GREEN).next_to(pr2[4], DOWN).shift(0.15*UP)
        circli = Line([0,0,0], [1, 0, 0], color = RED).next_to(pr3[5], DOWN).shift(0.15*UP)
        dis = 2.7
        pb = Surface(
            lambda u, v: [np.sin(u) * np.sin(v) - dis/2, np.sin(u)* np.cos(v) + dis*sqrt(3) ,  np.cos(u) - 1.5],
            u_range = [0, PI/2],
            v_range = [0, TAU],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.1,
            resolution = 12
        )
        p = Sphere(radius = 0.07)
        p.set_color(GREEN).shift( [ - dis /2 , dis * sqrt(3) , -0.5 ] )
        bet = ParametricFunction(lambda t : [ np.sin(t) *np.cos( PI/4 +   np.sin( 3 * t ) /2 ) - dis/2 , np.sin(t) * np.sin( PI/4 +  np.sin( 3 * t ) /2  ) + dis * sqrt(3) , np.cos(t) - 1.5  ] , t_range = [0, PI/4], color = PINK)
        bet2 = ParametricFunction(lambda t : [ np.sin(t) *np.cos( PI/4 +   np.sin( 3 * t ) /2 ) - dis/2 , np.sin(t) * np.sin( PI/4 +  np.sin( 3 * t ) /2  ) + dis * sqrt(3) , np.cos(t) - 1.5  ] , t_range = [0, PI/3], color = PINK)
        circ = ParametricFunction(lambda t : [ np.sin(PI/4) *np.cos(t) - dis/2 , np.sin(PI/4) * np.sin( t ) + dis * sqrt(3) , np.cos(PI/4) - 1.5  ] , t_range = [0, 2*PI], color = RED)
        gam = ParametricFunction(lambda t : [ np.sin(t) / sqrt(2) - dis/2 , np.sin(t) / sqrt(2) + dis * sqrt(3) , np.cos(t) -1.5 ] , t_range = [0, PI/3], color = BLUE)
        self.add_fixed_in_frame_mobjects(p0, p1, p2, p3, pr0, pr1, pr2, pr3, pr4, pr5, gamli, pcli, betli, pli, circli, pli0)
        self.remove( pr1, pr2, pr3, pr4, pr5, gamli, pcli, betli, pli, circli, pli0)
        self.wait()
        self.play(Write(pr1), Write(pr2), Create(gamli), Create(pcli), Create(pli0), Create(gam), Create(p), Create(pb))
        self.wait()
        self.play(Write(pr3), Create(betli), Create(circli), Create(bet), Create(circ), Create(pr4[0]), Create(pli))
        self.wait()
        self.play(Create(pr4[1]), Create(bet2), FadeOut(circ))
        self.wait()
        
        #self.begin_ambient_camera_rotation(rate=PI/8)
        



class sgcb(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/6)
        b0 = Tex(r"For a semigeodesic chart $s : $ ",r"$ U$", r" $ \to $ ", r"$\Sigma$", r", set $b : U \to \mathbb{R}$ as $b :  = \vert s_v \vert $.", font_size = 34).shift( 2*UP )
        b1 = Tex(r"$s (u,v) = (u,v,0)$, $b \equiv 1$.", font_size = 34).next_to(b0, DOWN).shift(0.3*DOWN)
        uli = Line([0,0,0], [0.4, 0,0], color = GREEN).next_to(b0[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.4, 0,0], color = PURPLE_B).next_to(b0[3], DOWN).shift(0.15*UP)
        len = 1.8
        sar = CurvedArrow( [-1 ,-0.5 ,0],[1 ,-0.5,0], radius = -4 , tip_length = 0.2)
        sla = Tex(r"$s$", font_size = 34).next_to(sar, UP)
        s = Surface(
            lambda u, v: [u - len , v + len * sqrt(3), -1.5],
            u_range = [-2,2],
            v_range = [-2,2],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.4,
            resolution = 12
        )
        u = Surface(
            lambda u, v: [u + len  , v - len*sqrt(3), -1.5],
            u_range = [-2,2],
            v_range = [-2,2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.4,
            resolution = 12
        )
        self.add_fixed_in_frame_mobjects(b0, sar, b1, sla, uli, sli)
        self.remove(b0, sar, b1, sla, uli, sli)
        self.play(Write(b0), Create(uli), Create(sli))
        self.wait()
        self.play(Write(b1), Create(sar), Create(s), Create(u), Create(sla))
        self.wait()
        self.play(FadeOut(b1), FadeOut(sar), FadeOut(sla), FadeOut(u), FadeOut(s))
        self.wait()
        

class sgcbp(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/6)
        b0 = Tex(r"For a semigeodesic chart $s : $ ",r"$ U$", r" $ \to $ ", r"$\Sigma$", r", set $b : U \to \mathbb{R}$ as $b :  = \vert s_v \vert $.", font_size = 34).shift( 2*UP )
        b1 = Tex(r"$s (u,v) = (\sin(u)  \cos (v), \sin(u) \sin (v), \cos (u) )$,", font_size = 34).next_to(b0, DOWN).shift(0.5*UP)
        b2 = Tex(r"$b = \vert s_v \vert = \vert ( -  \sin (u)  \sin (v) , \sin (u) \cos (v) , 0  ) \vert = \sin (u)$.", font_size = 34).next_to(b1, DOWN)
        b3 = Tex(r"$\Rightarrow$ ", r"$ b_{uu} < 0$", font_size = 34).next_to(b2, DOWN)
        uli = Line([0,0,0], [0.4, 0,0], color = GREEN).next_to(b0[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.4, 0,0], color = PURPLE_B).next_to(b0[3], DOWN).shift(0.15*UP)
        len = 1.8
        sar = CurvedArrow( [-0.5 ,-1 ,0],[1.5 ,-1,0], radius = -4 , tip_length = 0.2)
        sla = Tex(r"$s$", font_size = 34).next_to(sar, UP)
        s = Surface(
            lambda u, v: [np.sin(u) * np.cos(v) - len , np.sin(u) * np.sin(v) + len *sqrt(3)  , np.cos(u)  -2],
            u_range = [0, PI ],
            v_range = [0,2*PI],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.4,
            resolution = 24
        )
        u = Surface(
            lambda u, v: [u + len , v - len * sqrt(3), -2],
            u_range = [-2,2],
            v_range = [-2,2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.4,
            resolution = 12
        )
        self.add_fixed_in_frame_mobjects(b0, sar, b1, sla, uli, sli, b2, b3)
        self.remove(sar, b1, sla, b2, b3)
        self.play(b0.animate.shift(0.5*UP), uli.animate.shift(0.5*UP), sli.animate.shift(0.5*UP) )
        self.wait()
        self.play(Write(b1), Create(sar), Create(s), Create(u), Create(sla))
        self.wait()
        self.play(Write(b2))
        self.wait()
        self.play(Write(b3))
        self.wait()
        self.play(FadeOut(b1), FadeOut(sar), FadeOut(sla), FadeOut(u), FadeOut(s), FadeOut(b2), FadeOut(b3))
        self.wait()
        

class sgcbn(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/6)
        b0 = Tex(r"For a semigeodesic chart $s : $ ",r"$ U$", r" $ \to $ ", r"$\Sigma$", r", set $b : U \to \mathbb{R}$ as $b :  = \vert s_v \vert $.", font_size = 34).shift( 2.5*UP )
        b01 = Tex(r"Take  ", r"$(\gamma_1(t), \gamma_2(t) )$", r" with $\gamma_1 > 0$, $\gamma_1^{\prime \prime } > 0$.", font_size = 34).next_to(b0, DOWN)
        b1 = Tex(r"$s (u,v) = ( \gamma_1 (u)  \cos (v), \gamma_2(u) ,  \gamma_1(u) \sin (v) )$,", font_size = 34).next_to(b01, DOWN)
        b2 = Tex(r"$b = \vert s_v \vert = \vert ( -  \gamma_1 (u)   \sin (v) , 0 , \gamma_1 (u)  \cos (v)   ) \vert = \gamma_1 (u) $.", font_size = 34).next_to(b1, DOWN)
        b3 = Tex(r"$\Rightarrow$ ", r"$ b_{uu} > 0$", font_size = 34).next_to(b2, DOWN)
        uli = Line([0,0,0], [0.4, 0,0], color = GREEN).next_to(b0[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.4, 0,0], color = PURPLE_B).next_to(b0[3], DOWN).shift(0.15*UP)
        gli = Line([0,0,0], [1.5,0,0], color = BLUE).next_to(b01[1], DOWN).shift(0.15*UP)
        len = 1.8
        bas = 2
        sar = CurvedArrow( [-1 ,-1 ,0],[1 ,-1,0], radius = -4 , tip_length = 0.2)
        sla = Tex(r"$s$", font_size = 34).next_to(sar, UP)
        s = Surface(
            lambda u, v: [ 5* bas**(-u) * np.cos(v) /16 - len ,  u  + len *sqrt(3)  ,  5*  bas**(-u)  * np.sin(v) /16 -2],
            u_range = [-2,2 ],
            v_range = [0,2*PI],
            checkerboard_colors = [PURPLE_B, PURPLE_B],
            fill_opacity = 0.4,
            resolution = 24
        )
        gam = ParametricFunction(lambda t : [ - len ,  t  + len *sqrt(3)  ,  5* bas**(-t) /16  -2 ] , t_range = [-2,2], color = BLUE)
        u = Surface(
            lambda u, v: [u + len , v - len * sqrt(3), -2],
            u_range = [-2,2],
            v_range = [-2,2],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.4,
            resolution = 12
        )
        self.add_fixed_in_frame_mobjects(b0, sar, b1, sla, uli, sli, b2, b3, b01, gli)
        self.remove(sar, b1, sla, b2, b3, b01, gli)
        self.play(Write(b01), Write(b1), Create(sar), Create(s), Create(u), Create(sla), Create(gam), Create(gli))
        self.wait()
        self.play(Write(b2))
        self.wait()
        self.play(Write(b3))
        self.wait()
        self.play(FadeOut(b1), FadeOut(sar), FadeOut(sla), FadeOut(u), FadeOut(s), FadeOut(b2), FadeOut(b3), FadeOut(gam), FadeOut(b01), FadeOut(gli))
        self.wait()
        

class sgcbb(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=75 * DEGREES, theta=  PI/6)
        b0 = Tex(r"For a semigeodesic chart $s : $ ",r"$ U$", r" $ \to $ ", r"$\Sigma$", r", set $b : U \to \mathbb{R}$ as $b :  = \vert s_v \vert $.", font_size = 34).shift( 2.5*UP )
        bp1 = Tex(r"$K > 0$", font_size = 34).next_to(b0, DOWN).shift(3.6*LEFT + DOWN)
        b01 = Tex(r"$K \equiv 0$", font_size = 34).next_to(b0, DOWN).shift(DOWN)
        bn1 = Tex(r"$K < 0$", font_size = 34).next_to(b0, DOWN).shift(3.6*RIGHT + DOWN)
        bp2 = Tex(r"$b_{uu} < 0$", font_size = 34).next_to(bp1, DOWN)
        b02 = Tex(r"$b_{uu} = 0$", font_size = 34).next_to(b01, DOWN)
        bn2 = Tex(r"$b_{uu} > 0$", font_size = 34).next_to(bn1, DOWN)
        p = VGroup(bp1, bp2)
        z = VGroup(b01, b02)
        n = VGroup(bn1, bn2)
        boxp = SurroundingRectangle(p, buff = .2, color = ORANGE)
        box0 = SurroundingRectangle(z, buff = .2, color = BLUE)
        boxn = SurroundingRectangle(n, buff = .2, color = YELLOW)
        uli = Line([0,0,0], [0.4, 0,0], color = GREEN).next_to(b0[1], DOWN).shift(0.15*UP)
        sli = Line([0,0,0], [0.4, 0,0], color = PURPLE_B).next_to(b0[3], DOWN).shift(0.15*UP)
        len = 1.8
        bas = 2
        sp = Surface(
            lambda u, v: [ np.sin(u)*np.cos(v) + len , np.sin(u)*np.sin(v) - len*sqrt(3), np.cos(u) - 1.5 ],
            u_range = [0 , PI ],
            v_range = [0,2*PI],
            checkerboard_colors = [ORANGE, ORANGE],
            fill_opacity = 0.4,
            resolution = 24
        )
        s0 = Surface(
            lambda u, v: [ u , v , -1.5 ],
            u_range = [-1 , 1 ],
            v_range = [-1 , 1 ],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.4,
            resolution = 24
        )
        sn = Surface(
            lambda u, v: [  bas**(-u) * np.cos(v) / 6 - 1.1*len , 3*u/4 + 1.1*len*sqrt(3) , bas**(-u)*np.sin(v) /6 - 1.5   ],
            u_range = [ -2, 2 ],
            v_range = [0,2*PI],
            checkerboard_colors = [YELLOW, YELLOW],
            fill_opacity = 0.4,
            resolution = 24
        )
        self.add_fixed_in_frame_mobjects(b0, bp1, bp2, b01, b02, bn1, bn2, uli, sli, boxp, box0, boxn)
        self.remove(bp1, bp2, b01, b02, bn1, bn2, boxp, box0, boxn)
        self.play(Write(bp1), Write(bp2), Create(sp), Write(b01), Write(b02), Create(s0), Write(bn1), Write(bn2), Create(sn), Create(boxp),Create(box0),Create(boxn) )
        self.wait()


class jac(Scene):
    def construct(self):
        d0 = Tex(r"Proposition (Jacobi equation)", font_size = 34, color = BLUE).to_edge(UL).shift(RIGHT + 0.8* DOWN)
        d1 = Tex(r"For a semigeodesic chart $s : U \to \Sigma$, if $b = \vert s _ v \vert$, then", font_size = 34).next_to(d0, DOWN)
        d2 = Tex(r"$bK + b_{uu} = 0.$", font_size = 34).next_to(d1, DOWN)
        e0 = Tex(r"Proof", font_size = 34, color = BLUE).next_to(d2, DOWN)
        e1 = Tex(r"Set $X = s_u / \vert s_u \vert = s_u$, $ $  $Y = s_v / \vert s_v \vert = s_v / b $. ", font_size = 34).next_to(e0, DOWN)
        e2 = Tex(r"Shape $ = \begin{pmatrix} \ell & m\\ m & n \end{pmatrix}	$", font_size = 34).next_to(e1, DOWN)
        e3 = Tex(r"$ X_u =$ ", r"$s_{uu}$ = ",  r" ?", font_size = 34).next_to(e2, DOWN)
        e3p = Tex(r"$ \ell N $", font_size = 34).next_to(e3[1], RIGHT).shift(1.25*LEFT + 0.05*UP)
        e4 = Tex(r"$ s_{uu} \cdot N = $ ", r"$  \partial_u (s_u \cdot N)  $ ", r"$- s_{u} \cdot N_u$ ", r"$ = s_u  \cdot  $ Shape$(s_u)=$ ", r"$ X \cdot $Shape$(X) = \ell$.", font_size = 34).next_to(e3, DOWN).shift(0.1*DOWN)
        l = d0.get_left()[0]
        d1.shift((l-d1.get_left()[0])*RIGHT)
        d2.shift((-d2.get_center()[0])*RIGHT)
        e0.shift((l-e0.get_left()[0])*RIGHT)
        e1.shift((l-e1.get_left()[0])*RIGHT)
        e2.shift((-e2.get_center()[0])*RIGHT)
        e3.shift((l-e3.get_left()[0])*RIGHT)
        e4.shift((l-e4.get_left()[0])*RIGHT)
        xu = VGroup(e3[0], e3p)
        can = Line([0,0,0],[1.5, 0.3, 0], color = RED).next_to(e4[1], DOWN).shift(0.6*UP)
        boxxu = SurroundingRectangle(xu, buff = .15, color = YELLOW)
        self.play(Write(d0), Write(d1), Write(d2))
        self.wait()
        self.play(Write(e0), Write(e1))
        self.wait()
        self.play(Write(e2))
        self.wait()
        self.play(Write(e3[0]))
        self.wait()
        self.play(Write(e3[1]), Write(e3[2]))
        self.wait()
        self.play(Write(e4[0]))
        self.wait()
        self.play(Write(e4[1]), Write(e4[2]))
        self.wait()
        self.play(Create(can))
        self.wait()
        self.play(Write(e4[3]))
        self.wait()
        self.play(Write(e4[4]))
        self.wait()
        self.play(FadeOut(e3[2]), Write(e3p), Create(boxxu))
        self.wait()
        #self.begin_ambient_camera_rotation(rate=PI/8)
        

class jac2(Scene):
    def construct(self):
        e1 = Tex(r"$Y_u \cdot X = $ ", r"$\partial_u ( Y \cdot X)$", r" $ - Y \cdot X_u  = $ ", r"$- \ell $ $ Y \cdot N  =$ " ,  r"$ 0$." , font_size = 34).to_edge(UL).shift(RIGHT + 0.5 * DOWN)
        e2 = Tex(r"$Y_u \cdot Y =$ ", r"$ \frac{1}{2} \partial _u ( Y \cdot Y ) $ ", r"$=0.$", font_size = 34).next_to(e1, DOWN)
        e3 = Tex(r"$Y_u \cdot N = $ ",r"$\partial_u(Y \cdot N)$", r" $- Y \cdot N_u=$", r" $Y \cdot $Shape$(X)=$ ", r"$m$." , font_size = 34).next_to(e2, DOWN)
        e4 = Tex(r"$ Y_u = m N$", font_size = 34).next_to(e3, DOWN).shift(0.1*DOWN)
        l = e1.get_left()[0]
        e2.shift((l-e2.get_left()[0])*RIGHT)
        e3.shift((l-e3.get_left()[0])*RIGHT)
        e4.shift((-e4.get_center()[0])*RIGHT)
        can = Line([0,0,0],[1.5, 0.3, 0], color = RED).next_to(e1[1], DOWN).shift(0.6*UP)
        can2 = Line([0,0,0],[1.5, 0.3, 0], color = RED).next_to(e3[1], DOWN).shift(0.6*UP)
        boxyu = SurroundingRectangle(e4, buff = .15, color = YELLOW)
        self.play(Write(e1[0]))
        self.wait()
        self.play(Write(e1[1]), Write(e1[2]))
        self.wait()
        self.play(Create(can))
        self.wait()
        self.play(Write(e1[3]))
        self.wait()
        self.play(Write(e1[4]))
        self.wait()
        self.play(Write(e2[0]))
        self.wait()
        self.play(Write(e2[1]))
        self.wait()
        self.play(Write(e2[2]))
        self.wait()
        self.play(Write(e3[0]))
        self.wait()
        self.play(Write(e3[1]), Write(e3[2]))
        self.wait()
        self.play(Create(can2))
        self.wait()
        self.play(Write(e3[3]))
        self.wait()
        self.play(Write(e3[4]))
        self.wait()
        self.play(Write(e4), Create(boxyu))
        self.wait()
        #self.begin_ambient_camera_rotation(rate=PI/8)
        

class jac3(Scene):
    def construct(self):
        e1 = Tex(r"$Y_u \cdot X = $ ", r"$\partial_u ( Y \cdot X)$", r" $ - Y \cdot X_u  = $ ", r"$- \ell $ $ Y \cdot N  =$ " ,  r"$ 0$." , font_size = 34).to_edge(UL).shift(RIGHT + 0.5 * DOWN)
        e2 = Tex(r"$Y_u \cdot Y =$ ", r"$ \frac{1}{2} \partial _u ( Y \cdot Y ) $ ", r"$=0.$", font_size = 34).next_to(e1, DOWN)
        e3 = Tex(r"$Y_u \cdot N = $ ",r"$\partial_u(Y \cdot N)$", r" $- Y \cdot N_u=$", r" $Y \cdot $Shape$(X)=$ ", r"$m$." , font_size = 34).next_to(e2, DOWN)
        e4 = Tex(r"$ Y_u = m N$", font_size = 34).next_to(e3, DOWN).shift(0.1*DOWN)
        e5 = Tex(r"$X_v = $ ", r"$ s_{uv} =$ ", r"$ (s_v)_u = $ ", r"$ (b Y)_u = $ ", r"$b_u Y + b Y_u $" , font_size = 34).next_to(e4, DOWN).shift(0.1*DOWN)
        e6 = Tex(r"$ X_v = b_u Y + bm N$", font_size = 34).next_to(e5, DOWN).shift(0.1*DOWN)
        l = e1.get_left()[0]
        e2.shift((l-e2.get_left()[0])*RIGHT)
        e3.shift((l-e3.get_left()[0])*RIGHT)
        e4.shift((-e4.get_center()[0])*RIGHT)
        e5.shift((l-e5.get_left()[0])*RIGHT)
        e6.shift((-e6.get_center()[0])*RIGHT)
        can = Line([0,0,0],[1.5, 0.3, 0], color = RED).next_to(e1[1], DOWN).shift(0.6*UP)
        can2 = Line([0,0,0],[1.5, 0.3, 0], color = RED).next_to(e3[1], DOWN).shift(0.6*UP)
        boxyu = SurroundingRectangle(e4, buff = .15, color = YELLOW)
        boxxv = SurroundingRectangle(e6, buff = .15, color = YELLOW)
        self.add(e1, e2, e3, e4, can, can2, boxyu)
        self.wait()
        self.play(Write(e5[0]))
        self.wait()
        self.play(Write(e5[1]))
        self.wait()
        self.play(Write(e5[2]))
        self.wait()
        self.play(Write(e5[3]))
        self.wait()
        self.play(Write(e5[4]))
        self.wait()
        self.play(Write(e6), Create(boxxv))
        self.wait()


class jac4(Scene):
    def construct(self):
        e1 = Tex(r"$Y_u \cdot X = $ ", r"$\partial_u ( Y \cdot X)$", r" $ - Y \cdot X_u  = $ ", r"$- \ell $ $ Y \cdot N  =$ " ,  r"$ 0$." , font_size = 34).to_edge(UL).shift(RIGHT + 0.5 * DOWN)
        e2 = Tex(r"$Y_u \cdot Y =$ ", r"$ \frac{1}{2} \partial _u ( Y \cdot Y ) $ ", r"$=0.$", font_size = 34).next_to(e1, DOWN)
        e3 = Tex(r"$Y_u \cdot N = $ ",r"$\partial_u(Y \cdot N)$", r" $- Y \cdot N_u=$", r" $Y \cdot $Shape$(X)=$ ", r"$m$." , font_size = 34).next_to(e2, DOWN)
        e4 = Tex(r"$ Y_u = m N$", font_size = 34).next_to(e3, DOWN).shift(0.1*DOWN)
        e5 = Tex(r"$X_v = $ ", r"$ s_{uv} =$ ", r"$ (s_v)_u = $ ", r"$ (b Y)_u = $ ", r"$b_u Y + b Y_u $" , font_size = 34).next_to(e4, DOWN).shift(0.1*DOWN)
        e6 = Tex(r"$ X_v = b_u Y + bm N$", font_size = 34).next_to(e5, DOWN).shift(0.1*DOWN)
        e7 = Tex(r"$Y_v =$ ", r"$ (N \times X)_v = $ ", r"$ N_v \times X + N \times X_v =$ " , r"$-$ Shape$(s_v) \times X  + N \times X_v =$ ", font_size = 34).next_to(e6, DOWN).shift(0.1*DOWN)
        e8 = Tex(r"$= - b$ Shape$(Y ) \times X  + N \times (b_u Y + $ ", r"$bm N$", r"$)=$ ",r" $ - b n Y \times X + b_u N \times Y  $" , font_size = 34).next_to(e7, DOWN)
        e9 = Tex(r"$Y_v = -b_u X +bnN $", font_size = 34).next_to(e8, DOWN).shift(0.1*DOWN)
        l = e1.get_left()[0]
        e2.shift((l-e2.get_left()[0])*RIGHT)
        e3.shift((l-e3.get_left()[0])*RIGHT)
        e4.shift((-e4.get_center()[0])*RIGHT)
        e5.shift((l-e5.get_left()[0])*RIGHT)
        e6.shift((-e6.get_center()[0])*RIGHT)
        e7.shift((l-e7.get_left()[0])*RIGHT)
        e8.shift((l-e8.get_left()[0])*RIGHT)
        e9.shift((-e9.get_center()[0])*RIGHT)
        can = Line([0,0,0],[1.5, 0.3, 0], color = RED).next_to(e1[1], DOWN).shift(0.6*UP)
        can2 = Line([0,0,0],[1.5, 0.3, 0], color = RED).next_to(e3[1], DOWN).shift(0.6*UP)
        can3 = Line([0,0,0],[0.8, 0.3, 0], color = RED).next_to(e8[1], DOWN).shift(0.5*UP)
        boxyu = SurroundingRectangle(e4, buff = .15, color = YELLOW)
        boxxv = SurroundingRectangle(e6, buff = .15, color = YELLOW)
        boxyv = SurroundingRectangle(e9, buff = .15, color = YELLOW)
        self.add(e1, e2, e3, e4, can, can2, boxyu, e5, e6, boxxv)
        self.wait()
        self.play(Write(e7[0]))
        self.wait()
        self.play(Write(e7[1]))
        self.wait()
        self.play(Write(e7[2]))
        self.wait()
        self.play(Write(e7[3]))
        self.wait()
        self.play(Write(e8[0]), Write(e8[1]), Write(e8[2]))
        self.wait()
        self.play(Create(can3))
        self.wait()
        self.play(Write(e8[3]))
        self.wait()
        self.play(Write(e9), Create(boxyv))
        self.wait()


class jac5(Scene):
    def construct(self):
        e1 = Tex(r"$X_u  = \ell N $", font_size = 34).to_edge(UL).shift(RIGHT +1.5* DOWN)
        e2 = Tex(r"$Y_u = m N$", font_size = 34).next_to(e1, DOWN).shift(0.1*DOWN)
        e3 = Tex(r"$ X_v = b_u Y + bm N$", font_size = 34).next_to(e1, RIGHT).shift(RIGHT)
        e4 = Tex(r"$Y_v = -b_u X +bnN $", font_size = 34).next_to(e2, RIGHT).shift(RIGHT)
        bk = Tex(r"$bK = b\ell n - b m^2=$ ", r"$ X_u \cdot Y_v - X_v \cdot Y_u$", font_size = 34).shift(1.4*LEFT )
        bk2 = Tex(r"$= \left[ \partial_u(X  \cdot Y_v) - X \cdot Y_{uv } \right]  - \left[ \partial_v (X \cdot Y_u) - X \cdot Y_{uv } \right]$", font_size = 34).next_to(bk, DOWN)
        bk3 = Tex(r"$= - \partial_u( b_u ) $ ", font_size = 34).next_to(bk2, DOWN)
        bk4 = Tex(r"$ bK + b_{uu} = 0 $", font_size = 34).next_to(bk3, DOWN)
        l = bk.get_left()[0]
        bk2.shift((l - bk2.get_left()[0]+ 0.58)*RIGHT)
        bk3.shift((l - bk3.get_left()[0] + 0.58)*RIGHT)
        bk4.shift(( - bk4.get_center()[0])*RIGHT)
        eq = VGroup(e1, e2, e3, e4)
        eq.shift(eq.get_center()[0]*LEFT)
        box = SurroundingRectangle(eq, buff = .15, color = YELLOW)
        boxj = SurroundingRectangle(bk4, buff = .15, color = BLUE)
        can = Line([0,0,0],[1, 0.3, 0], color = RED).next_to(bk2, DOWN).shift(0.6*UP+0.8*LEFT)
        can2 = Line([0,0,0],[1, 0.3, 0], color = RED).next_to(bk2, DOWN).shift(0.6*UP+3*RIGHT)
        can3 = Line([0,0,0],[1.5, 0.3, 0], color = RED).next_to(bk2, DOWN).shift(0.6*UP + 1.2* RIGHT)
        self.play(Write(e1), Write(e2), Write(e3), Write(e4), Create(box))
        self.wait()
        self.play(Write(bk[0]))
        self.wait()
        self.play(Write(bk[1]))
        self.wait()
        self.play(Write(bk2))
        self.wait()
        self.play(Create(can),Create(can2))
        self.wait()
        self.play(Create(can3))
        self.wait()
        self.play(Write(bk3))
        self.wait()
        self.play(Write(bk4), Create(boxj))
        self.wait()


class egregium(Scene):
    def construct(self):
        j0 = Tex(r"Proposition (Jacobi equation)", font_size = 34, color = BLUE).to_edge(UL).shift(RIGHT + 0.8* DOWN)
        j1 = Tex(r"For a semigeodesic chart $s : U \to \Sigma$, if $b = \vert s _ v \vert$, then", font_size = 34).next_to(j0, DOWN)
        j2 = Tex(r"$bK + b_{uu} = 0.$", font_size = 34).next_to(j1, DOWN)
        t0 = Tex(r"Theorem (Egregium)", font_size = 34, color = PINK).next_to(j2, DOWN).shift(0.2*DOWN)
        t1 = Tex(r"Let  $\phi : \Sigma _ 1 \to \Sigma_2 $ be such that for all $p \in \Sigma_1, $", font_size = 34).next_to(t0, DOWN)
        t2 = Tex(r"$ \vert d_p \phi (X) \vert = \vert X \vert $ for all $X \in T_p \Sigma _1 $,", font_size = 34).next_to(t1, DOWN)
        t3 = Tex(r"Then", font_size = 34).next_to(t2, DOWN)
        t4 = Tex(r"$K_{\Sigma_2} (\phi (p))  = K_{\Sigma _1} ( p)$ for all $p \in \Sigma _1$.", font_size = 34).next_to(t3, DOWN).shift(0.1*DOWN)
        l = j0.get_left()[0]
        j1.shift((l-j1.get_left()[0])*RIGHT)
        j2.shift((-j2.get_center()[0])*RIGHT)
        t0.shift((l-t0.get_left()[0])*RIGHT)
        t1.shift((l-t1.get_left()[0])*RIGHT)
        t2.shift((-t2.get_center()[0])*RIGHT)
        t3.shift((l-t3.get_left()[0])*RIGHT)
        t4.shift((-t4.get_center()[0])*RIGHT)
        boxj = SurroundingRectangle(j2, buff = .15, color = BLUE)
        boxe = SurroundingRectangle(t4, buff = .15, color = PINK)
        self.play(Create(j0), Create(j1), Create(j2), Create(boxj))
        self.wait()
        self.play(Write(t0), Write(t1), Write(t2), Write(t3), Write(t4), Create(boxe))
        self.wait()





class egex(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=77 * DEGREES, theta=  PI/6)
        t0 = Tex(r"Theorem (Egregium)", font_size = 34, color = PINK).to_edge(UL).shift(RIGHT + 0.2* DOWN)
        t1 = Tex(r"Let  $\phi : $ ", r"$\Sigma _ 1$", r" $ \to $ ", r"$\Sigma_2 $", r" be such that for all $p \in \Sigma_1, $", font_size = 34).next_to(t0, DOWN)
        t2 = Tex(r"$ \vert d_p \phi (X) \vert = \vert X \vert $ for all $X \in T_p \Sigma _1 $,", font_size = 34).next_to(t1, DOWN)
        t3 = Tex(r"Then", font_size = 34).next_to(t2, DOWN)
        t4 = Tex(r"$K_{\Sigma_2} (\phi (p))  = K_{\Sigma _1} ( p)$ for all $p \in \Sigma _1$.", font_size = 34).next_to(t3, DOWN).shift(0.1*DOWN)
        left = t0.get_left()[0]
        t1.shift((left-t1.get_left()[0])*RIGHT)
        t2.shift((-t2.get_center()[0])*RIGHT)
        t3.shift((left-t3.get_left()[0])*RIGHT)
        t4.shift((-t4.get_center()[0])*RIGHT)
        s1li = Line([0,0,0], [0.6, 0,0], color = BLUE).next_to(t1[1], DOWN).shift(0.15*UP)
        s2li = Line([0,0,0], [0.6, 0,0], color = GREEN).next_to(t1[3], DOWN).shift(0.15*UP)
        boxe = SurroundingRectangle(t4, buff = .15, color = PINK)
        self.add_fixed_in_frame_mobjects(t0, t1, t2, t3, t4, boxe, s1li, s2li)
        r = 1
        s1 = Surface(
            lambda u, v: [r*u,r*v,-1.8 ],
            u_range = [-2*PI/3, 2*PI/3],
            v_range = [-2*PI/3, 2*PI/3],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.4,
            resolution = 24
        )
        s2 = Surface(
            lambda u, v: [ r*u , r*np.sin(v) ,r* np.cos(v) -1.8 ],
            u_range = [-2*PI/3, 2*PI/3 ],
            v_range = [-2*PI/3, 2*PI/3],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.4,
            resolution = 24
        )
        self.wait()
        self.play(Create(s1))
        self.wait()
        self.play(Transform(s1, s2))
        self.wait()



class egex2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi=72 * DEGREES, theta=  PI/6)
        t0 = Tex(r"Theorem (Egregium)", font_size = 34, color = PINK).to_edge(UL).shift(RIGHT + 0.2* DOWN)
        t1 = Tex(r"Let  $\phi : $ ", r"$\Sigma _ 1$", r" $ \to $ ", r"$\Sigma_2 $", r" be such that for all $p \in \Sigma_1, $", font_size = 34).next_to(t0, DOWN)
        t2 = Tex(r"$ \vert d_p \phi (X) \vert = \vert X \vert $ for all $X \in T_p \Sigma _1 $,", font_size = 34).next_to(t1, DOWN)
        t3 = Tex(r"Then", font_size = 34).next_to(t2, DOWN)
        t4 = Tex(r"$K_{\Sigma_2} (\phi (p))  = K_{\Sigma _1} ( p)$ for all $p \in \Sigma _1$.", font_size = 34).next_to(t3, DOWN).shift(0.1*DOWN)
        left = t0.get_left()[0]
        t1.shift((left-t1.get_left()[0])*RIGHT)
        t2.shift((-t2.get_center()[0])*RIGHT)
        t3.shift((left-t3.get_left()[0])*RIGHT)
        t4.shift((-t4.get_center()[0])*RIGHT)
        s1li = Line([0,0,0], [0.6, 0,0], color = BLUE).next_to(t1[1], DOWN).shift(0.15*UP)
        s2li = Line([0,0,0], [0.6, 0,0], color = GREEN).next_to(t1[3], DOWN).shift(0.15*UP)
        boxe = SurroundingRectangle(t4, buff = .15, color = PINK)
        self.add_fixed_in_frame_mobjects(t0, t1, t2, t3, t4, boxe, s1li, s2li)
        r = 3.1
        l = 1.5
        s1 = Surface(
            lambda u, v: [r* np.cos(u) * np.sin(v) ,r*  np.sin(u), r* np.cos(u)*np.cos(v) - 3.8 ],
            u_range = [-0.6,0.6 ],
            v_range = [-0.6,0.6],
            checkerboard_colors = [BLUE, BLUE],
            fill_opacity = 0.4,
            resolution = 24
        )
        s2 = Surface(
            lambda u, v: [ r* l* np.cos(u) * np.sin(v/l) , r*  u - r* (l**2)*u**3/6, r* l* np.cos(u)*np.cos(v/l) - r*(l-1) - 3.8 ],
            u_range = [-0.6,0.6 ],
            v_range = [-0.6,0.6],
            checkerboard_colors = [GREEN, GREEN],
            fill_opacity = 0.4,
            resolution = 24
        )
        self.wait()
        self.play(Create(s1))
        self.wait()
        self.play(Transform(s1, s2))
        self.wait()



class egco(Scene):
    def construct(self):
        t0 = Tex(r"Theorem (Egregium)", font_size = 34, color = PINK).to_edge(UL).shift(RIGHT + 0.2* DOWN)
        t1 = Tex(r"Let  $\phi : $ ", r"$\Sigma _ 1$", r" $ \to $ ", r"$\Sigma_2 $", r" be such that for all $p \in \Sigma_1, $", font_size = 34).next_to(t0, DOWN)
        t2 = Tex(r"$ \vert d_p \phi (X) \vert = \vert X \vert $ for all $X \in T_p \Sigma _1 $,", font_size = 34).next_to(t1, DOWN)
        t3 = Tex(r"Then", font_size = 34).next_to(t2, DOWN)
        t4 = Tex(r"$K_{\Sigma_2} (\phi (p))  = K_{\Sigma _1} ( p)$ for all $p \in \Sigma _1$.", font_size = 34).next_to(t3, DOWN).shift(0.1*DOWN)
        c0 = Tex(r"Corollary", font_size = 34, color = ORANGE).next_to(t4, DOWN).shift(0.3*DOWN)
        c1 = Tex(r"There is no flat map that represents a portion of the earth", font_size = 34).next_to(c0, DOWN)
        c2 = Tex(r"without distorting the distances.", font_size = 34).next_to(c1, DOWN)
        left = t0.get_left()[0]
        t1.shift((left-t1.get_left()[0])*RIGHT)
        t2.shift((-t2.get_center()[0])*RIGHT)
        t3.shift((left-t3.get_left()[0])*RIGHT)
        t4.shift((-t4.get_center()[0])*RIGHT)
        c0.shift((left-c0.get_left()[0])*RIGHT)
        c1.shift((left-c1.get_left()[0])*RIGHT)
        c2.shift((left-c2.get_left()[0])*RIGHT)
        s1li = Line([0,0,0], [0.6, 0,0], color = BLUE).next_to(t1[1], DOWN).shift(0.15*UP)
        s2li = Line([0,0,0], [0.6, 0,0], color = GREEN).next_to(t1[3], DOWN).shift(0.15*UP)
        boxe = SurroundingRectangle(t4, buff = .15, color = PINK)
        self.add(t0, t1, t2, t3, t4, boxe, s1li, s2li)
        self.wait()
        self.play(Create(c0), Create(c1), Create(c2))
        self.wait()

