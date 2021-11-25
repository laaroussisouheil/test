from manim import *
from numpy.lib.function_base import select
class pith (Scene):
    def contructur(self):
        sq=Square(side_length=2,stroke_color=GREEN,fill_color=BLUE,fill_opacity=0.75)
        self.play(Create(sq),run_time=3)
        self.wait()