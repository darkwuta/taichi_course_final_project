import numpy as np

import taichi as ti

class gguiScene:
    paused = False
    particles_radius = 0.2
    material_colors = [(0.1, 0.6, 0.9), (0.93, 0.33, 0.23), (1.0, 1.0, 1.0)]
    res = (1280, 720)
    curr_preset_id = 0
    colors_used = ([0.1, 0.6, 0.9])

    def __init__(self, res, objects):
        self.res = res
        self.frame_id = 0
        self.window = ti.ui.Window("PBF", res)
        self.canvas = self.window.get_canvas()
        self.scene = ti.ui.Scene()

        self.camera = ti.ui.make_camera()
        self.camera.position(30, 30, 30)
        self.camera.lookat(0.5, 0.3, 0.5)
        self.camera.fov(55)
        self.objects = objects

    def show_options(self):
        self.window.GUI.begin("Options", 0.05, 0.45, 0.2, 0.4)
        if self.window.GUI.button("restart"):
            self.objects.init()
        if self.paused:
            if self.window.GUI.button("Continue"):
                paused = False
        else:
            if self.window.GUI.button("Pause"):
                self.paused = True
        self.window.GUI.end()

    def render(self):
        self.camera.track_user_inputs(self.window, movement_speed=0.03,hold_key=ti.ui.RMB)
        self.scene.set_camera(self.camera)

        self.scene.ambient_light((1, 1, 1))

        self.scene.particles(self.objects.positions,  radius = self.particles_radius)

        self.scene.point_light(pos=(0.5, 1.5, 0.5), color=(0.5, 0.5, 0.5))
        self.scene.point_light(pos=(0.5, 1.5, 1.5), color=(0.5, 0.5, 0.5))

        self.canvas.scene(self.scene)

        self.show_options()

        self.window.show()