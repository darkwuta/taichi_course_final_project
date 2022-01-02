from Fluid import *
from Gui.ggui_Scene.ggui_Scene import *

def main():
    screen_resolution = 1280, 720
    #max_frame = 50000

    pbfFluid = PBFFluidSim3D()
    scene = gguiScene(res =screen_resolution, objects = pbfFluid)
    pbfFluid.init()


    frame = 0
    while scene.window.running:
        frame += 1
        pbfFluid.run_pbf(frame)
        pbfFluid.move_board(frame)
        scene.render()


if __name__ == "__main__":
    main()