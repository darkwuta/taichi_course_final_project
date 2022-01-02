import math

import taichi as ti

#import Exporter as ex
import numpy as np
ti.init(arch=ti.gpu)
ti.init(device_memory_fraction=0.9)



@ti.data_oriented
class PBFFluidSim3D:
    def __init__(self,
                 num_particles = 10240,
                 max_timesteps = 10,
                 time_delta=1.0 / 20.0,
                 BaseSpaceReslution = 400,
                 #boundary = (400, 400, 400),
                 gui = None,
                 do_render = False,
                 render_save_dir = "./output"
                 ):


        self.isPause = False
        self.base_size_factor = 400
        self.scaling_size_factor = 1
        self.dim = 3
        self.particles_color = 0x068587
        self.boundary_color = 0xebaca2
        self.background_color = 0x112f41

        self.gui = gui
        self.do_render = do_render
        self.render_save_dir = render_save_dir

        self.SpaceResolution = np.array([1, 1, 1]) * self.base_size_factor * self.scaling_size_factor
        self.screen_to_world_ratio = 10.0 * self.scaling_size_factor
        self.boundary = self.SpaceResolution / self.screen_to_world_ratio
        self.SpaceResolution_reciprocal = 1/self.SpaceResolution[0], 1/self.SpaceResolution[1], 1/self.SpaceResolution[2]

        self.cell_size = 2.51 / self.scaling_size_factor
        self.cell_reciprocal = 1.0 / self.cell_size

        # 向上取整
        def round_up(f, s):
            return (math.floor(f * self.cell_reciprocal / s)+1)*s

        print(self.boundary)

        self.grid_size = (round_up(self.boundary[0], 1), round_up(self.boundary[1], 1), round_up(self.boundary[2], 1))
        ####
        self.max_timesteps = max_timesteps
        self.num_particles = num_particles
        self.max_num_particles_per_cell = 200
        self.max_num_neighbors = 200
        self.time_delta = time_delta
        self.epsilon = 1e-5

        # TODO 这里原来是
        # self.particle_radius_in_world = 0.3
        # self.particle_radius = self.particle_radius_in_world * render_scaling
        # 可能会有问题
        self.particles_radius = 3.0
        self.particles_radius_in_world = self.particles_radius / self.screen_to_world_ratio


        # PBF parameter
        self.h = 1.1
        self.mass = 1.0
        self.rho0 = 1.0
        self.lambda_epsilon = 100
        self.vorticity_epsilon = 0.01
        self.viscosity_c = 0.01
        self.pbd_num_iters = 5
        self.corr_deltaQ_coeff = 0.3
        self.corrK = 0.001

        self.neighbor_radius = self.h * 1.05

        self.poly6_factor = 315.0 / 64.0 / np.pi
        self.spiky_grad_factor = -45.0 / np.pi

        #self.target = ti.Vector.field(self.dim, dtype = ti.f32)

        #self.total_pos_delta = ti.Vector.field(self.dim, dtype = ti.f32)
        self.positions = ti.Vector.field(self.dim, dtype = ti.f32)
        self.old_positions = ti.Vector.field(self.dim, dtype = ti.f32)
        self.velocities = ti.Vector.field(self.dim, dtype = ti.f32)

        

        self.grid_num_particles = ti.field(dtype = ti.i32)
        self.grid2particles = ti.field(dtype = ti.i32)
        self.particles_num_neighbors = ti.field(ti.f32)
        self.particles_neighbors = ti.field(ti.i32)
        self.lambdas = ti.field(ti.f32)
        self.lambdas_grad_i = ti.Vector.field(self.dim, ti.f32)
        self.lambdas_sum_gradient_sqr_i = ti.field(ti.f32)
        self.lambdas_density_constraints_i = ti.field(ti.f32)

        self.position_deltas = ti.Vector.field(self.dim, dtype = ti.f32)

        self.board_states = ti.Vector.field(2, dtype = ti.f32)


        self.place_vars()

        print(f'boundary={self.boundary} grid{self.grid_size} cell_size={self.cell_size}')

    # 初始化碰撞检测类的顶点信息和顶点坐标信息
    def init_surface(self):
        # 将numpy数组转为field
        self.fluid_surface_solver.init_field()
        # 初始化之后的赋值
        self.fluid_surface_solver.numpy_to_field()

    def place_vars(self):
        ti.root.dense(ti.i, self.num_particles).place(self.old_positions, self.positions, self.velocities)
        ti.root.dense(ti.i, self.num_particles).place(self.lambdas, self.position_deltas)

        grid_snode = ti.root.dense(ti.ijk, self.grid_size)
        grid_snode.place(self.grid_num_particles)
        grid_snode.dense(ti.l, self.max_num_particles_per_cell).place(self.grid2particles)

        nb_node = ti.root.dense(ti.i, self.num_particles)
        nb_node.place(self.particles_num_neighbors)
        nb_node.dense(ti.j, self.max_num_neighbors).place(self.particles_neighbors)#particles_neighbors存储每个粒子的200个邻居的索引

        ti.root.dense(ti.i, self.num_particles).place(self.lambdas_grad_i)
        ti.root.dense(ti.i, self.num_particles).place(self.lambdas_sum_gradient_sqr_i)
        ti.root.dense(ti.i, self.num_particles).place(self.lambdas_density_constraints_i)
        ti.root.place(self.board_states)

        ti.root.lazy_grad()


    def init(self):
        self.init_particles()
        self.init_board()

    @ti.kernel
    def init_particles(self):
        ratio = self.SpaceResolution[0] / self.screen_to_world_ratio
        for i in range(self.num_particles):
            self.old_positions[i] = ti.Vector([(ti.random() * 0.4 + 0.3) * ratio,
                                               (ti.random() * 0.4 + 0.3) * ratio,
                                               (ti.random() * 0.4 + 0.3) * ratio])
            self.positions[i] = ti.Vector([(ti.random() * 0.4 + 0.3) * ratio,
                                               (ti.random() * 0.4 + 0.3) * ratio,
                                               (ti.random() * 0.4 + 0.3) * ratio])
            self.velocities[i] = ti.Vector([0.0, -1.0, 0.0])
            print(self.positions[i])
    @ti.kernel
    def init_board(self):
        self.board_states[None] = ti.Vector([self.boundary[0] - self.epsilon, -0.0])

    @ti.kernel
    def swap_buffer(self, frame: ti.i32):
        for i in range(self.num_particles):
            self.old_positions[i] = self.positions[i]

    @ti.kernel
    def apply_gravity_within_boundary(self, frame:ti.i32):
        for i in self.positions:
            g = ti.Vector([0.0, -9.8, 0.0])
            pos, vel = self.positions[i], self.velocities[i]
            vel += g * self.time_delta
            pos += vel * self.time_delta
            self.positions[i] = self.confine_position_to_boundary(pos)
            # print("pos:",pos,"position:",self.positions[i])

    @ti.func
    def get_cell(self, pos):
        return (pos * self.cell_reciprocal).cast(int)

    @ti.func
    def is_in_grid(self, c):
        return 0 <= c[0] and c[0] < self.grid_size[0] and 0 <= c[1] and c[1
        ] < self.grid_size[1] and 0<= c[2] and c[2] < self.grid_size[2]

    @ti.func
    def poly6_value(self, s, h):
        result = 0.0
        if 0 < s and s < h:
            x = (h * h - s * s) / (h * h * h)
            result = self.poly6_factor * x * x * x
        return result

    @ti.func
    def spiky_gradient(self, r, h):
        result = ti.Vector([0.0, 0.0, 0.0])
        r_len = r.norm()
        if 0 < r_len and r_len < h:
            x = (h - r_len) / (h * h * h)
            g_factor = self.spiky_grad_factor * x * x
            result = r * g_factor / r_len
        return result

    @ti.func
    def compute_scorr(self, pos_ji):
        # Eq (13)
        x = self.poly6_value(pos_ji.norm(), self.h) / self.poly6_value(self.corr_deltaQ_coeff * self.h, self.h)

        x = x * x
        x = x * x
        return (-self.corrK) * x

    @ti.func
    def confine_position_to_boundary(self, p):
        bmin = self.particles_radius_in_world
        # First coordinate is for the x position of the board, which only moves in x
        bmax = ti.Vector([self.board_states[None][0], self.boundary[1], self.boundary[2]]) - self.particles_radius_in_world

        for i in ti.static(range(self.dim)):
            # Use randomness to prevent particles from sticking into each other after clamping
            if p[i] <= bmin:
                p[i] = bmin + self.epsilon * ti.random()
            elif bmax[i] <= p[i]:
                p[i] = bmax[i] - self.epsilon * ti.random()
        return p




    @ti.kernel
    def update_grid(self, frame:ti.i32):
        # TODO 有疑惑
        for p_i in self.positions:
            cell = self.get_cell(self.positions[p_i])
            offs = self.grid_num_particles[cell].atomic_add(1)
            self.grid2particles[cell, offs] = p_i

    @ti.kernel
    def find_particle_neighbors(self, frame:ti.i32):
        for p_i in self.positions:
            pos_i = self.positions[p_i]
            cell = self.get_cell(pos_i)
            nb_i = 0
            for offs in ti.static(ti.grouped(ti.ndrange((-1, 2), (-1, 2), (-1, 2)))):
                cell_to_check = cell + offs
                if self.is_in_grid(cell_to_check):
                    for j in range(self.grid_num_particles[cell_to_check]):
                        p_j = self.grid2particles[cell_to_check, j]
                        if nb_i < self.max_num_neighbors and p_j !=p_i and (pos_i - self.positions[p_j]).norm() < self.neighbor_radius:
                            self.particles_neighbors[p_i, nb_i] = p_j
                            nb_i += 1
            self.particles_num_neighbors[p_i] = nb_i


    @ti.kernel
    def compute_lambdas(self,frame:ti.i32):
        # Eq (8) ~ (11)
        for p_i in self.positions:
            pos_i = self.positions[p_i]

            grad_i = ti.Vector([0.0, 0.0, 0.0])
            sum_gradient_sqr = 0.0
            density_constraint = 0.0

            for j in range(self.particles_num_neighbors[p_i]):
                p_j = self.particles_neighbors[p_i, j]
                if p_j >= 0:
                    pos_ji = pos_i - self.positions[p_j]
                    grad_j = self.spiky_gradient(pos_ji, self.h)
                    grad_i += grad_j
                    sum_gradient_sqr += grad_j.dot(grad_j)
                    density_constraint += self.poly6_value(pos_ji.norm(), self.h)
            density_constraint = (self.mass * density_constraint / self.rho0) - 1.0

            sum_gradient_sqr += grad_i.dot(grad_i)
            self.lambdas[p_i] = (-density_constraint)/(sum_gradient_sqr+self.lambda_epsilon)


    @ti.kernel
    def compute_position_deltas(self, frame: ti.i32):
        # Eq(12), (14)
        for p_i in self.positions:
            pos_i = self.positions[p_i]
            lambda_i = self.lambdas[p_i]

            pos_delta_i = ti.Vector([0.0, 0.0, 0.0])
            for j in range(self.particles_num_neighbors[p_i]):
                p_j = self.particles_neighbors[p_i, j]
                if p_j >=0:
                    lambda_j = self.lambdas[p_j]
                    pos_ji = pos_i - self.positions[p_j]
                    scorr_ij = self.compute_scorr(pos_ji)
                    pos_delta_i += (lambda_i + lambda_j + scorr_ij) * self.spiky_gradient(pos_ji, self.h)
            pos_delta_i /= self.rho0
            self.position_deltas[p_i] = pos_delta_i

    @ti.kernel
    def confine_to_boundary(self):
        for i in self.positions:
            pos = self.positions[i]
            self.positions[i] = self.confine_position_to_boundary(pos)

    @ti.kernel
    def apply_position_deltas(self,frame:ti.i32):
        for i in self.positions:
            self.positions[i] += self.position_deltas[i]

    @ti.kernel
    def update_velocities(self, frame: ti.i32):
        for i in self.positions:
            self.velocities[i] = (self.positions[i]-self.old_positions[i])/self.time_delta

    @ti.kernel
    def vorticity_confinement(self,frame:ti.i32):
        # TODO 这里还没有明白为什么
        for p_i in self.positions:
            pos_i = self.positions[p_i]
            vel_i = self.velocities[p_i]
            omega_i = ti.Vector([0.0, 0.0, 0.0])

            omega_mag_i_grad = ti.Vector([0.0, 0.0, 0.0])

            for j in range(self.particles_num_neighbors[p_i]):
                p_j = self.particles_neighbors[p_i, j]
                if p_j >= 0:
                    pos_ji = pos_i - self.positions[p_j]
                    vel_ij = self.velocities[p_j] - vel_i
                    grad_j = -self.spiky_gradient(pos_ji, self.h)
                    z = vel_ij.cross(grad_j)
                    z_grad =ti.Matrix([[0.0,vel_ij[2], -vel_ij[1]],[-vel_ij[2],0.0,vel_ij[0]],[vel_ij[1], -vel_ij[0], 0.0]])
                    omega_mag_i_grad += z_grad @ z # @ 是矩阵乘法运算
                    omega_i +=vel_ij.cross(z)

            omega_mag_i = omega_i.norm()

            omega_mag_i_grad /= omega_mag_i

            location_vector = omega_mag_i_grad.normalized()
            f_vorticity = self.vorticity_epsilon * (location_vector.cross(omega_i))
            self.velocities[p_i] += self.time_delta * f_vorticity


    @ti.kernel
    def apply_XSPH_viscosity(self,frame:ti.i32):
        # TODO 这里还没有明白为什么
        # Eq (17)
        for v_i in self.velocities:
            pos_i = self.positions[v_i]
            vel_i = self.velocities[v_i]
            v_delta_i = ti.Vector([0.0, 0.0, 0.0])

            for j in range(self.particles_num_neighbors[v_i]):
                p_j = self.particles_neighbors[v_i, j]
                pos_ji = pos_i - self.positions[p_j]
                vel_ij = self.velocities[p_j] - vel_i
                if p_j >= 0:
                    pos_ji = pos_i - self.positions[p_j]
                    v_delta_i += vel_ij * self.poly6_value(pos_ji.norm(), self.h)

            self.velocities[v_i] += self.viscosity_c * v_delta_i

    @ti.kernel
    def move_board(self, frame:ti.i32):
        b = self.board_states[None]
        period = 90
        vel_strength = 8.0
        b[1] += 1.0
        if b[1] >= 2*period:
            b[1] = 0
        b[0] += -ti.sin(b[1] * np.pi / period) * vel_strength * self.time_delta
        self.board_states[None] = b
        #print("DEBUG::self.board_states[None]", self.board_states[None])


    def gui_show(self,frame):
        self.gui.clear(self.background_color)
        ratio = self.screen_to_world_ratio*self.SpaceResolution_reciprocal[0]
        # boundary , x-axis
        start = np.array([self.board_states[None][0]*ratio, 0.0])
        end = np.array([self.board_states[None][0]*ratio, 1.0])
        self.gui.line(start, end, radius=5, color=0xFFFFFF)
        # TODO 画出粒子位置

        pos_np = self.positions.to_numpy()
        for pos in pos_np:
            for j in range(self.dim):
                    pos[j] *= ratio

        self.gui.circles(pos_np[:, [0, 1]], radius=3, color=self.particles_color)
        self.gui.show()

    def run_pbf(self, frame):
        if not self.isPause:
            self.swap_buffer(frame)
            self.apply_gravity_within_boundary(frame)

            # clear grid
            self.grid_num_particles.fill(0)
            self.particles_neighbors.fill(-1)

            self.update_grid(frame)
            self.find_particle_neighbors(frame)

            for _ in range(self.pbd_num_iters):
                self.compute_lambdas(frame)
                self.compute_position_deltas(frame)
                self.apply_position_deltas(frame)

            self.confine_to_boundary()

            self.update_velocities(frame)

            #self.vorticity_confinement(frame)

            self.apply_XSPH_viscosity(frame)

