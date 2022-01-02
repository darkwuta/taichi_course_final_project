# 太极图形课S1-大作业-PBF

## 作业来源
​		该大作业本来是Games101的大作业（不知道为什么smart chain没有交上去QAQ）复现了PBF的算法，当时没有ggui，是先在taichi里进行模拟，然后保存顶点坐标，然后使用marching cubes方法对粒子表面重建，最后把obj序列导入blender中进行渲染完成的

本来作业地址在这里[here](https://github.com/darkwuta/GAMES101-Assignment/blob/main/FinalProject/FinalProject.md)：

本次大作业在写好的PBF的基础上增加ggui

因为期末周实在QAQ没有空就只能用这个交大作业了( 本来想加个marching cube)，而且太久没有看很多关于PBF的细节已经忘掉了，代码很乱，推荐可以看一下代码的主要参考[here](https://github.com/ben441318936/PBF3D_taichi)

 ps：其实就是把example的PBF改成3D加了ggui，没时间了，对不起，我忏悔，我有罪

参考论文：

[1] Position Based Fluids

## 运行方式
#### 运行环境：
> [Taichi] version 0.8.5, llvm 10.0.0, commit 45c6ad48, win, python 3.8.1

#### 运行：
> 运行PositionBasedFluid/PDF3D.py

## 效果展示
> 

## 整体结构
```
-LICENSE
-PositionBasedFluid(主要代码)
 -Fluid.py
 -PBF3D.py
-Gui
 -ggui_Scene
  -ggui_Scene.py
-References(参考论文)
 -pbf_sig_preprint.pdf
 -Position Based Fluids.md(个人翻译)
 -Position Based Fluids.pdf
-README.MD
```

## 实现细节：

`Fluid.py`是项目源代码

### 整体流程

算法伪代码：

![image-20220102222143388](image/image-20220102222143388.png)

**step1(1-4)：遍历所有粒子，更新速度和位置**

**step2(4,7)：遍历所有粒子，找到相邻粒子**

​	寻找相邻粒子要提前构建粒子的相邻粒子表

```python
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
```

**step3(8,19)：迭代求出位置偏移量并将其应用于当前位置**

**step4(20,24)：遍历所有粒子，利用偏移后的位置求出速度，然后更新位置**

公式的计算可以看References文件夹中的论文以及个人翻译

