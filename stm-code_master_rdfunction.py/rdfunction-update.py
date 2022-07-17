from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from matplotlib import pyplot as plt
import numpy as np
import warnings
from mpl_toolkits.mplot3d import Axes3D

try:
    from mayavi import mlab
    ENABLE_MAYAVI = True
except:
    ENABLE_MAYAVI = False
    warnings.warn("mayavi not found, show disabled")
    pass
import math
import csv


def nested_list_to_tuple_recur(nested_list):  # 嵌套列表转化为嵌套元组便于集合筛选出重复序列
    return tuple(
        nested_list_to_tuple_recur(l) if isinstance(l, list)
        else l for l in nested_list
    )


def points_height_matrix(smi, resolution, show=False):
    mol = Chem.MolFromSmiles(smi)
    # 在结构中添加H原子
    m3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(m3d, randomSeed=3)
    AllChem.MMFFOptimizeMolecule(m3d)  # 2D分子转为3D结构，并用力学函数进行角度校正及优化
    moleblock = Chem.MolToMolBlock(m3d)  # 返回分子中各个原子的三维坐标
    if show:
        Draw.ShowMol(m3d, size=(550, 550), kekulize=False)
    conf = m3d.GetConformer()
    sub = m3d.GetSubstructMatch(m3d)  # 这一步是统计你所需查找返回分子的个数，m3d：是查询分子的smile结构,我们这里使用自查询
    # 记录每一个原子的位置
    atom_position = []
    for s in sub:
        # list(conf.GetAtomPosition(s))#编号为s的原子坐标
        # m3d.GetAtoms()[s].GetSymbol()#编号为s对应的符号,str
        # print(m3d.GetAtoms()[s].GetSymbol(), list(conf.GetAtomPosition(s)))
        x = list(conf.GetAtomPosition(s))
        atom_symbol = m3d.GetAtoms()[s].GetSymbol()
        if atom_symbol == 'C':
            x.append(0.77)
        elif atom_symbol == 'H':
            x.append(0.37)
        elif atom_symbol == 'O':
            x.append(0.73)
        elif atom_symbol == 'S':
            x.append(1.02)
        elif atom_symbol == 'N':
            x.append(0.75)
        elif atom_symbol == 'Cl':
            x.append(0.99)
        else:
            raise NotImplementedError("only support C, H, O, S, N, Cl")
        atom_position.append(x)
    # print(atom_position)
    # 在三维坐标系绘制得到的原子坐标，可通过原子坐标调整三维坐标轴大小
    if show and ENABLE_MAYAVI:
        plt.figure()  # 得到画面
        ax1 = plt.axes(projection='3d')
        ax1.set_xlim(0, 10)  # X轴，横向向右方向
        ax1.set_ylim(10, 0)  # Y轴,左向与X,Z轴互为垂直
        ax1.set_zlim(0, 10)  # 竖向为Z轴
        # color1 = ['r', 'g', 'b', 'k', 'm']
        # marker1 = ['o', 'v', '1', 's', 'H']

        for x in atom_position:
            ax1.scatter(x[0], x[1], x[2], c='r', marker='o', linewidths=4)  # 用散点函数画点
        plt.show()

    # refrrnce_plane = -4  # 确定基平面为（0，0，-4）
    points = np.array(atom_position)  # 转化为数组操作

    # 计算基平面上的投影的x，y的坐标范围，确定生成的坐标点阵的区域
    x, y, z, r = points[:, 0], points[:, 1], points[:, 2], points[:, 3]
    r_max = max(r)
    x_max, x_min, y_max, y_min = int(max(x) + r_max) + 1, int(min(x) - r_max) - 1, int(max(y) + r_max) + 1, int(
        min(y) - r_max) - 1
    z_min = int(min(z) - r_max) - 1  # 方便后续对基平面进行平移操作

    # 画出3D分子的比例模型（堆积球模型）,直观感受实际大小
    if show and ENABLE_MAYAVI:
        mlab.points3d(x, y, z, r * 2, scale_factor=1, resolution=30, mode="sphere")
        mlab.outline()
        mlab.axes()
        mlab.show()

    # 建立矩阵开始记录高度
    # todo: change tuple operation to numpy operation
    points_intial = np.zeros([resolution, resolution])  # 建立网格点的初始化高度矩阵，与生成点阵坐标一一对应
    # print(x_max, x_min, y_max, y_min)
    x_axes = np.linspace(x_min, x_max, resolution)
    y_axes = np.linspace(y_min, y_max, resolution)
    X, Y = np.meshgrid(x_axes, y_axes)
    coordinate_points = np.array([X.ravel(), Y.ravel()])
    # print(coordinate_points)
    coordinate_points_array = coordinate_points.T
    # print(coordinate_points_array)
    coordinate_points_array = np.array(coordinate_points_array)  # 这里得再去定义一下类型，我也不知道为什么，不然后面无法操作
    coordinate_points_list = coordinate_points_array.tolist()
    coordinate_points_tuple = nested_list_to_tuple_recur(coordinate_points_list)  # 转换为嵌套元组
    v = [[] for i in range(resolution)]
    for i in range(resolution):
        for j in range(resolution):
            v[i].append(coordinate_points_tuple[resolution * i + j])

    # todo: calculate per pixel label of each atom in the image
    # 开始计算分子中各个原子距离基平面的高度，此时我们基平面设为xoy平面，虚拟为通过点光源簇从上向下进行投影，即以接触到的最高的原子的高度作为记录的高度
    for r in range(resolution):
        for t in range(resolution):
            x1 = v[r][t][0]
            y1 = v[r][t][1]  # 生成的坐标点的坐标
            height = []
            for a in atom_position:
                x2, y2, r2 = a[0], a[1], a[3]
                z2 = a[2] + abs(z_min)  # 我们将基平面平移至最底部，基平面此时的方程为：z=z_min
                changdu = math.sqrt(math.pow((x1 - x2), 2) + math.pow((y1 - y2), 2))
                if changdu <= r2:
                    point_height = math.sqrt(r2 * r2 - changdu * changdu) + z2
                    height.append(point_height)
                else:
                    continue
            if len(height) == 0:
                continue
            else:
                point_height_max = max(height)
                # print(point_height_max)
                points_intial[r][t] = point_height_max
    if show:
        plt.imshow(points_intial)
        plt.show()
    with open("points_height.csv", "w+", newline='', encoding='GBK') as f:
        writer = csv.writer(f, delimiter=',')
        for i in points_intial:  # 对于每一行的，将这一行的每个元素分别写在对应的列中
            writer.writerow(i)


if __name__ == "__main__":
    # 由于原子半径的限制，目前只支持以下例子
    smi = 'CNC(=O)N(N(CCCl)S(C)(=O)=O)S(C)(=O)=O'
    resolution = 400
    points_height_matrix(smi, resolution, show=True)
