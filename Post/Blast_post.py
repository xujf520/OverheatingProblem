import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def read_tecplot_data(filename):
    """
    读取 Tecplot 格式的二维数据文件
    """
    print(f"正在读取文件: {filename}")
    
    # 读取文件头部信息
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 解析参数
    I, J = 200, 200  # 从文件中得知的网格尺寸
    data_start = 0
    
    # 找到数据开始位置
    for i, line in enumerate(lines):
        if line.strip().startswith('0.0'):  # 数据行以数字开始
            data_start = i
            break
    
    # 读取数据
    data = []
    for line in lines[data_start:]:
        values = list(map(float, line.strip().split()))
        data.append(values)
    
    data = np.array(data)
    
    # 提取各变量
    x = data[:, 0].reshape(J, I)
    y = data[:, 1].reshape(J, I)
    rho = data[:, 2].reshape(J, I)
    u = data[:, 3].reshape(J, I)
    v = data[:, 4].reshape(J, I)
    p = data[:, 5].reshape(J, I)
    T = data[:, 6].reshape(J, I)
    rhou = data[:, 7].reshape(J, I)
    rhov = data[:, 8].reshape(J, I)
    rhoE = data[:, 9].reshape(J, I)
    
    return x, y, rho, u, v, p, T, rhou, rhov, rhoE

def plot_3d_surface(x, y, z, title, zlabel, colormap='viridis'):
    """
    创建三维曲面图
    """
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # 创建曲面图
    surf = ax.plot_surface(x, y, z, cmap=colormap, 
                          alpha=0.8, linewidth=0, 
                          antialiased=True)
    
    # 设置标签
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_zlabel(zlabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # 添加颜色条
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label=zlabel)
    
    # 设置视角
    ax.view_init(elev=30, azim=45)
    
    return fig, ax

def plot_contour_surface(x, y, z, title, zlabel, colormap='jet'):
    """
    创建等高线填充的三维曲面图
    """
    fig = plt.figure(figsize=(14, 10))
    
    # 三维曲面图
    ax1 = fig.add_subplot(121, projection='3d')
    surf1 = ax1.plot_surface(x, y, z, cmap=colormap, 
                            alpha=0.8, linewidth=0)
    ax1.set_xlabel('X', fontsize=12)
    ax1.set_ylabel('Y', fontsize=12)
    ax1.set_zlabel(zlabel, fontsize=12)
    ax1.set_title(f'3D Surface - {title}', fontsize=12)
    fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10)
    ax1.view_init(elev=30, azim=45)
    
    # 二维等高线图
    ax2 = fig.add_subplot(122)
    contour = ax2.contourf(x, y, z, 20, cmap=colormap)
    ax2.set_xlabel('X', fontsize=12)
    ax2.set_ylabel('Y', fontsize=12)
    ax2.set_title(f'Contour - {title}', fontsize=12)
    ax2.set_aspect('equal')
    fig.colorbar(contour, ax=ax2, shrink=0.8, aspect=10, label=zlabel)
    
    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    return fig

def create_interactive_plots(x, y, variables, titles, labels):
    """
    创建交互式多变量可视化
    """
    fig = plt.figure(figsize=(16, 12))
    
    # 选择要显示的变量（密度、压力、速度、温度）
    selected_vars = [
        (variables[0], titles[0], labels[0]),  # 密度
        (variables[5], titles[5], labels[5]),  # 压力
        (np.sqrt(variables[3]**2 + variables[4]**2), 'Velocity Magnitude', '|V|'),  # 速度大小
        (variables[6], titles[6], labels[6]),  # 温度
    ]
    
    for idx, (z, title, zlabel) in enumerate(selected_vars, 1):
        ax = fig.add_subplot(2, 2, idx, projection='3d')
        surf = ax.plot_surface(x, y, z, cmap='plasma' if idx%2==0 else 'viridis',
                              alpha=0.85, linewidth=0.1)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel(zlabel)
        ax.set_title(title, fontsize=11)
        ax.view_init(elev=30, azim=45)
        fig.colorbar(surf, ax=ax, shrink=0.6, aspect=10)
    
    plt.suptitle('BLAST HLL - 2D Fluid Dynamics Variables', fontsize=16, fontweight='bold')
    plt.tight_layout()
    return fig

# 主程序
def main():
    # 文件路径
    filename = "/mnt/d/Desktop/RP_FVM/data/BLAST_HLL_output_data_0.038000.plt"
    
    # 检查文件是否存在
    if not os.path.exists(filename):
        print(f"错误：文件不存在 - {filename}")
        return
    
    try:
        # 读取数据
        print("读取数据中...")
        x, y, rho, u, v, p, T, rhou, rhov, rhoE = read_tecplot_data(filename)
        
        # 准备变量列表
        variables = [rho, u, v, p, T, rhou, rhov, rhoE]
        titles = [
            'Density (ρ)',
            'X-Velocity (u)',
            'Y-Velocity (v)',
            'Pressure (p)',
            'Temperature (T)',
            'X-Momentum (ρu)',
            'Y-Momentum (ρv)',
            'Total Energy (ρE)'
        ]
        labels = ['ρ [kg/m³]', 'u [m/s]', 'v [m/s]', 'p [Pa]', 'T [K]', 'ρu', 'ρv', 'ρE']
        
        print(f"数据形状: X={x.shape}, Y={y.shape}")
        print(f"数据范围: X=[{x.min():.4f}, {x.max():.4f}], Y=[{y.min():.4f}, {y.max():.4f}]")
        print(f"密度范围: [{rho.min():.4f}, {rho.max():.4f}]")
        print(f"压力范围: [{p.min():.4f}, {p.max():.4f}]")
        
        # 1. 显示密度和压力的详细对比
        print("\n创建密度和压力对比图...")
        fig1 = plot_contour_surface(x, y, rho, 
                                   'Density Distribution at t=0.038s', 
                                   'Density [kg/m³]', 
                                   colormap='viridis')
        
        fig2 = plot_contour_surface(x, y, p, 
                                   'Pressure Distribution at t=0.038s', 
                                   'Pressure [Pa]', 
                                   colormap='hot')
        
        # 2. 创建交互式多变量图
        print("创建多变量交互图...")
        fig3 = create_interactive_plots(x, y, variables, titles, labels)
        
        # 3. 创建速度场可视化
        print("创建速度场可视化...")
        fig4 = plt.figure(figsize=(15, 6))
        
        # 速度大小
        vel_magnitude = np.sqrt(u**2 + v**2)
        
        ax1 = fig4.add_subplot(121, projection='3d')
        surf1 = ax1.plot_surface(x, y, vel_magnitude, cmap='coolwarm', 
                                alpha=0.8, linewidth=0)
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('Velocity Magnitude [m/s]')
        ax1.set_title('Velocity Magnitude Distribution')
        ax1.view_init(elev=30, azim=45)
        fig4.colorbar(surf1, ax=ax1, shrink=0.6, aspect=10)
        
        # 矢量图投影
        ax2 = fig4.add_subplot(122)
        # 下采样以避免过于密集的箭头
        stride = 5
        ax2.quiver(x[::stride, ::stride], y[::stride, ::stride], 
                  u[::stride, ::stride], v[::stride, ::stride],
                  vel_magnitude[::stride, ::stride], cmap='coolwarm',
                  scale=30, width=0.003, headwidth=3)
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.set_title('Velocity Vector Field (colored by magnitude)')
        ax2.set_aspect('equal')
        fig4.colorbar(plt.cm.ScalarMappable(cmap='coolwarm'), ax=ax2, 
                     label='Velocity Magnitude [m/s]')
        
        plt.suptitle('Velocity Field Analysis', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        # 显示所有图形
        print("\n显示图形中...")
        plt.show()
        
        print("\n数据处理和可视化完成！")
        
        # 保存数据统计
        print("\n=== 数据统计 ===")
        for title, var in zip(titles, variables):
            print(f"{title}:")
            print(f"  最小值: {var.min():.6f}")
            print(f"  最大值: {var.max():.6f}")
            print(f"  平均值: {var.mean():.6f}")
            print(f"  标准差: {var.std():.6f}")
            print()
            
    except Exception as e:
        print(f"处理数据时发生错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()