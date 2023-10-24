import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import subprocess as sp
import matplotlib.colors as mcolors

nk = 1200
wc = 1.0
nf = 7
n_kappa = 21
BASIS = "RAD"
a_0 = 4
g_wc_array =  [ 0.25]
# g_wc_array = [0.103]
y_max_array = [5.0, 5.0]
# y_max_array = [5]
nticks = [6,6]
n_graph_array = [  20, 20]
# y_max_array = [0.4]
dark = False
color = True
k_shift = 0.0
both = True

def create_smooth_dark_colormap():
    # Define the color transitions for the colormap
    colors = [(0.5, 0.5, 0.5),  # Dark background color
              (0.2, 0.4, 0.8),  # Light blue
              (0.0, 0.7, 0.0),  # Green
              (0.8, 0.1, 0.1)]  # Red
    
    # Create a custom colormap using the defined color transitions
    cmap = mcolors.LinearSegmentedColormap.from_list('smooth_dark_colormap', colors, N=256)
    
    return cmap

# Create the colormap
dark_colormap = create_smooth_dark_colormap()

label = 'bright'
black = 'k'
if dark:
    plt.style.use('dark_background')
    label = 'dark'
    black = '0.95'

k_points = np.linspace(-np.pi / a_0 - k_shift, np.pi / a_0 - k_shift, nk)

file_location = "./disp/"

e_min = 0
for ijk in np.flip(range(len(g_wc_array))):
    g_wc = g_wc_array[ijk]
    y_max = y_max_array[ijk]

    # DIR = f"{file_location}gather_data/"
    IMG_DIR = f"./disp/"
    # sp.call(f"mkdir -p {DIR}", shell=True)
    #######
    NPol = nf * n_kappa

    EPol = np.zeros(( len(k_points), NPol , 2))

    if both:
        BASIS = "RAD"
        EPol[:,:,0] = np.loadtxt(f"./data/E_{BASIS}_nk{nk}_nf{nf}_nkappa{n_kappa}_gwc{g_wc:.7f}_wc{wc:.4f}_kshift{k_shift:.2f}.csv", delimiter = ',')[:,1:]
        EPol[:,:,1] = np.loadtxt(f"./data/E_{BASIS}_nk{nk}_nf{nf}_nkappa{n_kappa}_gwc{g_wc:.7f}_wc{wc:.4f}_kshift{k_shift:.2f}_color.csv", delimiter = ',')
    else:
        EPol[:,:,0] = np.loadtxt(f"./data/E_{BASIS}_nk{nk}_nf{nf}_nkappa{n_kappa}_gwc{g_wc:.7f}_wc{wc:.4f}_kshift{k_shift:.2f}.csv", delimiter = ',')[:,1:]
        EPol[:,:,1] = np.loadtxt(f"./data/E_{BASIS}_nk{nk}_nf{nf}_nkappa{n_kappa}_gwc{g_wc:.7f}_wc{wc:.4f}_kshift{k_shift:.2f}_color.csv", delimiter = ',')

    # plt.style.use('dark_background')
    fs = 16
    n_states = n_graph_array[ijk]
    plot_states = np.arange( n_states )
    fig, ax = plt.subplots()

    if color:
        cols = np.ndarray.flatten(((EPol[1:,:n_states,1] + EPol[:-1,:n_states,1]) / 2 ).T)

        for lmn in range(n_states):
            x = k_points
            y = EPol[:,lmn,0] - e_min

            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            if lmn > 0:
                all_segments=np.concatenate((all_segments, segments), axis=0)
            else:
                all_segments = segments

        for lmn in range(all_segments.shape[0]):
            if all_segments[lmn,0,1] > y_max:
                cols[lmn] = np.min(cols)

        if dark:
            lc = LineCollection(all_segments, cmap=dark_colormap)
        else:
            lc = LineCollection(all_segments, cmap='jet')
        lc.set_array(cols)
        lc.set(capstyle = "round")
        lc.set_clim(0,4)
        # lc.set_linewidth(3)
        line = ax.add_collection(lc)
        cbar = fig.colorbar(line,ax=ax)
        cbar.ax.tick_params(labelsize=fs) 
    else:
        for state in plot_states:
            plt.plot( k_points, EPol[:,state], color=black)

    if both:
        BASIS = "pA"
        n_kappa = 501
        nk = 30
        nf = 14
        NPol = nf * n_kappa
        k_points = np.linspace(-np.pi / a_0 - k_shift, np.pi / a_0 - k_shift, nk)
        E_Pol2 = np.zeros(( len(k_points), NPol , 2))
        E_Pol2[:,:,0] = np.loadtxt(f"./data/E_{BASIS}_nk{nk}_nf{nf}_nkappa{n_kappa}_gwc{g_wc:.7f}_wc{wc:.4f}_kshift{k_shift:.2f}.csv", delimiter = ',')[:,1:]
        E_Pol2[:,:,1] = np.loadtxt(f"./data/E_{BASIS}_nk{nk}_nf{nf}_nkappa{n_kappa}_gwc{g_wc:.7f}_wc{wc:.4f}_kshift{k_shift:.2f}_color.csv", delimiter = ',')
        e_min = np.min(E_Pol2[:,0,0])

        n = 1
        x = k_points[::n]
        y = E_Pol2[::n,0,0] - e_min
        z = E_Pol2[::n,0,1]
        
        for lmn in range(n_states):
            x = np.append(x,k_points[::n])
            y = np.append(y, E_Pol2[::n,lmn,0] - e_min)
            z = np.append(z, E_Pol2[::n,lmn,1])

        plt.scatter(x,y,c=z, s = 15, cmap = "jet", marker='.')
            
        # cbar = fig.colorbar(line,ax=ax)

        # cols = np.ndarray.flatten(((E_Pol2[1:,:n_states,1] + E_Pol2[:-1,:n_states,1]) / 2 ).T)

        # for lmn in range(n_states):
        #     x = k_points
        #     y = (E_Pol2[:,lmn,0] - e_min)

        #     points = np.array([x, y]).T.reshape(-1, 1, 2)
        #     segments = np.concatenate([points[:-1], points[1:]], axis=1)
        #     if lmn > 0:
        #         all_segments=np.concatenate((all_segments, segments), axis=0)
        #     else:
        #         all_segments = segments

        # for lmn in range(all_segments.shape[0]):
        #     if all_segments[lmn,0,1] > y_max:
        #         cols[lmn] = np.min(cols)

        # if dark:
        #     lc = LineCollection(all_segments, cmap=dark_colormap)
        # else:
        #     lc = LineCollection(all_segments, cmap='jet')
        # lc.set_array(cols)
        # lc.set(capstyle = "round")
        # # lc.set_linewidth(3)
        # line = ax.add_collection(lc)



    plt.ylim(0.0,y_max)
    # plt.yscale('log')
    # ax.set_ylim(top=y_max)
    plt.xlim(min(k_points),max(k_points))
    # plt.rcParams["figure.figsize"] = (2,1.5)
    plt.xlabel(r"$k = k_\beta$ (a.u.)",fontsize=fs)
    plt.yticks(fontsize = fs)#, ticks = np.linspace(0 , y_max, nticks[ijk], True))
    plt.xticks(ticks = [- np.pi / 4,0, np.pi / 4], labels = ["$-\pi/a_0$", "0", "$\pi/a_0$"],fontsize = fs)
    plt.title(f"$\gamma_0 / \omega_0 =$ {'%s' % float('%.3g' % g_wc)}",fontsize=fs)
    plt.ylabel("Energy (a.u.)",fontsize=fs)
    plt.subplots_adjust(left=0.17,
                    bottom=0.18,
                    right=0.93,
                    top=0.87,
                    wspace=0.2,
                    hspace=0.2)

    if dark:
        # plt.title(f"$g_0 / \omega_0 =$ {g_wc}",fontsize=fs, pad = 15)
        plt.ylabel("Energy (a.u.)",fontsize=fs)
        plt.xlabel("$k$",fontsize=fs)
        # cbar.set_label("$<a^+a>$", fontsize = fs)
        plt.subplots_adjust(left=0.17,
                    bottom=0.18,
                    right=0.93,
                    top=0.87,
                    wspace=0.2,
                    hspace=0.2)

    plt.savefig(f"{IMG_DIR}disp_plot_g{np.round(g_wc,3)}_nf{nf}_{label}.jpg",dpi=600)
    plt.savefig(f"{IMG_DIR}/disp_plot_g{np.round(g_wc,3)}_nf{nf}_{label}.svg",format='svg')

    # np.save( f"{IMG_DIR}plot_data_{g_wc}_{nk}_{a_0}_{n_kappa}_{nf}_{wc}.npy", EPol)

    plt.figure().clear()