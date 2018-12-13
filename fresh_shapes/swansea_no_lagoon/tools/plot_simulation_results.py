from thetis import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpi4py import MPI
from matplotlib import rc
import matplotlib as mpl


import h5py

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

comm = MPI.COMM_WORLD
size = comm.size

ticks0 = np.linspace(0, 5,11 , endpoint=True)
ticks1 = np.linspace(0, 360, 13, endpoint=True)

def plot_cotidal_chart(mesh2d, amplitude, phase, fontsize = 10):
    P1_2d = FunctionSpace(mesh2d, 'DG', 1)
    x = SpatialCoordinate(mesh2d)
    x_vector, y_vector = interpolate(x[0], Function(P1_2d)).dat.data, interpolate(x[1], Function(P1_2d)).dat.data

    x_0 = comm.gather(x_vector, root=0)
    if comm.rank == 0: print([x.shape for x in x_0])
    y_0 = comm.gather(y_vector, root=0)
    amp_0 = comm.gather(amplitude.dat.data, root=0)
    phase_0 = comm.gather(phase.dat.data, root=0)

    if comm.rank == 0:
        x_f, y_f, amp_f, phase_f = np.concatenate(x_0), np.concatenate(y_0), \
                                   np.concatenate(amp_0), np.concatenate(phase_0)
        points = np.vstack([x_f,y_f]).T
        print(points.shape)

        fig_ratio = (np.amax(y_f)-np.amin(y_f))/(np.amax(x_f)-np.amin(x_f))
        amp_f[0], amp_f[-1], phase_f[0], phase_f[-1] = -0.5, 5.5, -30,390
        triangulation = np.vstack(np.split(np.arange(x_f.shape[0]),x_f.shape[0]/3))
        #
        fig = plt.figure(figsize=(1.25 * 5 / fig_ratio ,5))
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        CS = plt.triplot(points[:,0],points[:,1], triangulation, lw =0.1, color ="black", alpha  =0.1)
        CS0 = plt.tricontourf(points[:, 0], points[:, 1], triangulation, amp_f, levels=ticks0, cmap=plt.cm.Blues,
                            linewidths = 0.75,  norm = mpl.colors.Normalize(vmin=0, vmax=5, clip=False , ))

        CS1 = plt.tricontour(points[:, 0], points[:, 1], triangulation, phase_f, levels=ticks1, cmap=plt.cm.Greys,
                             linestyles = "solid",linewidths = 3.0,
                             norm = mpl.colors.Normalize(vmin=-60, vmax=360, clip=False, ))

        plt.xticks(np.arange(200000, 1000000 + 1, 50000))
        plt.yticks(np.arange(5000000, 8000000 + 1, 50000))
        plt.xlim([np.amin(x_f) - 0.02 *(np.amax(x_f)-np.amin(x_f)),np.amax(x_f) + 0.02 *(np.amax(x_f)-np.amin(x_f))])
        plt.ylim([np.amin(y_f) - 0.05 * (np.amax(y_f) - np.amin(y_f)), np.amax(y_f) + 0.05 * (np.amax(y_f) - np.amin(y_f))])
        plt.xlabel("x (km)")
        plt.ylabel("y (km)")

        ax = plt.gca()
        # Change only ax1
        scale_x = 1e3
        scale_y = 1e3
        ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:}'.format(x / scale_x))
        ax.xaxis.set_major_formatter(ticks_x)
        ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:}'.format(x / scale_y))
        ax.yaxis.set_major_formatter(ticks_y)

        #ax.set_aspect(0.9)
        pos1 = ax.get_position()  # get the original position
        pos2 = [pos1.x0 , pos1.y0 , pos1.width/1.15 , pos1.height* 1.05]
        ax.set_position(pos2)

        cbaxes = fig.add_axes([0.92, pos1.y0, 0.012, pos1.height* 1.05])
        cbar0 = plt.colorbar(CS0, ticks=ticks0, cax=cbaxes)
        cbar0.ax.tick_params(labelsize=fontsize)
        cbar0.set_label("$M_2$ Amplitude (m)", size=fontsize)

        cbaxes2 = fig.add_axes([0.85, pos1.y0, 0.012, pos1.height* 1.05])
        cbar1 = plt.colorbar(CS1, ticks=ticks1, cax=cbaxes2)
        cbar1.ax.tick_params(labelsize=fontsize)
        cbar1.set_label("$M_2$ phase ($^o$)", size=fontsize)

        plt.setp(CS1.collections, linewidth=0.75)
        plt.savefig("cotidal.png", dpi = 800)


if __name__ == "__main__":
    import os
    import sys

    sys.path.insert(0,'../')
    import inputs.input_file_paths

    os.chdir("../")
    mesh2d = Mesh(inputs.input_file_paths.mesh_file)
    P1_2d = FunctionSpace(mesh2d, 'CG', 1)
    DG_2d = FunctionSpace(mesh2d, 'DG', 1)

    chk = DumbCheckpoint("outputs/M2_amp", mode=FILE_READ)
    M2_amp = Function(P1_2d, name="M2_amp")
    chk.load(M2_amp); chk.close()

    chk = DumbCheckpoint("outputs/M2_phase", mode=FILE_READ)
    M2_phase = Function(P1_2d, name="M2_phase")
    chk.load(M2_phase); chk.close()

    M2_amp_DG = Function(DG_2d, name="M2_amp").interpolate(M2_amp)
    M2_phase_DG = Function(DG_2d, name="M2_phase").interpolate(M2_phase)
    plot_cotidal_chart(mesh2d, M2_amp_DG, M2_phase_DG)