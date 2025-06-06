import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.patches import Circle, Rectangle
import casadi as cas



# Connect the slider to the update function
global mass_force, muscle_length_isovolume_old, muscle_width_isovolume_old, perpendicular_muscle_force_old, parallel_muscle_force_old, theta_isovolume_old, muscle_fiber_length_isovolume_old
muscle_length_isovolume_old = None
muscle_width_isovolume_old = None
perpendicular_muscle_force_old = None
parallel_muscle_force_old = None
theta_isovolume_old = None
muscle_fiber_length_isovolume_old = None
# 5.5 * np.cos(np.arctan(2*tendon_length / (muscle_length-2*attachment)))  # So that the initial configuration work
mass_force = 5.444722215136416


def degroote_equation(muscle_force):
    # degroote_force = np.cos(theta_degroote) * muscle_force
    # The muscle length making the mass_force = degroote_force
    muscle_fiber_length_degroote = 2 * tendon_length / np.tan(np.arccos(mass_force / muscle_force)) + 2 * attachment
    theta_degroote = np.arctan(2 * tendon_length / (muscle_fiber_length_degroote - 2 * attachment))
    return muscle_fiber_length_degroote, theta_degroote


def isovolume_equation(muscle_force, muscle_length_isovolume_old, muscle_width_isovolume_old,
                       perpendicular_muscle_force_old, parallel_muscle_force_old, theta_isovolume_old, muscle_fiber_length_isovolume_old):
    volume = muscle_length * tendon_length  # Initial configuration

    muscle_length_isovolume = cas.SX.sym("muscle_length_isovolume", 1)
    muscle_width_isovolume = cas.SX.sym("muscle_width_isovolume", 1)
    perpendicular_muscle_force = cas.SX.sym("perpendicular_muscle_force", 1)
    parallel_muscle_force = cas.SX.sym("parallel_muscle_force", 1)
    theta_isovolume = cas.SX.sym("theta_isovolume", 1)
    muscle_fiber_length_isovolume = cas.SX.sym("muscle_fiber_length_isovolume", 1)
    w = cas.vertcat(muscle_length_isovolume, muscle_width_isovolume, perpendicular_muscle_force,
                    parallel_muscle_force, theta_isovolume,
                    muscle_fiber_length_isovolume)

    g = []
    g += [muscle_length_isovolume * muscle_width_isovolume - volume]
    g += [parallel_muscle_force - cas.cos(theta_isovolume) * muscle_force]
    g += [perpendicular_muscle_force / muscle_width_isovolume - parallel_muscle_force / muscle_length_isovolume]
    g += [perpendicular_muscle_force - cas.sin(theta_isovolume) * muscle_force]
    g += [muscle_length_isovolume - 2 * 0.05 + cas.cos(theta_isovolume) * muscle_fiber_length_isovolume]
    g += [muscle_width_isovolume - cas.sin(theta_isovolume) * muscle_fiber_length_isovolume]

    constraints = cas.Function("constraints", [w], [cas.vertcat(*g)])
    opts = {
        "print_level": 3
    }
    root_finder = cas.rootfinder('G', 'newton', constraints, opts)

    if muscle_length_isovolume_old is None:
        muscle_length_init = 0.8
        muscle_width_init = 0.1
        theta_init = np.arctan(2*tendon_length / (muscle_length_init - 2*attachment))
        parallel_muscle_force_init = muscle_force * np.cos(theta_init)
        perpendicular_muscle_force_init = muscle_force * np.sin(theta_init)
        muscle_fiber_length_init = muscle_width_init / np.sin(theta_init)
        w0 = [muscle_length_init,
              muscle_width_init,
              perpendicular_muscle_force_init,
              parallel_muscle_force_init,
              theta_init,
              muscle_fiber_length_init]
    else:
        w0 = [muscle_length_isovolume_old, muscle_width_isovolume_old, perpendicular_muscle_force_old,
              parallel_muscle_force_old,
              theta_isovolume_old, muscle_fiber_length_isovolume_old]
    w_opt = root_finder(w0)

    muscle_length_isovolume_opt = float(w_opt[0])
    muscle_width_isovolume_opt = float(w_opt[1])
    perpendicular_muscle_force_opt = float(w_opt[2])
    parallel_muscle_force_opt = float(w_opt[3])
    theta_isovolume_opt = float(w_opt[4])
    muscle_fiber_length_isovolume_opt = float(w_opt[5])

    return muscle_length_isovolume_opt, muscle_width_isovolume_opt, perpendicular_muscle_force_opt, parallel_muscle_force_opt, theta_isovolume_opt, muscle_fiber_length_isovolume_opt



# Create a figure and axis
fig, axs = plt.subplots(2, 2, figsize=(7, 7))
plt.subplots_adjust(bottom=0.25)
axs = axs.ravel()

# Set up the subplots
axs[0].set_title('DeGroote')
axs[1].set_title('IsoVolume')
axs[0].set_ylim(-1.1, 0)
axs[1].set_ylim(-1.1, 0)
axs[0].set_xlim(-0.5, 0.5)
axs[1].set_xlim(-0.5, 0.5)
axs[0].set_aspect('equal')
axs[1].set_aspect('equal')
axs[2].set_xlabel("Muscle force")
axs[2].set_ylabel("Muscle fiber length")
axs[2].set_ylim(0.2, 1.0)
fig.delaxes(axs[3])

# Not moving components
mass_radius = 0.05
tendon_length = 0.1
attachment = 0.05
muscle_length = 0.8
vertical_length = 0.5
axs[0].plot(np.array([0, 0]), np.array([0, -tendon_length]), '-k')
axs[0].plot(np.array([-tendon_length, 0]), np.array([-tendon_length, -tendon_length]), '-k')
axs[0].plot(np.array([-tendon_length, -tendon_length]), np.array([-tendon_length, -tendon_length-muscle_length]), '-k')
axs[1].plot(np.array([0, 0]), np.array([0, -0.1]), '-k')

x_array = np.linspace(10, 25, 30)
y_degroote = np.zeros((30, ))
y_isovolume = np.zeros((30, ))
for i in range(30):
    y_degroote[i], _ = degroote_equation(x_array[i])
    if i == 0:
        (muscle_length_isovolume_old, muscle_width_isovolume_old, perpendicular_muscle_force_old, parallel_muscle_force_old, theta_isovolume_old, y_isovolume[i]) = isovolume_equation(x_array[i], muscle_length_isovolume_old, muscle_width_isovolume_old, perpendicular_muscle_force_old, parallel_muscle_force_old, theta_isovolume_old, 0.8)
    else:
        (muscle_length_isovolume_old, muscle_width_isovolume_old, perpendicular_muscle_force_old, parallel_muscle_force_old, theta_isovolume_old, y_isovolume[i]) = isovolume_equation(x_array[i], muscle_length_isovolume_old, muscle_width_isovolume_old, perpendicular_muscle_force_old, parallel_muscle_force_old, theta_isovolume_old, 0.8)

    muscle_fiber_length_isovolume_old = y_isovolume[i]

axs[2].plot(x_array, y_degroote, "-g", label="DeGroote")
axs[2].plot(x_array, y_isovolume, "-m", label="IsoVolume")
axs[2].legend(bbox_to_anchor=(1.05, 0.5))


# Moving components
muscle_width = 2*tendon_length
muscle_top_isovolume = axs[1].plot(np.array([-muscle_width/2, 0]), np.array([-tendon_length, -tendon_length]), '-k')
muscle_left_isovolume = axs[1].plot(np.array([-muscle_width/2, -muscle_width/2]), np.array([-tendon_length, -tendon_length-vertical_length]), '-k')
muscle_right_isovolume = axs[1].plot(np.array([muscle_width/2, muscle_width/2]), np.array([-tendon_length, -tendon_length-muscle_length]), '-k')
muscle_bottom_isovolume = axs[1].plot(np.array([0, muscle_width/2]), np.array([-muscle_length-tendon_length, -muscle_length-tendon_length]), '-k')
tendon_bottom_isovolume = axs[1].plot(np.array([0, 0]), np.array([-muscle_length-2*tendon_length, -muscle_length-tendon_length]), '-k')
muscle_isovolume = axs[1].plot(np.array([-muscle_width/2, muscle_width/2]), np.array([-tendon_length-attachment, -tendon_length-muscle_length+attachment]), '-r')

muscle_right_degroote = axs[0].plot(np.array([tendon_length, tendon_length]), np.array([-tendon_length-muscle_length, -tendon_length-muscle_length+vertical_length]), '-k')
muscle_bottom_degroote = axs[0].plot(np.array([0, tendon_length]), np.array([-tendon_length-muscle_length, -tendon_length-muscle_length]), '-k')
tendon_bottom_degroote = axs[0].plot(np.array([0, 0]), np.array([-muscle_length-2*tendon_length, -muscle_length-tendon_length]), '-k')
muscle_degroote = axs[0].plot(np.array([-tendon_length, tendon_length]), np.array([-tendon_length-attachment, -tendon_length-muscle_length+attachment]), '-r')

mass_degroote = Circle((0, -muscle_length-3*tendon_length), mass_radius, edgecolor='blue', facecolor='lightblue', alpha=0.6)
axs[0].add_patch(mass_degroote)

mass_isovolume = Circle((0, -muscle_length-3*tendon_length), mass_radius, edgecolor='blue', facecolor='lightblue', alpha=0.6)
axs[1].add_patch(mass_isovolume)
rectangle_isovolume = Rectangle((-muscle_width/2, -tendon_length), width=muscle_width, height=-muscle_length, color='black', alpha=0.2)
axs[1].add_patch(rectangle_isovolume)

vertical_line = axs[2].plot(np.array([10, 10]), np.array([-0, 1]), "--k")

# Set up the slider
ax_slider = plt.axes([0.25, 0.01, 0.65, 0.03], facecolor='lightgoldenrodyellow')
slider = Slider(ax_slider, 'Muscle force', 10.0, 25.0, valinit=10.0, valstep=0.1)

# Update function for the slider
def update(val):
    global muscle_length_isovolume_old, muscle_width_isovolume_old, perpendicular_muscle_force_old, parallel_muscle_force_old, theta_isovolume_old, muscle_fiber_length_isovolume_old

    def update_visual_degroote(muscle_fiber_length_degroote, theta_degroote):
        vertical_length_degroote = np.cos(theta_degroote) * muscle_fiber_length_degroote
        muscle_right_degroote[0].set_ydata(np.array([-tendon_length-2*attachment-vertical_length_degroote, -tendon_length-2*attachment]))
        muscle_bottom_degroote[0].set_ydata(np.array([-tendon_length-2*attachment-vertical_length_degroote,
                                                   -tendon_length-2*attachment-vertical_length_degroote]))
        tendon_bottom_degroote[0].set_ydata(np.array([-tendon_length-2*attachment-vertical_length_degroote,
                                                   -tendon_length-2*attachment-vertical_length_degroote-tendon_length]))
        muscle_degroote[0].set_ydata(np.array([-tendon_length - attachment,
                                            -tendon_length-attachment-vertical_length_degroote]))

        mass_degroote.set_center((0, -3*tendon_length-2*attachment-vertical_length_degroote))
        return

    def update_visual_isovolume(muscle_length_isovolume, muscle_width_isovolume):

        muscle_top_isovolume[0].set_xdata(np.array([-muscle_width_isovolume / 2, 0]))
        muscle_left_isovolume[0].set_xdata(np.array([-muscle_width_isovolume / 2, -muscle_width_isovolume / 2]))
        muscle_right_isovolume[0].set_xdata(np.array([muscle_width_isovolume / 2, muscle_width_isovolume / 2]))
        muscle_right_isovolume[0].set_ydata(np.array([-tendon_length -muscle_length_isovolume, -tendon_length -muscle_length_isovolume+vertical_length]))
        muscle_bottom_isovolume[0].set_xdata(np.array([0, muscle_width_isovolume / 2]))
        muscle_bottom_isovolume[0].set_ydata(np.array([-tendon_length -muscle_length_isovolume, -tendon_length -muscle_length_isovolume]))
        tendon_bottom_isovolume[0].set_ydata(np.array([-muscle_length_isovolume - 2 * tendon_length, -muscle_length_isovolume - tendon_length]))
        muscle_isovolume[0].set_xdata(np.array([-muscle_width_isovolume / 2, muscle_width_isovolume / 2]))
        muscle_isovolume[0].set_ydata(np.array([-tendon_length - attachment, -tendon_length - muscle_length_isovolume + attachment]))

        mass_isovolume.set_center((0, -muscle_length_isovolume - 3 * tendon_length))
        rectangle_isovolume.set_xy((-muscle_width_isovolume / 2, -tendon_length))
        rectangle_isovolume.set_width(muscle_width_isovolume)
        rectangle_isovolume.set_height(-muscle_length_isovolume)

        return


    muscle_force = slider.val

    muscle_fiber_length_degroote, theta_degroote = degroote_equation(muscle_force)
    update_visual_degroote(muscle_fiber_length_degroote, theta_degroote)

    muscle_length_isovolume_opt, muscle_width_isovolume_opt, perpendicular_muscle_force_opt, parallel_muscle_force_opt, theta_isovolume_opt, muscle_fiber_length_isovolume_opt = isovolume_equation(muscle_force, muscle_length_isovolume_old, muscle_width_isovolume_old, perpendicular_muscle_force_old, theta_isovolume_old, muscle_fiber_length_isovolume_old)
    muscle_length_isovolume_old = muscle_length_isovolume_opt
    muscle_width_isovolume_old = muscle_width_isovolume_opt
    perpendicular_muscle_force_old = perpendicular_muscle_force_opt
    parallel_muscle_force_old = parallel_muscle_force_opt
    theta_isovolume_old = theta_isovolume_opt
    muscle_fiber_length_isovolume_old = muscle_fiber_length_isovolume_opt

    update_visual_isovolume(muscle_length_isovolume_opt, muscle_width_isovolume_opt)

    vertical_line[0].set_xdata(np.array([muscle_force, muscle_force]))

    fig.canvas.draw_idle()  # Redraw the canvas


slider.on_changed(update)

# Show the plot
plt.show()


