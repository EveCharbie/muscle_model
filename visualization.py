import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.patches import Circle, Rectangle
import casadi as cas

# Create a figure and axis
fig, axs = plt.subplots(1, 2)
# plt.subplots_adjust(bottom=0.25)

# Set up the subplots
axs[0].set_title('DeGroote')
axs[1].set_title('IsoVolume')
axs[0].set_ylim(-1.1, 0)
axs[1].set_ylim(-1.1, 0)
axs[0].set_xlim(-0.5, 0.5)
axs[1].set_xlim(-0.5, 0.5)
axs[0].set_aspect('equal')
axs[1].set_aspect('equal')

# Not moving components
mass_radius = 0.1
tendon_length = 0.1
attachment = 0.05
axs[0].plot(np.array([0, 0]), np.array([0, -tendon_length]), '-k')
axs[0].plot(np.array([-tendon_length, 0]), np.array([-tendon_length, -tendon_length]), '-k')
axs[0].plot(np.array([-tendon_length, -tendon_length]), np.array([-tendon_length, -tendon_length-0.3]), '-k')
axs[1].plot(np.array([0, 0]), np.array([0, -0.1]), '-k')

# Moving components
muscle_width = 0.2
muscle_length = 0.5
muscle_top_isovolume = axs[1].plot(np.array([-muscle_width/2, 0]), np.array([-tendon_length, -tendon_length]), '-k')
muscle_left_isovolume = axs[1].plot(np.array([-muscle_width/2, -muscle_width/2]), np.array([-tendon_length, -tendon_length-0.3]), '-k')
muscle_right_isovolume = axs[1].plot(np.array([muscle_width/2, muscle_width/2]), np.array([-tendon_length, -tendon_length-0.5]), '-k')
muscle_bottom_isovolume = axs[1].plot(np.array([0, muscle_width/2]), np.array([-0.5-tendon_length, -0.5-tendon_length]), '-k')
tendon_bottom_isovolume = axs[1].plot(np.array([0, 0]), np.array([-0.5-2*tendon_length, -0.5-tendon_length]), '-k')
muscle_isovolume = axs[1].plot(np.array([-muscle_width/2, muscle_width/2]), np.array([-tendon_length-attachment, -tendon_length-muscle_length+attachment]), '-r')

muscle_right_degroote = axs[0].plot(np.array([tendon_length, tendon_length]), np.array([-tendon_length-0.5, -tendon_length-0.2]), '-k')
muscle_bottom_degroote = axs[0].plot(np.array([0, tendon_length]), np.array([-tendon_length-0.5, -tendon_length-0.5]), '-k')
tendon_bottom_degroote = axs[0].plot(np.array([0, 0]), np.array([-0.5-2*tendon_length, -0.5-tendon_length]), '-k')
muscle_degroote = axs[0].plot(np.array([-tendon_length, tendon_length]), np.array([-tendon_length-attachment, -tendon_length-0.5+attachment]), '-r')

mass_degroote = Circle((0, -0.5-3*tendon_length), mass_radius, edgecolor='blue', facecolor='lightblue', alpha=0.6)
axs[0].add_patch(mass_degroote)

mass_isovolume = Circle((0, -0.5-3*tendon_length), mass_radius, edgecolor='blue', facecolor='lightblue', alpha=0.6)
axs[1].add_patch(mass_isovolume)
rectangle_isovolume = Rectangle((-muscle_width/2, -tendon_length), width=muscle_width, height=-muscle_length, color='black', alpha=0.2)
axs[1].add_patch(rectangle_isovolume)

# Set up the slider
ax_slider = plt.axes([0.25, 0.01, 0.65, 0.03], facecolor='lightgoldenrodyellow')
slider = Slider(ax_slider, 'Muscle force', 5.1, 15.0, valinit=5.5, valstep=0.1)

# Update function for the slider
def update(val):

    def degroote_equation(muscle_force):
        # degroote_force = np.cos(theta_degroote) * muscle_force
        # The muscle length making the mass_force = degroote_force
        muscle_fiber_length_degroote = 2 * tendon_length / np.tan(np.arccos(mass_force / muscle_force)) + 2 * attachment
        theta_degroote = np.arctan(2 * tendon_length / (muscle_fiber_length_degroote - 2 * attachment))
        return muscle_fiber_length_degroote, theta_degroote

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

    def isovolume_equation(muscle_force):
        volume = 0.5 * 0.2  # Initial configuration

        muscle_length_isovolume = cas.SX.sym("muscle_length_isovolume", 1)
        muscle_width_isovolume = cas.SX.sym("muscle_width_isovolume", 1)
        perpendicular_muscle_force = cas.SX.sym("perpendicular_muscle_force", 1)
        theta_isovolume = cas.SX.sym("theta_isovolume", 1)
        muscle_fiber_length_isovolume = cas.SX.sym("muscle_fiber_length_isovolume", 1)
        w = cas.vertcat(muscle_length_isovolume, muscle_width_isovolume, perpendicular_muscle_force, theta_isovolume, muscle_fiber_length_isovolume)

        g = []
        g += [muscle_length_isovolume * muscle_width_isovolume - volume]
        g += [perpendicular_muscle_force / muscle_length_isovolume - mass_force / muscle_width_isovolume]
        g += [perpendicular_muscle_force - np.sin(theta_isovolume) * muscle_force]
        g += [muscle_length_isovolume - 2 * 0.05 + np.cos(theta_isovolume) * muscle_fiber_length_isovolume]
        g += [muscle_width_isovolume - np.sin(theta_isovolume) * muscle_fiber_length_isovolume]

        constraints = cas.Function("constraints", [w], [cas.vertcat(*g)])
        root_finder = cas.rootfinder('G', 'newton', constraints)
        w_opt = root_finder([1, 1, 1, 1, 1])

        muscle_length_isovolume_opt = float(w_opt[0])
        muscle_width_isovolume_opt = float(w_opt[1])
        perpendicular_muscle_force_opt = float(w_opt[2])
        theta_isovolume_opt = float(w_opt[3])
        muscle_fiber_length_isovolume_opt = float(w_opt[4])

        return muscle_fiber_length_isovolume_opt, muscle_length_isovolume_opt, muscle_width_isovolume_opt, theta_isovolume_opt

    def update_visual_isovolume(muscle_length_isovolume, muscle_width_isovolume):

        muscle_top_isovolume[0].set_xdata(np.array([-muscle_width_isovolume / 2, 0]))
        muscle_left_isovolume[0].set_xdata(np.array([-muscle_width_isovolume / 2, -muscle_width_isovolume / 2]))
        muscle_right_isovolume[0].set_xdata(np.array([muscle_width_isovolume / 2, muscle_width_isovolume / 2]))
        muscle_right_isovolume[0].set_ydata(np.array([-tendon_length -muscle_length_isovolume, tendon_length -muscle_length_isovolume+0.3]))
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
    # 5.5 * np.cos(np.arctan(2*tendon_length / (0.5-2*attachment)))  # So that the initial configuration work
    mass_force = 4.919349550499537

    muscle_fiber_length_degroote, theta_degroote = degroote_equation(muscle_force)
    update_visual_degroote(muscle_fiber_length_degroote, theta_degroote)

    muscle_fiber_length_isovolume, muscle_length_isovolume, muscle_width_isovolume, theta_isovolume = isovolume_equation(muscle_force)
    update_visual_isovolume(muscle_length_isovolume, muscle_width_isovolume)

    fig.canvas.draw_idle()  # Redraw the canvas

# Connect the slider to the update function
slider.on_changed(update)

# Show the plot
plt.show()


