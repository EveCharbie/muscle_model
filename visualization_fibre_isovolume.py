import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.patches import Circle


# Connect the slider to the update function
# 5.5 * np.cos(np.arctan(2*tendon_length / (muscle_length-2*attachment)))  # So that the initial configuration work
mass_force = 5.288381712024528


def degroote_equation(muscle_force):
    # degroote_force = np.cos(theta_degroote) * muscle_force
    # The muscle length making the mass_force = degroote_force
    muscle_fiber_length_degroote = 2 * tendon_length / np.tan(np.arccos(mass_force / muscle_force)) + 2 * attachment
    theta_degroote = np.arctan(2 * tendon_length / (muscle_fiber_length_degroote - 2 * attachment))
    return muscle_fiber_length_degroote, theta_degroote


# def isovolume_equation(muscle_force_normalize):
#     """
#     Inverse behavior from what I want
#     """
#     volume = 0.5 * 0.01  # Initial configuration
#     muscle_force = muscle_force_normalize / 2
#     coeff = 0.1
#
#     # muscle_force / muscle_fiber_length_isovolume = Fw * coeff / muscle_fiber_width_isovolume
#     # phi = np.pi/2 - theta_isovolume
#     # muscle_fiber_length_isovolume * muscle_fiber_width_isovolume = volume
#     # np.sin(theta_isovolume) * muscle_force = np.sin(phi) * Fw
#     # np.cos(theta_isovolume) + np.cos(phi) * Fw = mass_force
#
#     theta_isovolume = np.arccos(muscle_force / mass_force)
#     muscle_fiber_length_isovolume = np.sqrt(volume / (np.tan(theta_isovolume) * coeff))
#     muscle_fiber_width_isovolume = volume / muscle_fiber_length_isovolume
#     phi = np.pi/2 - theta_isovolume
#
#     return muscle_fiber_length_isovolume, muscle_fiber_width_isovolume, theta_isovolume


def isovolume_equation(muscle_force_normalize, volume):
    muscle_force = muscle_force_normalize * 0.8
    if volume is None:
        volume = 0.08  # Initial configuration

    theta_isovolume = np.arccos(muscle_force / mass_force)
    muscle_fiber_length_isovolume = np.sqrt(volume / np.tan(theta_isovolume))
    muscle_fiber_width_isovolume = volume / muscle_fiber_length_isovolume

    return muscle_fiber_length_isovolume, muscle_fiber_width_isovolume, theta_isovolume


def isovolume_equation(muscle_force, volume):
    """
    Do not consider the vertical component of Fw since it is applied in both directions.
    """
    if volume is None:
        volume = 0.08  # Initial configuration

    theta_isovolume = np.arccos(mass_force / muscle_force)
    muscle_fiber_length_isovolume = np.sqrt(volume / np.tan(theta_isovolume))
    muscle_fiber_width_isovolume = volume / muscle_fiber_length_isovolume

    return muscle_fiber_length_isovolume, muscle_fiber_width_isovolume, theta_isovolume

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
axs[2].set_ylim(0.0, 1.0)
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

x_array = np.linspace(5.5, 10, 30)
y_degroote = np.zeros((30, ))
y_isovolume = np.zeros((30, ))
for i in range(30):
    y_degroote[i], _ = degroote_equation(x_array[i])
    y_isovolume[i], _, _ = isovolume_equation(x_array[i], None)

axs[2].plot(x_array, y_degroote, "-g", label="DeGroote")
plot_isovolume = axs[2].plot(x_array, y_isovolume, "-m", label="IsoVolume")
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

vertical_line = axs[2].plot(np.array([5.5, 5.5]), np.array([-0, 1]), "--k")

# Set up the slider
ax_slider_1 = plt.axes([0.25, 0.01, 0.65, 0.03], facecolor='lightgoldenrodyellow')
ax_slider_2 = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor='lightgoldenrodyellow')
slider_force = Slider(ax_slider_1, 'Muscle force', 5.0, 10.0, valinit=5.5, valstep=0.1)
slider_volume = Slider(ax_slider_2, 'Volume', 0.0001, 0.2, valinit=0.08, valstep=0.01)

# Update function for the slider
def update(val):

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

    def update_visual_isovolume(muscle_fiber_length_isovolume, muscle_fiber_width_isovolume, theta_isovolume, volume):
        muscle_length_isovolume = 2*attachment + muscle_fiber_length_isovolume * np.cos(theta_isovolume)
        muscle_width_isovolume = muscle_fiber_length_isovolume * np.sin(theta_isovolume)

        muscle_top_isovolume[0].set_xdata(np.array([-muscle_width_isovolume / 2, 0]))
        muscle_left_isovolume[0].set_xdata(np.array([-muscle_width_isovolume / 2, -muscle_width_isovolume / 2]))
        muscle_right_isovolume[0].set_xdata(np.array([muscle_width_isovolume / 2, muscle_width_isovolume / 2]))
        muscle_right_isovolume[0].set_ydata(np.array([-tendon_length -muscle_length_isovolume, -tendon_length -muscle_length_isovolume+vertical_length]))
        muscle_bottom_isovolume[0].set_xdata(np.array([0, muscle_width_isovolume / 2]))
        muscle_bottom_isovolume[0].set_ydata(np.array([-tendon_length -muscle_length_isovolume, -tendon_length -muscle_length_isovolume]))
        tendon_bottom_isovolume[0].set_ydata(np.array([-muscle_length_isovolume - 2 * tendon_length, -muscle_length_isovolume - tendon_length]))
        muscle_isovolume[0].set_xdata(np.array([-muscle_width_isovolume / 2, muscle_width_isovolume / 2]))
        muscle_isovolume[0].set_ydata(np.array([-tendon_length - attachment, -tendon_length - muscle_length_isovolume + attachment]))
        muscle_isovolume[0].set_linewidth(muscle_fiber_width_isovolume * 10)

        mass_isovolume.set_center((0, -muscle_length_isovolume - 3 * tendon_length))

        x_array = np.linspace(5.5, 10, 30)
        y_isovolume = np.zeros((30,))
        for i in range(30):
            y_isovolume[i], _, _ = isovolume_equation(x_array[i], volume)
        plot_isovolume[0].set_ydata(y_isovolume)

        return


    muscle_force = slider_force.val
    volume = slider_volume.val

    muscle_fiber_length_degroote, theta_degroote = degroote_equation(muscle_force)
    update_visual_degroote(muscle_fiber_length_degroote, theta_degroote)

    muscle_fiber_length_isovolume, muscle_fiber_width_isovolume, theta_isovolume = isovolume_equation(muscle_force, volume)
    update_visual_isovolume(muscle_fiber_length_isovolume, muscle_fiber_width_isovolume, theta_isovolume, volume)

    vertical_line[0].set_xdata(np.array([muscle_force, muscle_force]))

    fig.canvas.draw_idle()  # Redraw the canvas


slider_force.on_changed(update)
slider_volume.on_changed(update)

# Show the plot
plt.show()


