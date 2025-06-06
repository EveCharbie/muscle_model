import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.patches import Circle, Rectangle

# Create a figure and axis
fig, axs = plt.subplots(1, 2)
# plt.subplots_adjust(bottom=0.25)

# Set up the subplots
axs[0].set_title('Degroote')
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

    muscle_force = slider.val
    mass_force = 5

    # degroote_force = np.cos(theta_degroote) * muscle_force
    # The muscle length making the mass_force = degroote_force
    muscle_length_degroote = 2 * tendon_length / np.tan(np.arccos(mass_force / muscle_force)) + 2*attachment
    theta_degroote = np.arctan(2 * tendon_length / (muscle_length_degroote - 2 * attachment))

    def update_visual_degroote(muscle_length_degroote, theta_degroote):
        vertical_length_degroote = np.cos(theta_degroote) * muscle_length_degroote
        muscle_right_degroote[0].set_ydata(np.array([-tendon_length-2*attachment-vertical_length_degroote, -tendon_length-2*attachment]))
        muscle_bottom_degroote[0].set_ydata(np.array([-tendon_length-2*attachment-vertical_length_degroote,
                                                   -tendon_length-2*attachment-vertical_length_degroote]))
        tendon_bottom_degroote[0].set_ydata(np.array([-tendon_length-2*attachment-vertical_length_degroote,
                                                   -tendon_length-2*attachment-vertical_length_degroote-tendon_length]))
        muscle_degroote[0].set_ydata(np.array([-tendon_length - attachment,
                                            -tendon_length-attachment-vertical_length_degroote]))

        mass_degroote.set_center((0, -3*tendon_length-2*attachment-vertical_length_degroote))

    update_visual_degroote(muscle_length_degroote, theta_degroote)
    fig.canvas.draw_idle()  # Redraw the canvas

# Connect the slider to the update function
slider.on_changed(update)

# Show the plot
plt.show()


