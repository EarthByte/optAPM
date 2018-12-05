# Copyright (C) 2016  Michael G. Tetley
# EarthByte Group, University of Sydney / Seismological Laboratory, California Institute of Technology
# Contact email: michael.tetley@sydney.edu.au / mtetley@caltech.edu
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


# Required libraries
import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.basemap import Basemap


class Plotting():


    @staticmethod
    def plotBlockMetrics(misfit_lat, misfit_lon, misfit_rot, misfit_vel, misfit_azi, poles_age, poles_A95, timesteps, time_start, title):

        plt.figure(figsize=(15,8),dpi=150)

        plt.xlabel('Age [Myr]', size=20)
        plt.ylabel('Misfit' + u'\xb0', size=20)

        plt.axhline(0, color='k', linewidth=3)
        plt.plot(poles_age, np.zeros(len(poles_age)), 'o-', markersize=8, color='k', linewidth=3, label='Paleomag model')
        plt.errorbar(poles_age, np.zeros(len(poles_age)), yerr=poles_A95, linestyle="None", linewidth='1', color='k', 
                     alpha=0.7, label="Pmag uncertainty")

        plt.plot(timesteps, misfit_lat, color='b', label='Latitude')
        plt.plot(timesteps, misfit_lon, color='r', label='Longitude')

        plt.xticks(np.arange(0, time_start, 10))
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.xlim(0, time_start)
        plt.gca().xaxis.grid(True)

        plt.title(str(title) + ' palaeomagnetic misfit', size=20)

        plt.plot(timesteps, -(np.array(misfit_rot)), linestyle='dashed', linewidth=2, color='g', label='Rotation angle')
        plt.plot(timesteps, misfit_azi, linestyle='solid', linewidth=1, color='orange', label='Velocity azimuth')

        ylim = plt.ylim()

        plt.legend(bbox_to_anchor=(1.32, 1), loc=1, borderaxespad=0., fontsize=16)

        plt2 = plt.twinx()
        plt2.set_ylabel('Misfit [cm/yr]', color='purple', fontsize=20)
        plt2.tick_params(axis='y', which='major', labelsize=12)

        plt2.plot(timesteps, misfit_vel, linestyle='solid', linewidth=1, color='purple', label='Velocity magnitude')

        plt2.set_ylim(ylim[0] / 2, ylim[1] / 2)

        for tl in plt2.get_yticklabels():
            tl.set_color('purple')

        plt2.legend(bbox_to_anchor=(1.336, 0.66), borderaxespad=0., fontsize=16)

        plt.show()
