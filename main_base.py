

import math
import xml.etree.ElementTree as ET
import re
from scipy.interpolate import interp1d, splprep, splev
from matplotlib import pyplot as plt
import numpy as np
from scipy.spatial import cKDTree
import copy
from shapely.geometry import Point, Polygon

def latlon_to_meters(lat, lon, lat_origin, lon_origin):
    # Radius of the Earth in meters
    R = 6378137  # WGS84 ellipsoid major axis
    # Convert latitude and longitude from degrees to radians
    lat = math.radians(lat)
    lon = math.radians(lon)
    lat_origin = math.radians(lat_origin)
    lon_origin = math.radians(lon_origin)
    
    # Latitude conversion to meters
    lat_m = R * (lat - lat_origin)
    
    # Longitude conversion to meters (adjusted by the latitude)
    lon_m = R * math.cos(lat_origin) * (lon - lon_origin)
    
    return lat_m, lon_m


def meters_to_latlon(lat_m, lon_m, lat_origin, lon_origin):
    # Radius of the Earth in meters (WGS84 ellipsoid major axis)
    R = 6378137
    
    # Convert the origin coordinates to radians
    lat_origin = math.radians(lat_origin)
    lon_origin = math.radians(lon_origin)
    
    # Convert meters back to latitude (in radians)
    lat = lat_m / R + lat_origin
    
    # Convert meters back to longitude (in radians, adjusted by latitude)
    lon = lon_m / (R * math.cos(lat_origin)) + lon_origin
    
    # Convert the latitude and longitude from radians back to degrees
    lat = math.degrees(lat)
    lon = math.degrees(lon)
    
    return lat, lon

def load_track(file_path):
    track_center_lat = None
    track_center_lon = None
    tree = ET.parse(file_path)
    kml = tree.getroot()
    namespace = re.match(r'{.*}', kml.tag).group(0)
    document = kml.find(f'{namespace}Document')
    name = document.find(f'{namespace}name').text
    track = {
        'name': name,
    }

    for placemark in document.iter(f'{namespace}Placemark'):
        name = placemark.find(f'{namespace}name').text
        point = placemark.find(f'{namespace}Point')
        if point is not None:
            coordinates = point.find(f'{namespace}coordinates')
            line = coordinates.text.strip()
            coordinate_lst = line.split(',')
            if line != '' and len(coordinate_lst) == 3:
                track_center_lon = float(coordinate_lst[0])
                track_center_lat = float(coordinate_lst[1])
                break

    if track_center_lon is None or track_center_lat is None:
        raise Exception("There is not track center point")
    
    track[name] = {
        'lat': track_center_lat,
        'lon': track_center_lon
    }

    for placemark in document.iter(f'{namespace}Placemark'):
        name = placemark.find(f'{namespace}name').text
        line_string = placemark.find(f'{namespace}Polygon')
        if line_string is not None:
            track[name] = {}
            for line_name in ['inner', 'outer']:
                points = []
                boundary_is = line_string.find(f'{namespace}{line_name}BoundaryIs')
                linear_ring = boundary_is.find(f'{namespace}LinearRing')
                coordinates = linear_ring.find(f'{namespace}coordinates')
                for line in coordinates.text.split('\n'):
                    line = line.strip()
                    coordinate_lst = line.split(',')
                    if line != '' and len(coordinate_lst) == 3:
                        lon = float(coordinate_lst[0])
                        lat = float(coordinate_lst[1])
                        alt = float(coordinate_lst[2])
                        y, x = latlon_to_meters(lat, lon, track_center_lat, track_center_lon)
                        points.append([x, y])
                points.append(points[0])
                track[name][line_name] = points

    for placemark in document.iter(f'{namespace}Placemark'):
        name = placemark.find(f'{namespace}name').text
        line_string = placemark.find(f'{namespace}LineString')
        if line_string is not None:
            coordinates = line_string.find(f'{namespace}coordinates')
            points = []
            for line in coordinates.text.split('\n'):
                line = line.strip()
                coordinate_lst = line.split(',')
                if line != '' and len(coordinate_lst) == 3:
                    lon = float(coordinate_lst[0])
                    lat = float(coordinate_lst[1])
                    alt = float(coordinate_lst[2])
                    y, x = latlon_to_meters(lat, lon, track_center_lat, track_center_lon)
                    points.append([x, y])
            points.append(points[0])
            track[name] = points
    return track

def save_point(file, name, point):
    file.write('    <Placemark>\n')
    file.write(f'      <name>{name}</name>\n')
    file.write('      <Point>\n')
    file.write('        <coordinates>\n')
    file.write(f'          {point['lon']},{point['lat']},0\n')
    file.write('        </coordinates>\n')
    file.write('      </Point>\n')
    file.write('    </Placemark>\n')


def save_polygon(file, center, name, inside, outside):
    track_center_lat = center['lat']
    track_center_lon = center['lon']
    file.write('    <Placemark>\n')
    file.write(f'      <name>{name}</name>\n')
    file.write('      <Polygon>\n')
    file.write('        <innerBoundaryIs>\n')
    file.write('          <LinearRing>\n')
    file.write('            <tessellate>1</tessellate>\n')
    file.write('            <coordinates>\n')
    for x, y in inside:
        lat, lon = meters_to_latlon(y, x, track_center_lat, track_center_lon)
        file.write(f'              {lon},{lat},0\n')
    file.write('            </coordinates>\n')
    file.write('          </LinearRing>\n')
    file.write('        </innerBoundaryIs>\n')
    file.write('        <outerBoundaryIs>\n')
    file.write('          <LinearRing>\n')
    file.write('            <tessellate>1</tessellate>\n')
    file.write('            <coordinates>\n')
    for x, y in outside:
        lat, lon = meters_to_latlon(y, x, track_center_lat, track_center_lon)
        file.write(f'              {lon},{lat},0\n')
    file.write('            </coordinates>\n')
    file.write('          </LinearRing>\n')
    file.write('        </outerBoundaryIs>\n')
    file.write('      </Polygon>\n')
    file.write('    </Placemark>\n')



def save_line(file, center, name, line):
    track_center_lat = center['lat']
    track_center_lon = center['lon']
    file.write('    <Placemark>\n')
    file.write(f'      <name>{name}</name>\n')
    file.write('      <LineString>\n')
    file.write('        <coordinates>\n')
    for x, y in line:
        lat, lon = meters_to_latlon(y, x, track_center_lat, track_center_lon)
        file.write(f'          {lon},{lat},0\n')
    file.write('        </coordinates>\n')
    file.write('      </LineString>\n')
    file.write('    </Placemark>\n')

def save_track(file_path, name, center, inside, outside, middle, race, base):
    with open(file_path, 'wt') as file:
        file.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        file.write('<kml xmlns="http://www.opengis.net/kml/2.2">\n')
        file.write('  <Document>\n')
        file.write(f'    <name>{name}</name>\n')
        save_point(file, 'center', center)
        save_polygon(file, center, 'track', inside, outside)
        if middle is not None:
            save_line(file, center, 'middle', middle)
        if race is not None:
            save_line(file, center, 'race', race)
        if base is not None:
            save_line(file, center, 'base', base)
        file.write('  </Document>\n')
        file.write('</kml>\n')

# Function to compute cumulative arc length along a polyline (track points)
def cumulative_arc_length(points):
    # Compute distance between consecutive points
    deltas = np.diff(points, axis=0)
    distances = np.hypot(deltas[:, 0], deltas[:, 1])
    # Cumulative sum of distances to get the arc length at each point
    return np.concatenate([[0], np.cumsum(distances)])

# Function to resample a track based on cumulative arc length
def resample_track(points, num_samples):
    # Get cumulative arc length
    arc_length = cumulative_arc_length(points)
    total_length = arc_length[-1]
    
    # Define the interpolation function for x and y coordinates
    interp_func_x = interp1d(arc_length, points[:, 0], kind='linear')
    interp_func_y = interp1d(arc_length, points[:, 1], kind='linear')
    
    # Resample points uniformly spaced along the arc length
    new_arc_length = np.linspace(0, total_length, num_samples)
    
    # Interpolated x and y coordinates at resampled positions
    resampled_x = interp_func_x(new_arc_length)
    resampled_y = interp_func_y(new_arc_length)
    
    # Combine into new resampled array
    return np.vstack((resampled_x, resampled_y)).T



def add_plot(np_points):
    # Extract the x and y coordinates
    x_coords = np_points[:, 0]  # All rows, first column
    y_coords = np_points[:, 1]  # All rows, second column

    # Plot the points
    plt.scatter(x_coords, y_coords)  # Scatter plot for individual points
    plt.plot(x_coords, y_coords)     # Line plot to connect the points


# Function to interpolate the inner track
def interpolate_inner_track(inner_track, num_points=1000):
    # Parametric spline interpolation of the inner track
    tck, u = splprep(inner_track.T, s=0)  # No smoothing
    u_fine = np.linspace(0, 1, num_points)
    inner_interp = np.array(splev(u_fine, tck)).T
    return inner_interp

# Function to find the nearest point on the interpolated inner track
def find_nearest_points(outer_track, inner_track):
    # Interpolate the inner track
    inner_interp = interpolate_inner_track(inner_track)
    
    # Create a KDTree for fast nearest-neighbor search on the interpolated inner track
    inner_tree = cKDTree(inner_interp)
    
    # Find the nearest point for each outer track point
    distances, indices = inner_tree.query(outer_track)
    
    # Get the corresponding nearest points from the interpolated inner track
    nearest_inner_points = inner_interp[indices]
    
    return nearest_inner_points, distances


# From https://github.com/e-koch/ewky_scripts/blob/master/curvature.py

# The MIT License (MIT)
#
# Copyright (c) 2014 Eric Koch
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

def menger_curvature(pt1, pt2, pt3, atol=1e-3):

    vec21 = np.array([pt1[0]-pt2[0], pt1[1]-pt2[1]])
    vec23 = np.array([pt3[0]-pt2[0], pt3[1]-pt2[1]])

    norm21 = np.linalg.norm(vec21)
    norm23 = np.linalg.norm(vec23)

    theta = np.arccos(np.dot(vec21, vec23)/(norm21*norm23))
    if np.isclose(theta-np.pi, 0.0, atol=atol):
        theta = 0.0

    dist13 = np.linalg.norm(vec21-vec23)

    return 2*np.sin(theta) / dist13



def improve_race_line(old_line, inner_border, outer_border):
    # Number of times to iterate each new race line point
    # keep this at 3-8 for best balance of performance and desired result
    XI_ITERATIONS=4

    '''Use gradient descent, inspired by K1999, to find the racing line'''
    # start with the center line
    new_line = copy.deepcopy(old_line)
    ls_inner_border = Polygon(inner_border)
    ls_outer_border = Polygon(outer_border)
    for i in range(0,len(new_line)):
        xi = new_line[i]
        npoints = len(new_line)
        prevprev = (i - 2 + npoints) % npoints
        prev = (i - 1 + npoints) % npoints
        nexxt = (i + 1 + npoints) % npoints
        nexxtnexxt = (i + 2 + npoints) % npoints
        #print("%d: %d %d %d %d %d" % (npoints, prevprev, prev, i, nexxt, nexxtnexxt))
        ci = menger_curvature(new_line[prev], xi, new_line[nexxt])
        c1 = menger_curvature(new_line[prevprev], new_line[prev], xi)
        c2 = menger_curvature(xi, new_line[nexxt], new_line[nexxtnexxt])
        target_ci = (c1 + c2) / 2
        #print("i %d ci %f target_ci %f c1 %f c2 %f" % (i, ci, target_ci, c1, c2))

        # Calculate prospective new track position, start at half-way (curvature zero)
        xi_bound1 = copy.deepcopy(xi)
        xi_bound2 = ((new_line[nexxt][0] + new_line[prev][0]) / 2.0, (new_line[nexxt][1] + new_line[prev][1]) / 2.0)
        p_xi = copy.deepcopy(xi)
        for j in range(0,XI_ITERATIONS):
            p_ci = menger_curvature(new_line[prev], p_xi, new_line[nexxt])
            #print("i: {} iter {} p_ci {} p_xi {} b1 {} b2 {}".format(i,j,p_ci,p_xi,xi_bound1, xi_bound2))
            if np.isclose(p_ci, target_ci):
                break
            if p_ci < target_ci:
                # too flat, shrinking track too much
                xi_bound2 = copy.deepcopy(p_xi)
                new_p_xi = ((xi_bound1[0] + p_xi[0]) / 2.0, (xi_bound1[1] + p_xi[1]) / 2.0)
                if Point(new_p_xi).within(ls_inner_border) or not Point(new_p_xi).within(ls_outer_border):
                    xi_bound1 = copy.deepcopy(new_p_xi)
                else:
                    p_xi = new_p_xi
            else:
                # too curved, flatten it out
                xi_bound1 = copy.deepcopy(p_xi)
                new_p_xi = ((xi_bound2[0] + p_xi[0]) / 2.0, (xi_bound2[1] + p_xi[1]) / 2.0)

                # If iteration pushes the point beyond the border of the track,
                # just abandon the refinement at this point.  As adjacent
                # points are adjusted within the track the point should gradually
                # make its way to a new position.  A better way would be to use
                # a projection of the point on the border as the new bound.  Later.
                if Point(new_p_xi).within(ls_inner_border) or not Point(new_p_xi).within(ls_outer_border):
                    xi_bound2 = copy.deepcopy(new_p_xi)
                else:
                    p_xi = new_p_xi
        new_xi = p_xi
        # New point which has mid-curvature of prev and next points but may be outside of track
        #print((new_line[i], new_xi))
        new_line[i] = new_xi
    return new_line


# track_name = "Karting Altafulla"
track_name = "Karting Vendrell"
track = load_track(f'tracks/{track_name}.kml')
inner_side = np.array(track['track']['inner'])
outer_side = np.array(track['track']['outer'])
base_side = np.array(track['base'])
# add_plot(inner_side)
# add_plot(outer_side)
# add_plot(base_side)

num_samples = 200
resampled_outer = resample_track(outer_side, num_samples)
resampled_inner = resample_track(inner_side, num_samples)
base_line = resample_track(base_side, num_samples)

nearest_inner_points, distances = find_nearest_points(resampled_outer, resampled_inner)

middle_line = (resampled_outer + nearest_inner_points) / 2
inner_border = nearest_inner_points
outer_border = resampled_outer

save_track(f'tracks/{track_name} Sane.kml', track['name'], 
           track['center'], inner_border, outer_border, 
           middle_line, None, base_line)

# Number of times to scan the entire race track to iterate
# 500 will get a good start, 1500 will be closer to optimal result
LINE_ITERATIONS=400

print(len(base_line))
# start along centerline of track
race_line = copy.deepcopy(base_line[:-1])  # Use this for centerline being outer bound
for i in range(LINE_ITERATIONS):
    race_line = improve_race_line(race_line, inner_border, outer_border)
    if i % 20 == 0: print("Iteration %d" % i)

# need to put duplicate point race_line[0] at race_line[-1] to make a closed loops
loop_race_line = np.append(race_line, [race_line[0]], axis=0)


add_plot(nearest_inner_points)
add_plot(resampled_outer)
# add_plot(middle_line)
add_plot(base_line)
add_plot(loop_race_line)
plt.show()


save_track(f'tracks/{track_name} Race Line.kml', track['name'], 
           track['center'], inner_border, outer_border, 
           middle_line, loop_race_line, base_line)

pass

