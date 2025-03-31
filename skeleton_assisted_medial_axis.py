# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 10:35:34 2025

@author: Alexandros Papagiannakis, PhD
"""

import matplotlib.pyplot as plt
import skimage.morphology as morph
import numpy as np

def image_framing(image, frame_width):
    y_array = np.zeros((image.shape[0], frame_width))
    image = np.concatenate((y_array, image), axis=1)
    image = np.concatenate((image, y_array), axis=1)

    x_array = np.zeros((frame_width, image.shape[1]))
    image = np.concatenate((x_array, image), axis=0)
    image = np.concatenate((image, x_array), axis=0)
    return image

def get_intersection_points(num_neighbors):
    return np.nonzero(num_neighbors>4)

def get_merging_points(num_neighbors):
    return np.nonzero(num_neighbors>3)

def get_square(skeleton, x, y):
    return skeleton[y-1:y+2, x-1:x+2]

def skeletonize_mask(cropped_mask):
    return morph.medial_axis(cropped_mask, return_distance=True)

def get_pixel_neighbors(skeleton):
    num_neighbors = np.zeros(skeleton.shape, dtype=int)
    for i in range(0, skeleton.shape[0]):
        for j in range(0, skeleton.shape[1]):
            if skeleton[i,j]==True:
                region = get_square(skeleton, j, i)
                num_neighbors[i,j] = np.sum(region)
    return num_neighbors

def get_position_identity(nodes, x, y):
    if nodes[y][x] == 2:
        square = get_square(nodes, x, y)
        if np.max(square.ravel()) == 3:
            return 'edge'
        elif np.max(square.ravel()) == 4:
            return 'merge'
        elif np.max(square.ravel()) > 4:
            return 'cross'
        
def get_edge_coordinates(nodes):
    
    edge_coords = np.nonzero(nodes==2)
    start_x = []
    start_y = []
    for i in range(edge_coords[0].shape[0]):
        if get_position_identity(nodes, edge_coords[1][i], edge_coords[0][i]) == 'edge':
            start_x.append(edge_coords[1][i])
            start_y.append(edge_coords[0][i])
    return start_x, start_y
    

def get_next_node(coord_list, nodes):
    """
    Parameters
    ----------
    coord_list : list of (x,y) tuples
        A list of the central line edge coordinates that will be populated with recursion
    nodes : 2D numpy array
        A 2D numpy array returned from the get_pixel_neighbors function

    Returns
    -------
    list of (x,y) tuples
        The list of the ordered (x,y) coordinates between the two cell edges

    """
    recurs = 0
    last_coord = coord_list[-1]
    last_x = last_coord[0]
    last_y = last_coord[1]
    
    square = get_square(nodes, last_x, last_y)
    nonzero_coords = np.nonzero(square>=3)
    for cord_i in range(nonzero_coords[0].shape[0]):
        new_x = nonzero_coords[1][cord_i]+last_x-1
        new_y = nonzero_coords[0][cord_i]+last_y-1
        if (new_x, new_y) not in coord_list:
            coord_list.append((new_x, new_y))
            recurs = 1 
    if recurs == 0:
        return coord_list
    elif recurs == 1:
        return get_next_node(coord_list, nodes)
        

def scan_skeleton(skeleton):
    neigh = get_pixel_neighbors(skeleton)
    start_x, start_y = get_edge_coordinates(neigh)
    print(start_x, start_y)
    if len(start_x) == 2 and len(start_y) == 2:
        coord_list = [(start_x[0], start_y[0])]
        coord_list = get_next_node(coord_list, neigh)
    else:
        raise ValueError('Abnormal edges detected')
    coord_list.append((start_x[1], start_y[1]))
    return coord_list


def unwrap_coordinates(coord_list, frame):
    x_coords = []
    y_coords = []
    for crdi in coord_list:
        x_coords.append(crdi[0]-frame)
        y_coords.append(crdi[1]-frame)
    return x_coords, y_coords, list(np.arange(len(x_coords)))


def fit_bivariate_polynomials(i_coords, x_coords, y_coords, degree):
    
    x_fit = np.polyfit(i_coords, x_coords, degree)
    y_fit = np.polyfit(i_coords, y_coords, degree)
    
    x_hat = np.polyval(x_fit, i_coords)
    y_hat = np.polyval(y_fit, i_coords)
    
    return x_hat, y_hat


def skeleton_assisted_bivariate_axis(cropped_mask, degree):
    """
    This function can be used to implement the skeleton assisted central line estimation    
    
    Parameters
    ----------
    cropped_mask : 2D binary mumpy array
        The cropped cell mask where True is the cell pixels and False the background pixels

    Raises
    ------
    ValueError
        When there are cell crossing events

    Returns
    -------
    Pandas DataDrame
        columns: 
            i_coords: the index coorinates of the medial axis used for the bivariate fit
            x_coords: the raw x coordinates of the medial axis linearly interpolated for sub-pixel resolution
            y_coords: the raw y coordinates of the medial axis linearly interpolated for sub-pixel resolution
            x_hat: the bivariate-fitted x coordinates using numpy polyfit and interpoalted linearly
            y_hat: the bivariate-fitted y coordinates using numpy polyfit and interpoalted linearly
            x: the smoothed x coordinates using rol=3 pixels, min-periods=1 pixel and center=True interpolated linearly
            y: the smoothed y coordinates using rol=3 pixels, min-periods=1 pixel and center=True interpolated linearly

    """
    cropped_mask = image_framing(cropped_mask, 1)
    skel, dist = skeletonize_mask(cropped_mask)
    neigh = get_pixel_neighbors(skel)
    
    # plt.imshow(neigh)
    # plt.colorbar()
    # plt.show()
    # plt.imshow(cropped_mask)
    # plt.show()
    print(np.max(neigh))
    if np.max(neigh)<=4:
        coord_list = scan_skeleton(skel)
        x_coords, y_coords, i_coords = unwrap_coordinates(coord_list, 1)
        x_hat, y_hat = fit_bivariate_polynomials(i_coords, x_coords, y_coords, degree)
        
        axis_df=  pd.DataFrame()
        axis_df['i_coords'] = i_coords
        axis_df['x_coords'] = x_coords
        axis_df['y_coords'] = y_coords
        axis_df['x_hat'] = x_hat
        axis_df['y_hat'] = y_hat
        roll_df = axis_df.rolling(3, min_periods=1, center=True).mean()
        axis_df['x'] = roll_df.x_coords
        axis_df['y'] = roll_df.y_coords
        axis_df = axis_df.reindex(np.arange(0,axis_df.shape[0]-1+0.1,0.1)).interpolate()
        
        return axis_df
    
    else:
        raise ValueError('Cell crossing event detected')


