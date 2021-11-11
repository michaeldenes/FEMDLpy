import numpy as np

# Define new function
def equipts_on_sphere_box(numpts = int(4 * np.pi * (6378**2))+1, box=None):
    """
    Return a list of equidistributed points on a sphere contained inside a box,
    box must be in [min_lon, max_lon, min_lat, max_lat] format
    """
    # Box = [min_lon, max_lon, min_lat, max_lat]
    # e.g box = [-9,-1,-33,-26]
    if box is None:
        #lats = 180*np.arccos(1-2*keep_indices/numpts)/np.pi-90
        #lons = 180*(keep_indices*np.pi*(1+np.sqrt(5)) % (2*np.pi))/np.pi-180
        return "Missing box"
    else:
        # Get values and shift to[0,360]X[0,180]
        alpha_1_deg = box[0] + 180
        alpha_2_deg = box[1] + 180
        beta_1_deg = box[2] + 90
        beta_2_deg = box[3] + 90

        # Convert to radians
        alpha_1_rad = np.pi*alpha_1_deg/180
        alpha_2_rad = np.pi*alpha_2_deg/180
        beta_1_rad = np.pi*beta_1_deg/180
        beta_2_rad = np.pi*beta_2_deg/180

        # Find indices that fit in the latitude range
        # beta_a = arccos(1-2i/n)
        # beta_b = arccos(1-2j/n)
        # -> i = 0.5*(n-ncos(beta_a))
        # -> j = 0.5*(n-ncos(beta_b))
        i = int(0.5*(numpts-numpts*np.cos(beta_1_rad)))-1
        j = int(0.5*(numpts-numpts*np.cos(beta_2_rad)))+1

        #From indices in latitude range find the indices that fit the longitude range
        indices = np.arange(i,j)

        longitudes = indices*np.pi*(1+np.sqrt(5)) % (2*np.pi)
        keep = np.where(np.logical_and(longitudes>=alpha_1_rad, longitudes<=alpha_2_rad))

        keep_indices = indices[keep[0]]

        #Create longitude and latitude arrays
        lats = 180*np.arccos(1-2*keep_indices/numpts)/np.pi-90
        lons = 180*(keep_indices*np.pi*(1+np.sqrt(5)) % (2*np.pi))/np.pi-180
    return np.asarray([lons,lats])
