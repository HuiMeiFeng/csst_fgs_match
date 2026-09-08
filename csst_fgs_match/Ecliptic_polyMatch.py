"""
Identifier:     KSC-SJA-csst_fgs_match/Ecliptic_polyMatch.py
Name:           Ecliptic_polyMatch.py
Description:    Create a execution file of CSST FGS matching result
Author:         Huimei Feng
Created:        2024-01-05
Modified-History:
    2024-01-05, Huimei Feng, created
"""

import numpy as np
import pandas as pd
# import time
import os
import multiprocessing
from astropy.coordinates import SkyCoord
from shapely.geometry import Point, Polygon


def radMatch(lon, lat, A, D):
    rads = np.degrees(np.arccos(np.cos(np.radians(lat))
                                * np.cos(np.radians(D))
                                * np.cos(np.radians(lon-A))
                                + np.sin(np.radians(lat))
                                * np.sin(np.radians(D))))
    L = lon[rads < 1.8]
    B = lat[rads < 1.8]
    inx = np.where(rads < 1.8)[0]
    dicts = {'lon': L, 'lat': B, 'index': inx}
    return dicts


def gnoPro(lon, lat, A, D):
    dicts = radMatch(lon, lat, A, D)
    L = dicts['lon']
    B = dicts['lat']
    xi = []
    eta = []
    for i in range(len(L)):
        xi_i = np.degrees((np.cos(np.radians(B[i]))*np.sin(np.radians(L[i]-A)))
                          / (np.sin(np.radians(B[i])) * np.sin(np.radians(D))
                          + np.cos(np.radians(B[i])) * np.cos(np.radians(D))
                          * np.cos(np.radians(L[i]-A))))
        eta_i = np.degrees((np.sin(np.radians(B[i])) * np.cos(np.radians(D))
                           - np.cos(np.radians(B[i])) * np.sin(np.radians(D))
                           * np.cos(np.radians(L[i]-A)))
                           / (np.sin(np.radians(B[i])) * np.sin(np.radians(D))
                           + np.cos(np.radians(B[i])) * np.cos(np.radians(D))
                           * np.cos(np.radians(L[i]-A))))
        xi.append(xi_i)
        eta.append(eta_i)
    dicts_1 = {'xi': xi, 'eta': eta}
    dicts.update(dicts_1)
    return dicts


def isInPoly(xi, eta, polyg):
    points_xi_in_poly = []
    points_eta_in_poly = []
    for point in zip(xi, eta):
        points = Point(np.array(point))
        poly = Polygon(polyg)
        if poly.contains(points):
            points_xi_in_poly.append(list(points.coords)[0][0])
            points_eta_in_poly.append(list(points.coords)[0][1])
    return points_xi_in_poly, points_eta_in_poly


# match and write file
def matchWriteFile(ra, dec, mag, lon, lat, A, D, polyg, output_file):
    dicts = gnoPro(lon, lat, A, D)
    xi = dicts['xi']
    eta = dicts['eta']
    inx_rad = dicts['index']

    ass = []
    bss = []
    guide_ids = []
    for i in range(0, polyg.shape[0], 1):
        poly = polyg[i]
        ai, bi = isInPoly(xi, eta, poly)
        ass.append(ai)
        bss.append(bi)
        guide_ids.append([i+1])

    match_len = len([item for t in ass for item in t])

    # write information of matched sources that match with gaia catalogue
    # write information of objects without guiding
    if match_len == 0:
        df_lst_zip_1 = pd.DataFrame(zip([-9999], [-9999], [-9999], [-9999]))
        df_lst_zip_1.to_csv(output_file, mode='a', index=False,
                            header=False)
    else:
        ass_1 = [item for t in ass for item in t]

        xi = np.array(xi)
        inx_lst = []
        for i in ass_1:
            indx = np.where(xi == i)[0][0]
            inx_f = inx_rad[indx]
            inx_lst.append(inx_f)

        RA = []
        DEC = []
        MAG = []
        for i in inx_lst:
            RA.append(ra[i])
            DEC.append(dec[i])
            MAG.append(mag[i])

        if polyg.shape[0] == 8:
            guide_idss = [guide_ids[0]*len(ass[0]),
                          guide_ids[1]*len(ass[1]),
                          guide_ids[2]*len(ass[2]),
                          guide_ids[3]*len(ass[3]),
                          guide_ids[4]*len(ass[4]),
                          guide_ids[5]*len(ass[5]),
                          guide_ids[6]*len(ass[6]),
                          guide_ids[7]*len(ass[7]),]
        else:
            guide_idss = [guide_ids[0]*len(ass[0]),
                          guide_ids[1]*len(ass[1]),
                          guide_ids[2]*len(ass[2]),
                          guide_ids[3]*len(ass[3]),]
        guide_id = [item for t in guide_idss for item in t]
        df_lst_zip_2 = pd.DataFrame(zip(RA, DEC, MAG, guide_id))
        df_lst_zip_2.to_csv(output_file, mode='a', index=False,
                            header=False)


def core_poly_match(star_cata, ra_column, dec_column, mag_column,
                    pointing_file, pointing_lon, pointing_lat,
                    cmos_peak_pos, output_file, numPool):

    """
    Calculation of available stars in each cmos with one pointing in rectangle.

    The functions of matching available stars are still developing. Currently,
    we require ten parameters to execute the pipeline.

    Parameters
    ----------
    star_cata : str (Required)
        The input star catalogue file name with full directory.

    ra_column : str (Required)
        The name of right ascesion (R.A. in degree unit) in input star catalog.

    dec_column: str (Required)
        The name of declination (DEC in degree unit) in input star catalogue.

    mag_column: str (Required)
        The name of magnitude in input star catalogue.

    pointing_file: str (Required)
        The file path of input pointing catalogue.

    pointing_lon: str (Required)
        The name of ecliptic longtitude in pointing catalogue.

    pointing_lat: str (Required)
        The name of ecliptic latitude in pointing catalogue.

    cmos_pos: str (Required)
        The file path of 8 guiding cmos position file which in degree unit.

    output_file: str (Required)
        The output file path which will contain all matching results.

    numPool: int (Required)
        The number of pool is going to be used.

    Returns
    -------
    output_file
        Result containing all matching stars.

    Examples
    ----------
    >>> core_poly_match('./UCAC5_2289pointing_1deg_20221107.csv',
    'RAJ2000', 'DEJ2000', 'Gmag', './pointing_2289_20221107.csv', 'LON', 'LAT',
    './cmos_peak_position.csv', './output_star_catalog.csv', 20)
    >>> get the results file.
    """

    tab = pd.read_csv(star_cata)
    ra = np.array(tab[ra_column])
    dec = np.array(tab[dec_column])
    mag = np.array(tab[mag_column])
    c_equ = SkyCoord(ra, dec, frame='icrs', unit="deg")
    c_ecl = c_equ.barycentricmeanecliptic
    lon = np.array(c_ecl.lon)  # ecliptic longitude
    lat = np.array(c_ecl.lat)  # ecliptic latitude

    # read pointing
    data = pd.read_csv(pointing_file)
    A = np.array(data[pointing_lon])
    D = np.array(data[pointing_lat])

    # read position of CMOSs
    df = pd.read_csv(cmos_peak_pos)
    if df.shape[1] == 8:
        polyg = np.array([np.array(df[['cmos1_x', 'cmos1_y']]).tolist(),
                          np.array(df[['cmos2_x', 'cmos2_y']]).tolist(),
                          np.array(df[['cmos3_x', 'cmos3_y']]).tolist(),
                          np.array(df[['cmos4_x', 'cmos4_y']]).tolist(),])
    else:
        polyg = np.array([np.array(df[['cmos1_x', 'cmos1_y']]).tolist(),
                          np.array(df[['cmos2_x', 'cmos2_y']]).tolist(),
                          np.array(df[['cmos3_x', 'cmos3_y']]).tolist(),
                          np.array(df[['cmos4_x', 'cmos4_y']]).tolist(),
                          np.array(df[['cmos5_x', 'cmos5_y']]).tolist(),
                          np.array(df[['cmos6_x', 'cmos6_y']]).tolist(),
                          np.array(df[['cmos7_x', 'cmos7_y']]).tolist(),
                          np.array(df[['cmos8_x', 'cmos8_y']]).tolist(),])

    pool = multiprocessing.Pool(numPool)
    for i in range(A.size):
        pool.apply_async(func=matchWriteFile,
                         args=(ra, dec, mag, lon, lat, A[i], D[i],
                               polyg, output_file,))
    pool.close()
    pool.join()

    # write header of output file
    if os.path.exists(output_file):
        ef = pd.read_csv(output_file, header=None,
                         names=['guide_ra', 'guide_dec',
                                'phot_g_mean_mag', 'cmos_id'])
        ef.to_csv(output_file, index=False)


# if __name__ == '__main__':
#    print('Parent process %s.' % os.getpid())
#    start = time.time()

#    core_poly_match('./UCAC5_2289pointing_1deg_20221107.csv', 'RAJ2000',
#                    'DEJ2000', 'Gmag', './pointing_2289_20221107.csv',
#                    'LON', 'LAT', './cmos_peak_position_f4.csv',
#                    './output_star_catalog_f4.csv', 20)

#    print('all escape:', time.time() - start)
