import unittest
import os
from ..csst_fgs_match.Ecliptic_polyMatch import core_poly_match


class CsstfgsPolyMatch(unittest.TestCase):
    def test_poly_match(self):

        '''
        Aim
        --
        Test 'Ecliptic_polyMatch' model.

        Criteria
        ---
        The input file created.

        Details
        -----
        Use 'core_poly_match' function to output the file.
        If the file doesn't exist, then raise 'AssertionError'.
        '''

        star_cata = './data/input_star_catalog.csv'
        ra_column = 'RAJ2000'
        dec_column = 'DEJ2000'
        mag_column = 'Gmag'
        pointing_file = './data/input_pointing.csv'
        pointing_lon = 'LON'
        pointing_lat = 'LAT'
        cmos_peak_pos = './data/cmos_peak_position_msc.csv'
        output_file = './output_star_catalog_msc.csv'
        numPool = 10
        core_poly_match(star_cata, ra_column, dec_column,
                        mag_column, pointing_file, pointing_lon,
                        pointing_lat, cmos_peak_pos, output_file,
                        numPool)
        self.assertTrue(os.path.exists(output_file))


if __name__ == '__main__':
    unittest.main()
