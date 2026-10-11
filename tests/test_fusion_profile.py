import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ploidyspec.fusion_profile import junction


class TestJunction(unittest.TestCase):
    def test_boundary_between_two_components(self):
        labels = [3, 3, 3, 4, 3, 3, 4, 4, 4, 4]
        i, a, b, agree = junction(labels)
        self.assertEqual((i, a, b), (6, 3, 4))
        self.assertAlmostEqual(agree, 0.9)

    def test_one_component_has_no_junction(self):
        self.assertIsNone(junction([3, 3, 3, 3, 3]))


if __name__ == "__main__":
    unittest.main()
