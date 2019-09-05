import unittest


class TestScriptSSAReconstrution(unittest.TestCase):

    def test_execute(self):
        Aperture(8, 8, radius=4, data=self.test_data)


if __name__ == "__main__":
    unittest.main()
