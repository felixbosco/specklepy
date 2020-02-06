import unittest
from specklepy.core import bootstrap
from specklepy.utils.plot import imshow



class TestBootstrap(unittest.TestCase):

    def setUp(self):
        self.nFrames = 800
        self.nDraws = 10

    def test_get_draw_vectors(self):
        sample_draw_vectors = bootstrap.get_sample_draw_vectors(nDraws=self.nDraws, nFrames=self.nFrames)
        # imshow(sample_draw_vectors)
        self.assertEqual(len(sample_draw_vectors[0]), self.nFrames)



if __name__ == "__main__":
    unittest.main()
