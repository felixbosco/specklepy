import unittest
from specklepy.core import bootstrap


class TestBootstrap(unittest.TestCase):

    def setUp(self):
        self.nFrames = 800
        self.nDraws = 10

    def test_get_draw_vectors(self):
        sample_draw_vectors = bootstrap.random_draw_vectors(number_draws=self.nDraws, number_frames=self.nFrames)
        # imshow(sample_draw_vectors)
        self.assertEqual(len(sample_draw_vectors[0]), self.nFrames)



if __name__ == "__main__":
    unittest.main()
