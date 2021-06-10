import unittest
import random
import numpy as np
from mock import MagicMock


from model import logic

# TODO: attr_count seems like a really weak part of the design
# But does using NumPy arrays require this sort of design?
ATTR_COUNT = 7 # Number of attributes associated with a Particle

  
class TestGetEventParticles(unittest.TestCase):

    def test_all_active_and_1_subregion_returns_valid_list(self):
        test_length = 10
        num_particles = 3

        # Mock subregion
        mock_subregion = MagicMock()
        mock_subregion.leftBoundary.return_value = 0
        mock_subregion.rightBoundary.return_value = test_length
        mock_sub_list = [mock_subregion]

        # Entrainment
        entrainment_events = 3

        # Level limit. Shouldn't matter so take random
        level_limit = random.randint(0, 10)

        # Make fake model particle array:
        #   1) Flat
        #   2) Random positions
        #   3) All active
        model_particles = np.zeros((num_particles, ATTR_COUNT))
        model_particles[:,3] = np.arange(num_particles)
        model_particles[:,4] = np.ones(num_particles)
        model_particles[:,0] = np.random.randint(test_length, size=num_particles)

        list = logic.get_event_particles(entrainment_events, mock_sub_list, 
                                                    model_particles, level_limit)

        self.assertNotNone(list)
        self.assertCountEqual(list, model_particles[:,3])


if __name__ == '__main__':
    unittest.main()