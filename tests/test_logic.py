import unittest
import random
import numpy as np
from mock import MagicMock


from model import logic

# TODO: attr_count seems like a really weak part of the design
# But does using NumPy arrays require this sort of design?
ATTR_COUNT = 7 # Number of attributes associated with a Particle

  
class TestGetEventParticlesWithOneSubregion(unittest.TestCase):
    def setUp(self):
        self.test_length = 10
        self.num_particles = 3

        mock_subregion = MagicMock()
        mock_subregion.leftBoundary.return_value = 0
        mock_subregion.rightBoundary.return_value = self.test_length
        mock_subregion.getName.return_value = 'Mock_Subregion'
        self.mock_sub_list = [mock_subregion]

        self.entrainment_events = 3

        self.level_limit = random.randint(0, random.randint(2, 10))


    def test_all_active_returns_valid_list(self):
        
        model_particles = np.zeros((self.num_particles, ATTR_COUNT))
        model_particles[:,3] = np.arange(self.num_particles) # unique ids
        model_particles[:,4] = np.ones(self.num_particles) # all active
        model_particles[:,0] = np.random.randint(
                                                self.test_length, 
                                                size=self.num_particles ) # random placement

        list = logic.get_event_particles(
                                        self.entrainment_events, 
                                        self.mock_sub_list, 
                                        model_particles, 
                                        self.level_limit )
        self.assertCountEqual(list, model_particles[:,3])

        # Height dependancy should not effect list results here
        hp_list = list = logic.get_event_particles(
                                        self.entrainment_events, 
                                        self.mock_sub_list, 
                                        model_particles, 
                                        self.level_limit,
                                        height_dependant=True )
        self.assertCountEqual(hp_list, model_particles[:,3])
        self.assertCountEqual(hp_list, list)
    
    # TODO: Mock logging object and assert warning is logged
    def test_not_all_active_returns_valid_list(self):

        mp_one_inactive = np.zeros((self.num_particles, ATTR_COUNT))
        mp_one_inactive[:,3] = np.arange(self.num_particles) 
        mp_one_inactive[0][4] = 1
        mp_one_inactive[1][4] = 1
        mp_one_inactive[:,0] = np.random.randint(self.test_length, size=self.num_particles)

        list = logic.get_event_particles(
                                        self.entrainment_events, 
                                        self.mock_sub_list, 
                                        mp_one_inactive, 
                                        self.level_limit )
        self.assertEqual(len(list), 2)
        active_list = mp_one_inactive[mp_one_inactive[:,4] != 0]
        self.assertCountEqual(list, active_list[:,3])


if __name__ == '__main__':
    unittest.main()